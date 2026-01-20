#include <biofvm/agent_container.h>

#include "basic_agent_adapter.h"
#include "microenvironment_adapter.h"
#include "../../BioFVM/BioFVM_microenvironment_interface.h"
#include <stdexcept>

namespace BioFVM {

// ============================================================================
// Constructor
// ============================================================================

int physicore_basic_agent_adapter::ID_counter = 0;

physicore_basic_agent_adapter::physicore_basic_agent_adapter()
	: index(0)
	, microenvironment(dynamic_cast<physicore_microenvironment_adapter*>(get_microenvironment_i()))
	, ID(static_cast<int>(ID_counter++))
	, type(0)
	, cached_voxel_index(-1)
	, position_cache(3, 0.0)
	, velocity(3, 0.0)
	, previous_velocity(3, 0.0)
{
	auto physicore_microenvironment = microenvironment->get_physicore_microenvironment();

	agent = dynamic_cast<physicore::biofvm::agent_container*>(physicore_microenvironment->agents.get())->create();

	agent->volume() = 1.0; // Default volume
	agent->is_active() = 1; // Active by default
}

physicore_basic_agent_adapter::~physicore_basic_agent_adapter()
{
	auto physicore_microenvironment = microenvironment->get_physicore_microenvironment();

	dynamic_cast<physicore::biofvm::agent_container*>(physicore_microenvironment->agents.get())->remove_agent(agent);
}

// ============================================================================
// Volume methods
// ============================================================================

double& physicore_basic_agent_adapter::get_total_volume()
{
	return agent->volume();
}

void physicore_basic_agent_adapter::set_total_volume(double volume)
{
	agent->volume() = volume;
}

void physicore_basic_agent_adapter::update_voxel_index()
{
	if (microenvironment) {
        if(!microenvironment->get_mesh().is_position_valid(
            get_position()[0], get_position()[1], get_position()[2]))
        {	
            cached_voxel_index = -1;
            agent->is_active() = 0;
            return;
        }
		cached_voxel_index = microenvironment->nearest_voxel_index(get_position());
	}
}

// ============================================================================
// Microenvironment registration and access
// ============================================================================

void physicore_basic_agent_adapter::register_microenvironment(Microenvironment_Interface* microenv)
{
	if (microenv) {
		// Resize cached vectors to match substrate count
		unsigned int num_substrates = microenv->number_of_densities();
		// Note: substrate arrays in physicore SoA already have correct size
		// No need to resize them, just verify
		assert(agent->secretion_rates().size() == num_substrates);
	}
}

Microenvironment_Interface* physicore_basic_agent_adapter::get_microenvironment_interface()
{
	return microenvironment;
}

// ============================================================================
// ID and type accessors
// ============================================================================

int physicore_basic_agent_adapter::get_ID() const
{
	return ID;
}

void physicore_basic_agent_adapter::set_ID(int new_ID)
{
	ID = new_ID;
}

int physicore_basic_agent_adapter::get_index() const
{
	return index;
}

void physicore_basic_agent_adapter::set_index(int new_index)
{
	// Note: Underlying physicore agent has its own index in SoA
	// This index is just for BioFVM interface purposes
	index = new_index;
}

int physicore_basic_agent_adapter::get_type() const
{
	return type;
}

void physicore_basic_agent_adapter::set_type(int new_type)
{
	type = new_type;
}

// ============================================================================
// Position methods
// ============================================================================

bool physicore_basic_agent_adapter::assign_position(double x, double y, double z)
{
	if(!microenvironment->get_mesh().is_position_valid(x,y,z))
	{	
		// std::cout<<"Error: the new position for agent "<< ID << " is invalid: "<<x<<","<<y<<","<<"z"<<std::endl;
		return false;
	}

	// Write directly to physicore SoA via public interface
	auto pos = agent->position();
	pos[0] = x;
	pos[1] = y;
	if (pos.size() > 2) {
		pos[2] = z;
	}
	position_cache[2] = z; // Update cache
	
	return true;
}

bool physicore_basic_agent_adapter::assign_position(std::vector<double> new_position)
{
	return assign_position(new_position[0], new_position[1], new_position[2]);
}

double* physicore_basic_agent_adapter::get_position_internal()
{
	return agent->position().data();
}

const std::vector<double>& physicore_basic_agent_adapter::get_position() const
{
	// Lazy sync: load from physicore on first access or after dirty flag set
	for (size_t i = 0; i < agent->position().size(); ++i) {
		position_cache[i] = agent->position()[i];
	}
	
	return position_cache;
}

void physicore_basic_agent_adapter::update_position(double dt)
{

}

// ============================================================================
// Velocity methods
// ============================================================================

std::vector<double>& physicore_basic_agent_adapter::get_velocity()
{
	return velocity;
}

const std::vector<double>& physicore_basic_agent_adapter::get_velocity() const
{
	return velocity;
}

std::vector<double>& physicore_basic_agent_adapter::get_previous_velocity()
{
	return previous_velocity;
}

const std::vector<double>& physicore_basic_agent_adapter::get_previous_velocity() const
{
	return previous_velocity;
}

// ============================================================================
// Activity status
// ============================================================================

bool physicore_basic_agent_adapter::get_is_active() const
{
	return agent->is_active() == 1;
}

void physicore_basic_agent_adapter::set_is_active(bool active)
{
	agent->is_active() = active ? 1 : 0;
}

// ============================================================================
// Substrate interaction rates (direct pointer access to SoA)
// ============================================================================

double* physicore_basic_agent_adapter::get_secretion_rates()
{
	return agent->secretion_rates().data();
}

const double* physicore_basic_agent_adapter::get_secretion_rates() const
{
	return agent->secretion_rates().data();
}

double* physicore_basic_agent_adapter::get_saturation_densities()
{
	return agent->saturation_densities().data();
}

const double* physicore_basic_agent_adapter::get_saturation_densities() const
{
	return agent->saturation_densities().data();
}

double* physicore_basic_agent_adapter::get_uptake_rates()
{
	return agent->uptake_rates().data();
}

const double* physicore_basic_agent_adapter::get_uptake_rates() const
{
	return agent->uptake_rates().data();
}

double* physicore_basic_agent_adapter::get_net_export_rates()
{
	return agent->net_export_rates().data();
}

const double* physicore_basic_agent_adapter::get_net_export_rates() const
{
	return agent->net_export_rates().data();
}

// ============================================================================
// Internalized substrates
// ============================================================================

double* physicore_basic_agent_adapter::get_internalized_total_substrates()
{
	return agent->internalized_substrates().data();
}

const double* physicore_basic_agent_adapter::get_internalized_total_substrates() const
{
	return agent->internalized_substrates().data();
}

double* physicore_basic_agent_adapter::get_fraction_released_at_death()
{
	return agent->fraction_released_at_death().data();
}

const double* physicore_basic_agent_adapter::get_fraction_released_at_death() const
{
	return agent->fraction_released_at_death().data();
}

double* physicore_basic_agent_adapter::get_fraction_transferred_when_ingested()
{
	return agent->fraction_transferred_when_ingested().data();
}

const double* physicore_basic_agent_adapter::get_fraction_transferred_when_ingested() const
{
	return agent->fraction_transferred_when_ingested().data();
}

void physicore_basic_agent_adapter::release_internalized_substrates()
{
	for (int s = 0; s < static_cast<int>(agent->internalized_substrates().size()); ++s) {
		agent->internalized_substrates()[s] /= microenvironment->get_physicore_microenvironment()->mesh.voxel_volume(); // Convert to density
		agent->internalized_substrates()[s] *= agent->fraction_released_at_death()[s]; // Apply fraction to release

		microenvironment->nearest_density_vector(cached_voxel_index)[s] += agent->internalized_substrates()[s]; // Release to voxel
		agent->internalized_substrates()[s] = 0.0; //
	}
}

void physicore_basic_agent_adapter::set_internal_uptake_constants(double dt)
{
	// This is a no-op for physicore since uptake constants are directly used
	// in the centralized solver and do not depend on agent volume.
}

// ============================================================================
// Secretion and uptake simulation (NOT SUPPORTED)
// ============================================================================

void physicore_basic_agent_adapter::simulate_secretion_and_uptake(double dt)
{
	throw std::runtime_error(
		"simulate_secretion_and_uptake() not supported by physicore adapter. "
		"Physicore uses a centralized solver architecture for substrate updates. "
		"Per-agent secretion/uptake updates are incompatible with this design. "
		"Use physicore's solver methods instead.");
}

// ============================================================================
// Voxel access
// ============================================================================

int physicore_basic_agent_adapter::get_current_voxel_index()
{
	return cached_voxel_index;
}

// ============================================================================
// Density and gradient access
// ============================================================================

double* physicore_basic_agent_adapter::nearest_density_vector()
{
	if (!microenvironment) {
		throw std::runtime_error("Microenvironment not registered");
	}
	
	int voxel_idx = get_current_voxel_index();
	return microenvironment->nearest_density_vector(voxel_idx);
}

std::vector<double>& physicore_basic_agent_adapter::nearest_gradient(int substrate_index)
{
	int voxel_idx = get_current_voxel_index();
	auto& gradient = microenvironment->gradient_vector(voxel_idx);
	return gradient[substrate_index];
}

std::vector<std::vector<double>>& physicore_basic_agent_adapter::nearest_gradient_vector()
{
	if (!microenvironment) {
		throw std::runtime_error("Microenvironment not registered");
	}
	
	int voxel_idx = get_current_voxel_index();
	return microenvironment->gradient_vector(voxel_idx);
}

} // namespace BioFVM
