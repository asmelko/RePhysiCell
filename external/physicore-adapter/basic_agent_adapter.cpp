/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2025, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

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

Microenvironment_Interface* physicore_basic_agent_adapter::get_microenvironment()
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

void physicore_basic_agent_adapter::compute_gradient_at_voxel(
	int voxel_index, int substrate_index, std::vector<double>& gradient) const 
{
	gradient.resize(3, 0.0);
	
	auto indices = microenvironment->me->mesh.voxel_position(agent->position());
	int i = indices[0];
	int j = indices[1];
	int k = indices[2];
	
	const auto& grid = microenvironment->me->mesh.grid_shape;
	const auto& voxel_shape = microenvironment->me->mesh.voxel_shape;
	
	// Central finite difference for interior voxels, one-sided for boundaries
	
	// X-direction gradient
	if (i > 0 && i < static_cast<int>(grid[0]) - 1) {
		double rho_plus = microenvironment->me->solver->get_substrate_density(substrate_index, i + 1, j, k);
		double rho_minus = microenvironment->me->solver->get_substrate_density(substrate_index, i - 1, j, k);
		gradient[0] = (rho_plus - rho_minus) / (2.0 * voxel_shape[0]);
	} else if (i == 0 && grid[0] > 1) {
		double rho_0 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_1 = microenvironment->me->solver->get_substrate_density(substrate_index, i + 1, j, k);
		gradient[0] = (rho_1 - rho_0) / voxel_shape[0];
	} else if (i == static_cast<int>(grid[0]) - 1 && grid[0] > 1) {
		double rho_0 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_m1 = microenvironment->me->solver->get_substrate_density(substrate_index, i - 1, j, k);
		gradient[0] = (rho_0 - rho_m1) / voxel_shape[0];
	}
	
	// Y-direction gradient
	if (j > 0 && j < static_cast<int>(grid[1]) - 1) {
		double rho_plus = microenvironment->me->solver->get_substrate_density(substrate_index, i, j + 1, k);
		double rho_minus = microenvironment->me->solver->get_substrate_density(substrate_index, i, j - 1, k);
		gradient[1] = (rho_plus - rho_minus) / (2.0 * voxel_shape[1]);
	} else if (j == 0 && grid[1] > 1) {
		double rho_0 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_1 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j + 1, k);
		gradient[1] = (rho_1 - rho_0) / voxel_shape[1];
	} else if (j == static_cast<int>(grid[1]) - 1 && grid[1] > 1) {
		double rho_0 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_m1 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j - 1, k);
		gradient[1] = (rho_0 - rho_m1) / voxel_shape[1];
	}
	
	// Z-direction gradient
	if (k > 0 && k < static_cast<int>(grid[2]) - 1) {
		double rho_plus = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k + 1);
		double rho_minus = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k - 1);
		gradient[2] = (rho_plus - rho_minus) / (2.0 * voxel_shape[2]);
	} else if (k == 0 && grid[2] > 1) {
		double rho_0 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_1 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k + 1);
		gradient[2] = (rho_1 - rho_0) / voxel_shape[2];
	} else if (k == static_cast<int>(grid[2]) - 1 && grid[2] > 1) {
		double rho_0 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_m1 = microenvironment->me->solver->get_substrate_density(substrate_index, i, j, k - 1);
		gradient[2] = (rho_0 - rho_m1) / voxel_shape[2];
	}
}

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

	gradient_cache.resize(
		microenvironment->number_of_densities(),
		std::vector<double>(3, 0.0));
	
	compute_gradient_at_voxel(voxel_idx, substrate_index, gradient_cache[substrate_index]);

	return gradient_cache[substrate_index];
}

std::vector<std::vector<double>>& physicore_basic_agent_adapter::nearest_gradient_vector()
{
	if (!microenvironment) {
		throw std::runtime_error("Microenvironment not registered");
	}
	
	int voxel_idx = get_current_voxel_index();
	return microenvironment->gradient_vector(voxel_idx);
}

std::vector<Basic_Agent_Interface*> all_basic_agents;

std::vector<Basic_Agent_Interface*>* Agent_Container_Interface::get_all_basic_agents()
{
    return &all_basic_agents;
}

} // namespace BioFVM
