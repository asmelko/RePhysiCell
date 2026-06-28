#ifndef __PhysiCell_mechanics_agent_h__
#define __PhysiCell_mechanics_agent_h__

#include <vector>

#include "PhysiCell_mechanics_agent_interface.h"
#include "PhysiCell_mechanics_agent_PIMPL.h"
#include "PhysiCell_mechanics_functions.h"

#include "../core/PhysiCell_phenotype.h"

namespace PhysiCell {

class Cell;
class mechanics_environment;

/**
 * @brief Concrete mechanics agent that owns its mechanics data directly.
 *
 * Data members mirror the per-agent fields of mech_agent_data and its
 * sub-structs (base_membrane_data, base_motility_data, base_potential_data),
 * but held as single-agent values rather than SoA arrays.
 *
 * Every Mechanics_Agent_Interface getter is implemented by returning a
 * reference or pointer into the corresponding member.
 */
class Mechanics_Agent : public Mechanics_Agent_Interface
{
public:
	// ========================================================================
	// mech_agent_data — direct fields
	// ========================================================================

	std::vector<double> velocity;
	bool                is_movable            = true;
	std::vector<Mechanics_Agent*>    neighbors;


	std::vector<int>    neighbors2, springs;

	// ========================================================================
	// mechanics_data (membrane + potential data)
	// ========================================================================

	Mechanics_Data mechanics_data;

	// ========================================================================
	// motility_data
	// ========================================================================

	Motility_Data motility_data;

	Radius_Data radius_data;

	// ========================================================================
	// Remaining direct fields not captured by Mechanics_Data / Motility_Data
	// ========================================================================

	// std::function<void(double*)> update_migration_bias_direction;

	double           simple_pressure          = 0.0;
	std::vector<double> previous_velocity;

	// ========================================================================
	// Mechanics_Agent_Interface implementations
	// ========================================================================
	
	int type;
	int get_type() const override;
	void set_type(int new_type) override;

	// ---- mech_agent_data — direct fields -----------------------------------

	double* get_velocity() override
	{ return velocity.data(); }

	double& get_radius() override
	{ return radius_data.radius; }

	bool& get_is_movable() override
	{ return is_movable; }

	std::vector<int>& get_neighbors() override
	{ return neighbors2; }

	// ---- base_membrane_data ------------------------------------------------

	double& get_cell_BM_adhesion_strength() override
	{ return mechanics_data.cell_BM_adhesion_strength; }

	double& get_cell_BM_repulsion_strength() override
	{ return mechanics_data.cell_BM_repulsion_strength; }

	// ---- base_motility_data ------------------------------------------------

	bool& get_is_motile() override
	{ return motility_data.is_motile; }

	double& get_persistence_time() override
	{ return motility_data.persistence_time; }

	double& get_migration_speed() override
	{ return motility_data.migration_speed; }

	double* get_migration_bias_direction() override
	{ return motility_data.migration_bias_direction.data(); }

	double& get_migration_bias() override
	{ return motility_data.migration_bias; }

	double* get_motility_vector() override
	{ return motility_data.motility_vector.data(); }

	bool& get_restrict_to_2d() override
	{ return motility_data.restrict_to_2D; }

	int& get_chemotaxis_index() override
	{ return motility_data.chemotaxis_index; }

	int& get_chemotaxis_direction() override
	{ return motility_data.chemotaxis_direction; }

	double* get_chemotactic_sensitivities() override
	{ return motility_data.chemotactic_sensitivities.data(); }

	// std::function<void(double*)>& get_update_migration_bias_direction() override
	// { return update_migration_bias_direction; }

	// ---- base_potential_data -----------------------------------------------

	double& get_cell_cell_adhesion_strength() override
	{ return mechanics_data.cell_cell_adhesion_strength; }

	double& get_cell_cell_repulsion_strength() override
	{ return mechanics_data.cell_cell_repulsion_strength; }

	double* get_cell_adhesion_affinities() override
	{ return mechanics_data.cell_adhesion_affinities.data(); }

	double& get_relative_maximum_adhesion_distance() override
	{ return mechanics_data.relative_maximum_adhesion_distance; }

	int& get_maximum_number_of_attachments() override
	{ return mechanics_data.maximum_number_of_attachments; }

	double& get_attachment_elastic_constant() override
	{ return mechanics_data.attachment_elastic_constant; }

	double& get_attachment_rate() override
	{ return mechanics_data.attachment_rate; }

	double& get_detachment_rate() override
	{ return mechanics_data.detachment_rate; }

	double& get_simple_pressure() override
	{ return simple_pressure; }

	double* get_previous_velocity() override
	{ return previous_velocity.data(); }

	std::vector<int>& get_springs() override
	{ return springs; }
	
	
	
	double& get_relative_maximum_attachment_distance() override
	{ return mechanics_data.relative_maximum_attachment_distance; }

	double& get_relative_detachment_distance() override
	{ return mechanics_data.relative_detachment_distance; }

	double& get_maximum_attachment_rate() override
	{ return mechanics_data.maximum_attachment_rate; }

	Mechanics_Agent(Cell* pCell);

	// Non-owning pointer to the canonical position storage.
	// Initialised from the wrapped Basic_Agent; rebound by Cell to point at
	// the Cell's own Position_Entity subobject.
	BioFVM::Position_Entity* pos_entity = nullptr;

	void bind_position_entity(BioFVM::Position_Entity* pe) override;

	Position_Entity*  get_position_entity() noexcept override;

	// Owned fallback used when no external Position_Entity is supplied
	// (e.g. standalone Basic_Agent not embedded in a Cell).
	Position_Entity default_position;


	Mechanics_Agent_PIMPL* pOwner = nullptr;

	Mechanics_Functions functions;


	std::vector<Mechanics_Agent*> attached_cells; 
	std::vector<std::pair<Mechanics_Agent*, bool>> spring_attachments; 
	bool is_out_of_domain = false;
	int current_mechanics_voxel_index = -1;
	int updated_current_mechanics_voxel_index = 0; // keeps the updated voxel index for later adjusting of current voxel index

	void update_motility_vector( double dt_ );
	void add_potentials(Mechanics_Agent*);       // Add repulsive and adhesive forces.
	void set_previous_velocity(double xV, double yV, double zV);
	int get_current_mechanics_voxel_index() override;
	bool assign_position(const std::vector<double>& new_position) override;
	bool assign_position(double, double, double) override;

		// mechanics 
	void update_position( double dt ) override; //
	std::vector<double> displacement = {0.0, 0.0, 0.0}; // this should be moved to state, or made private  

	void update_voxel_in_container(void) override;

	void attach_cell( Mechanics_Agent* pAddMe ); // done 
	void detach_cell( Mechanics_Agent* pRemoveMe ); // done 

	void remove_self_from_all_neighbors( void ) override; 
	void remove_all_attached_cells( void ) override; // done 

	void attach_cell_as_spring( Mechanics_Agent* pAddMe, bool attacking_spring ); // done 
	void detach_cell_as_spring( Mechanics_Agent* pRemoveMe ); // done 
	void remove_all_spring_attachments( void ) override; // done 

	
	bool& get_is_out_of_domain() override
	{ return is_out_of_domain; }

	
	int number_of_attached_cells( void ); 
};



} // namespace PhysiCell

#endif // __PhysiCell_mechanics_agent_h__
