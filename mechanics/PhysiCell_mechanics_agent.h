#ifndef __PhysiCell_mechanics_agent_h__
#define __PhysiCell_mechanics_agent_h__

#include <vector>

#include "PhysiCell_mechanics_agent_interface.h"
#include "../core_interface/PhysiCell_core_entity_PIMPL.h"
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

	std::vector<double> velocity_;
	bool is_movable_ = true;
	std::vector<CoreEntityPIMPL*> neighbors_;
	std::vector<std::pair<CoreEntityPIMPL*, bool>> spring_attachments; 


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

	double simple_pressure_ = 0.0;
	std::vector<double> previous_velocity_;

	// ========================================================================
	// Mechanics_Agent_Interface implementations
	// ========================================================================
	
	int type_;
	int& type() override;

	// ---- mech_agent_data — direct fields -----------------------------------

	double* velocity() override
	{ return velocity_.data(); }

	double& radius() override
	{ return radius_data.radius; }

	bool& is_movable() override
	{ return is_movable_; }

	std::vector<CoreEntityPIMPL*>& neighbors() override
	{ return neighbors_; }

	// ---- base_membrane_data ------------------------------------------------

	double& cell_BM_adhesion_strength() override
	{ return mechanics_data.cell_BM_adhesion_strength; }

	double& cell_BM_repulsion_strength() override
	{ return mechanics_data.cell_BM_repulsion_strength; }

	// ---- base_motility_data ------------------------------------------------

	bool& is_motile() override
	{ return motility_data.is_motile; }

	double& persistence_time() override
	{ return motility_data.persistence_time; }

	double& migration_speed() override
	{ return motility_data.migration_speed; }

	double* migration_bias_direction() override
	{ return motility_data.migration_bias_direction.data(); }

	double& migration_bias() override
	{ return motility_data.migration_bias; }

	double* motility_vector() override
	{ return motility_data.motility_vector.data(); }

	bool& restrict_to_2d() override
	{ return motility_data.restrict_to_2D; }

	int& chemotaxis_index() override
	{ return motility_data.chemotaxis_index; }

	int& chemotaxis_direction() override
	{ return motility_data.chemotaxis_direction; }

	double* chemotactic_sensitivities() override
	{ return motility_data.chemotactic_sensitivities.data(); }

	// std::function<void(double*)>& get_update_migration_bias_direction() override
	// { return update_migration_bias_direction; }

	// ---- base_potential_data -----------------------------------------------

	double& cell_cell_adhesion_strength() override
	{ return mechanics_data.cell_cell_adhesion_strength; }

	double& cell_cell_repulsion_strength() override
	{ return mechanics_data.cell_cell_repulsion_strength; }

	double* cell_adhesion_affinities() override
	{ return mechanics_data.cell_adhesion_affinities.data(); }

	double& relative_maximum_adhesion_distance() override
	{ return mechanics_data.relative_maximum_adhesion_distance; }

	int& maximum_number_of_attachments() override
	{ return mechanics_data.maximum_number_of_attachments; }

	double& attachment_elastic_constant() override
	{ return mechanics_data.attachment_elastic_constant; }

	double& attachment_rate() override
	{ return mechanics_data.attachment_rate; }

	double& detachment_rate() override
	{ return mechanics_data.detachment_rate; }

	double& simple_pressure() override
	{ return simple_pressure_; }

	double* previous_velocity() override
	{ return previous_velocity_.data(); }

	std::vector<std::pair<CoreEntityPIMPL*, bool>>& springs() override
	{ return spring_attachments; }
	
	
	
	double& relative_maximum_attachment_distance() override
	{ return mechanics_data.relative_maximum_attachment_distance; }

	double& relative_detachment_distance() override
	{ return mechanics_data.relative_detachment_distance; }

	double& maximum_attachment_rate() override
	{ return mechanics_data.maximum_attachment_rate; }

	Mechanics_Agent(Cell* pCell, std::vector<double>& position);

private:

public:
	// The vector must outlive the agent and retain at least three coordinates.
	std::vector<double>& position;


	bool uses_position(const std::vector<double>& position) const
	{ return &position == &position; }

	CoreEntityPIMPL* pOwner = nullptr;

	Mechanics_Functions functions;


	bool is_out_of_domain_ = false;
	int current_mechanics_voxel_index_ = -1;
	int updated_current_mechanics_voxel_index = 0; // keeps the updated voxel index for later adjusting of current voxel index

	void update_motility_vector( double dt_ );
	void add_potentials(Mechanics_Agent*);       // Add repulsive and adhesive forces.
	int current_mechanics_voxel_index() override;
	bool assign_position(const std::vector<double>& new_position) override;
	bool assign_position(double, double, double) override;

		// mechanics 
	void update_position( double dt ); //
	std::vector<double> displacement = {0.0, 0.0, 0.0}; // this should be moved to state, or made private  

	void update_voxel_in_container(void) override;

	void remove_self_from_all_neighbors( void ) override; 

	void attach_cell_as_spring( Mechanics_Agent* pAddMe, bool attacking_spring ); // done 
	void detach_cell_as_spring( Mechanics_Agent* pRemoveMe ); // done 
	void remove_all_spring_attachments( void ) override; // done 

	
	bool& is_out_of_domain() override
	{ return is_out_of_domain_; }
};



} // namespace PhysiCell

#endif // __PhysiCell_mechanics_agent_h__
