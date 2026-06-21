#ifndef __PhysiCell_mechanics_agent_interface_h__
#define __PhysiCell_mechanics_agent_interface_h__

#include <vector>

namespace PhysiCell {

/**
 * @brief Abstract interface for agents that participate in mechanical interactions.
 *
 * Each pure-virtual method exposes a mutable reference (or raw pointer for
 * flat per-agent vector slices) to a single data field from mech_agent_data
 * and its sub-structs (base_membrane_data, base_motility_data,
 * base_potential_data).  No behaviour methods live here — only data accessors.
 */
class Mechanics_Agent_Interface
{
public:
	virtual ~Mechanics_Agent_Interface() = default;

	// ========================================================================
	// mech_agent_data — direct fields
	// ========================================================================

	/** @brief Pointer to the start of this agent's velocity components. */
	virtual double* get_velocity() = 0;

	virtual double& get_radius() = 0;

	virtual bool& get_is_movable() = 0;

	virtual bool& get_is_out_of_domain() = 0;

	/** @brief Neighbor agent indices (mechanics voxel adjacency list). */
	virtual std::vector<int>& get_neighbors() = 0;

	// ========================================================================
	// base_membrane_data
	// ========================================================================

	virtual double& get_cell_BM_adhesion_strength() = 0;
	virtual double& get_cell_BM_repulsion_strength() = 0;

	// ========================================================================
	// base_motility_data
	// ========================================================================

	virtual bool& get_is_motile() = 0;
	virtual double& get_persistence_time() = 0;
	virtual double& get_migration_speed() = 0;

	/** @brief Pointer to the start of this agent's migration-bias-direction components. */
	virtual double* get_migration_bias_direction() = 0;

	virtual double& get_migration_bias() = 0;

	/** @brief Pointer to the start of this agent's motility-vector components. */
	virtual double* get_motility_vector() = 0;

	virtual bool& get_restrict_to_2d() = 0;

	virtual int& get_chemotaxis_index() = 0;
	virtual int& get_chemotaxis_direction() = 0;

	/** @brief Pointer to the start of this agent's chemotactic-sensitivity values (one per substrate). */
	virtual double* get_chemotactic_sensitivities() = 0;

	// virtual std::function<void(double*)>& get_update_migration_bias_direction() = 0;

	// ========================================================================
	// base_potential_data
	// ========================================================================

	virtual double& get_cell_cell_adhesion_strength() = 0;
	virtual double& get_cell_cell_repulsion_strength() = 0;

	/** @brief Pointer to the start of this agent's per-cell-definition adhesion affinities. */
	virtual double* get_cell_adhesion_affinities() = 0;

	virtual double& get_relative_maximum_adhesion_distance() = 0;

	virtual int& get_maximum_number_of_attachments() = 0;
	virtual double& get_attachment_elastic_constant() = 0;
	virtual double& get_attachment_rate() = 0;
	virtual double& get_detachment_rate() = 0;

	virtual double& get_simple_pressure() = 0;

	/** @brief Pointer to the start of this agent's previous-velocity components. */
	virtual double* get_previous_velocity() = 0;

	/** @brief Spring-attachment partner indices for this agent. */
	virtual std::vector<int>& get_springs() = 0;
	
	virtual double& get_relative_maximum_attachment_distance() = 0; 
	virtual double& get_relative_detachment_distance() = 0; 
	virtual double& get_maximum_attachment_rate() = 0; 





	virtual int get_current_mechanics_voxel_index() = 0;


	virtual void remove_self_from_all_neighbors( void ) = 0; 
	virtual void remove_all_attached_cells( void ) = 0; // done 
	virtual void remove_all_spring_attachments( void ) = 0; // done 

	virtual void update_voxel_in_container(void) = 0;

};

} // namespace PhysiCell

#endif // __PhysiCell_mechanics_agent_interface_h__
