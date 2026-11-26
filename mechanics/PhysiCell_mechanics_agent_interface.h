#ifndef __PhysiCell_mechanics_agent_interface_h__
#define __PhysiCell_mechanics_agent_interface_h__

#include <vector>


namespace PhysiCell {

class CoreEntityPIMPL; // forward declaration

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
	
	virtual int& type() = 0;

	/** @brief Pointer to the start of this agent's velocity components. */
	virtual double* velocity() = 0;

	virtual double& radius() = 0;

	virtual bool& is_movable() = 0;

	virtual bool& is_out_of_domain() = 0;

	/** @brief Neighbor agent indices (mechanics voxel adjacency list). */
	virtual std::vector<CoreEntityPIMPL*>& neighbors() = 0;

	// ========================================================================
	// base_membrane_data
	// ========================================================================

	virtual double& cell_BM_adhesion_strength() = 0;
	virtual double& cell_BM_repulsion_strength() = 0;

	// ========================================================================
	// base_motility_data
	// ========================================================================

	virtual bool& is_motile() = 0;
	virtual double& persistence_time() = 0;
	virtual double& migration_speed() = 0;

	/** @brief Pointer to the start of this agent's migration-bias-direction components. */
	virtual double* migration_bias_direction() = 0;

	virtual double& migration_bias() = 0;

	/** @brief Pointer to the start of this agent's motility-vector components. */
	virtual double* motility_vector() = 0;

	virtual bool& restrict_to_2d() = 0;

	virtual int& chemotaxis_index() = 0;
	virtual int& chemotaxis_direction() = 0;

	/** @brief Pointer to the start of this agent's chemotactic-sensitivity values (one per substrate). */
	virtual double* chemotactic_sensitivities() = 0;

	// virtual std::function<void(double*)>& get_update_migration_bias_direction() = 0;

	// ========================================================================
	// base_potential_data
	// ========================================================================

	virtual double& cell_cell_adhesion_strength() = 0;
	virtual double& cell_cell_repulsion_strength() = 0;

	/** @brief Pointer to the start of this agent's per-cell-definition adhesion affinities. */
	virtual double* cell_adhesion_affinities() = 0;

	virtual double& relative_maximum_adhesion_distance() = 0;

	virtual int& maximum_number_of_attachments() = 0;
	virtual double& attachment_elastic_constant() = 0;
	virtual double& attachment_rate() = 0;
	virtual double& detachment_rate() = 0;

	virtual double& simple_pressure() = 0;

	/** @brief Pointer to the start of this agent's previous-velocity components. */
	virtual double* previous_velocity() = 0;

	/** @brief Spring-attachment partner indices for this agent. */
	virtual std::vector<std::pair<CoreEntityPIMPL*, bool>>& springs() = 0;
	
	virtual double& relative_maximum_attachment_distance() = 0; 
	virtual double& relative_detachment_distance() = 0; 
	virtual double& maximum_attachment_rate() = 0; 





	virtual int current_mechanics_voxel_index() = 0;


	virtual void remove_self_from_all_neighbors( void ) = 0; 
	virtual void remove_all_spring_attachments( void ) = 0; // done 

	virtual void update_voxel_in_container(void) = 0;

	virtual bool assign_position(double x, double y, double z) = 0;
	virtual bool assign_position(const std::vector<double>& new_position) = 0;
};

} // namespace PhysiCell

#endif // __PhysiCell_mechanics_agent_interface_h__
