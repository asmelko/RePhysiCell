#ifndef __PhysiCell_mechanics_environment_interface_h__
#define __PhysiCell_mechanics_environment_interface_h__

#include <vector>
#include "../BioFVM/BioFVM_mesh.h"
#include "PhysiCell_mechanics_agent_interface.h"
#include "PhysiCell_mechanics_standard_models_interface.h"

namespace PhysiCell {

/**
 * @brief Abstract interface for the mechanics environment.
 *
 * Analogous to BioFVM::Microenvironment_Interface for the diffusion solver,
 * Mechanics_Environment_Interface abstracts the spatial grid and all core
 * mechanical operations: neighbor creation, pairwise potentials, motility,
 * basement-membrane forces, and Adams-Bashforth position integration.
 *
 * The interface is intentionally narrow: it holds the mechanics mesh and
 * exposes the five required mechanical operations plus grid lifecycle helpers.
 * Concrete implementations (e.g., mechanics_environment) own the agent grid
 * and all associated data.
 */
class Mechanics_Environment_Interface : public Mechanics_Standard_Models_Interface
{
public:
	virtual ~Mechanics_Environment_Interface() = default;

	// ========================================================================
	// Mesh access
	// ========================================================================

	/** @brief The Cartesian mesh that defines the mechanics spatial grid. */
	virtual BioFVM::Cartesian_Mesh& get_mechanics_mesh() = 0;
	virtual const BioFVM::Cartesian_Mesh& get_mechanics_mesh() const = 0;

	/** @brief Per-voxel maximum interactive distance (used for neighbor voxel filtering). */
	virtual std::vector<double>& get_max_cell_interactive_distance_in_voxel() = 0;

	// ========================================================================
	// Initialization
	// ========================================================================

	virtual void initialize(double x_start, double x_end,
	                        double y_start, double y_end,
	                        double z_start, double z_end,
	                        double voxel_size) = 0;

	virtual void initialize(double x_start, double x_end,
	                        double y_start, double y_end,
	                        double z_start, double z_end,
	                        double dx, double dy, double dz) = 0;

	// ========================================================================
	// Core mechanical operations (the five required by the interface contract)
	// ========================================================================

	/**
	 * @brief Compute repulsive and adhesive pairwise potential between agents in neighboring voxels.
	 *
	 * Internally calls agent->add_potentials(other) for each neighbor pair.
	 * Also updates simple pressure and the list of neighbors for each agent.
	 * Also computes membrane interactions and motility contributions.
	 */
	virtual void compute_velocities(double dt) = 0;

	/**
	 * @brief Compute spring attachment forces and add them to the agent velocity.
	 *
	 * Updates attachments based on the attachment rates
	 * and adds spring forces to the velocity for each attached neighbor.
	 */
	virtual void compute_spring_attachments(double dt) = 0;

	/**
	 * @brief Adams-Bashforth position update for all agents.
	 *
	 * Integrates position using current and previous velocities, resets current
	 * velocity to zero, and updates the pending mechanics voxel index.
	 */
	virtual void update_positions(double dt) = 0;

	/**
	 * @brief Update the mechanics environment after all position updates are done.
	 *
	 * Updates the voxel lists of agents based on their updated positions.
	 */
	virtual void update_container() = 0;

	// ========================================================================
	// Grid lifecycle helpers
	// ========================================================================

	virtual void add_agent_to_voxel(Mechanics_Agent_PIMPL* agent, int voxel_index) = 0;
	virtual void remove_agent_from_voxel(Mechanics_Agent_PIMPL* agent, int voxel_index) = 0;
	virtual void add_agent_to_outer_voxel(Mechanics_Agent_PIMPL* agent) = 0;

	virtual void register_agent(Mechanics_Agent_PIMPL* agent) = 0;
	virtual void remove_agent(Mechanics_Agent_PIMPL* agent) = 0;
	virtual bool contain_any_cell(int voxel_index) = 0;




	virtual void attach_cells_as_spring( Mechanics_Agent_PIMPL* pCell_1, Mechanics_Agent_PIMPL* pCell_2, bool attacking_spring = false ) = 0;
	virtual void detach_cells_as_spring( Mechanics_Agent_PIMPL* pCell_1 , Mechanics_Agent_PIMPL* pCell_2 ) = 0;

	virtual bool is_neighbor_voxel(Mechanics_Agent_PIMPL* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index) = 0;

	virtual std::vector<Mechanics_Agent_PIMPL*> agents_in_my_container( Mechanics_Agent_PIMPL* pCell ) = 0;
	virtual std::vector<Mechanics_Agent_PIMPL*> nearby_agents( Mechanics_Agent_PIMPL* pCell ) = 0;
	virtual std::vector<Mechanics_Agent_PIMPL*> nearby_interacting_agents( Mechanics_Agent_PIMPL* pCell ) = 0;

};

} // namespace PhysiCell

#endif // __PhysiCell_mechanics_environment_interface_h__
