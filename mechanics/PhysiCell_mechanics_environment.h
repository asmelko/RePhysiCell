#ifndef __PhysiCell_mechanics_environment_h__
#define __PhysiCell_mechanics_environment_h__

#include "PhysiCell_mechanics_agent_interface.h"
#include "PhysiCell_mechanics_environment_interface.h"
#include "../BioFVM/BioFVM_mesh.h"
#include <vector>
#include "PhysiCell_mechanics_agent.h"
#include "PhysiCell_mechanics_standard_models.h"

namespace PhysiCell {

class Cell; // forward declaration — actual type used for the agent grid

/**
 * @brief Concrete mechanics environment.
 *
 * Owns the Cartesian mechanics mesh, the per-voxel agent grid, and implements
 * all five mechanical operations defined by Mechanics_Environment_Interface:
 *
 *   create_neighbors               — BM forces + all pairwise potentials + motility
 *   compute_pairwise_potential     — repulsion / adhesion between two agents
 *   compute_motility_potential     — stochastic directed migration
 *   compute_basement_membrane_potential — BM adhesion / repulsion
 *   update_position                — Adams-Bashforth time integration
 *
 * The agent grid stores Cell* pointers (the only concrete agent type in PhysiCell).
 * Cell_Container holds a mechanics_environment as a value member and exposes
 * backward-compatible public references to underlying_mesh, agent_grid, etc.
 */
class mechanics_environment : public Mechanics_Environment_Interface
{
private:
	bool simulate_2D_ = false;

public:
	// ---- Data owned by the mechanics environment ---------------------------

	std::vector<Cell*> *all_agents; // alias to Cell_Container::all_cells

	BioFVM::Cartesian_Mesh mechanics_mesh;

	/** @brief Per-voxel lists of cells. Grid is resized in initialize(). */
	std::vector<std::vector<Mechanics_Agent*>> agent_grid;

	/** @brief Lists of cells that have moved outside domain boundaries. */
	std::vector<std::vector<Mechanics_Agent*>> agents_in_outer_voxels;

	/** @brief Per-voxel maximum interaction radius (updated on volume changes). */
	std::vector<double> max_cell_interactive_distance_in_voxel;

	bool disable_automated_spring_adhesions;

	Mechanics_Standard_Models models;

	// ---- Mechanics_Environment_Interface -----------------------------------

	BioFVM::Cartesian_Mesh& get_mechanics_mesh() override { return mechanics_mesh; }
	const BioFVM::Cartesian_Mesh& get_mechanics_mesh() const override { return mechanics_mesh; }

	std::vector<double>& get_max_cell_interactive_distance_in_voxel() override
	{ return max_cell_interactive_distance_in_voxel; }

	void initialize(double x_start, double x_end,
	                double y_start, double y_end,
	                double z_start, double z_end,
	                double voxel_size) override;

	void initialize(double x_start, double x_end,
	                double y_start, double y_end,
	                double z_start, double z_end,
	                double dx, double dy, double dz) override;

	virtual void compute_velocities(double dt) override;
	virtual void compute_spring_attachments(double dt) override;
	virtual void update_positions(double dt) override;
	virtual void update_container() override;

	void add_agent_to_voxel(Mechanics_Agent_PIMPL* agent, int voxel_index) override;
	void add_agent_to_voxel(Mechanics_Agent* agent, int voxel_index);
	void remove_agent_from_voxel(Mechanics_Agent_PIMPL* agent, int voxel_index) override;
	void remove_agent_from_voxel(Mechanics_Agent* agent, int voxel_index);
	void add_agent_to_outer_voxel(Mechanics_Agent_PIMPL* agent) override;
	void add_agent_to_outer_voxel(Mechanics_Agent* agent);

	void register_agent( Mechanics_Agent_PIMPL* agent ) override;
	void register_agent( Mechanics_Agent* agent );
	void remove_agent(Mechanics_Agent_PIMPL* agent ) override;
	void remove_agent(Mechanics_Agent* agent );
	bool contain_any_cell(int voxel_index) override;

	int find_escaping_face_index(Mechanics_Agent* agent);


	void attach_cells_as_spring( Mechanics_Agent_PIMPL* pCell_1, Mechanics_Agent_PIMPL* pCell_2, bool attacking_spring ) override;
	void attach_cells_as_spring( Mechanics_Agent* pCell_1, Mechanics_Agent* pCell_2, bool attacking_spring );
	void detach_cells_as_spring( Mechanics_Agent_PIMPL* pCell_1 , Mechanics_Agent_PIMPL* pCell_2 ) override;
	void detach_cells_as_spring( Mechanics_Agent* pCell_1 , Mechanics_Agent* pCell_2 );

	bool is_neighbor_voxel(Mechanics_Agent_PIMPL* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index) override;
	bool is_neighbor_voxel(Mechanics_Agent* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index);

	void standard_add_basement_membrane_interactions( Mechanics_Agent_PIMPL* pCell, double dt ) override;
    void standard_domain_edge_avoidance_interactions( Mechanics_Agent_PIMPL* pCell, double dt ) override;
    double distance_to_domain_edge(Mechanics_Agent_PIMPL* pCell) override;

    void dynamic_spring_attachments( Mechanics_Agent_PIMPL* pCell , double dt ) override;
    void dynamic_spring_attachments( Mechanics_Agent* pCell , double dt );
    void standard_elastic_contact_function( Mechanics_Agent_PIMPL* pC1, Mechanics_Agent_PIMPL* pC2, double dt ) override;
    void standard_elastic_contact_function( Mechanics_Agent* pC1, Mechanics_Agent* pC2, double dt );
    void standard_elastic_contact_function_confluent_rest_length( Mechanics_Agent_PIMPL* pC1I, Mechanics_Agent_PIMPL* pC2I, double dt ) override;

    void chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt ) override;
    void advanced_chemotaxis_function_normalized( Mechanics_Agent_PIMPL* pCellI, double dt ) override;
    void advanced_chemotaxis_function( Mechanics_Agent_PIMPL* pCellI, double dt ) override; 

	std::vector<Mechanics_Agent_PIMPL*> agents_in_my_container( Mechanics_Agent_PIMPL* pCell ) override;
	std::vector<Mechanics_Agent_PIMPL*> nearby_agents( Mechanics_Agent_PIMPL* pCell ) override;
	std::vector<Mechanics_Agent_PIMPL*> nearby_interacting_agents( Mechanics_Agent_PIMPL* pCell ) override;

	std::vector<Mechanics_Agent_PIMPL*> agents_in_my_container( Mechanics_Agent* pCell );
	std::vector<Mechanics_Agent_PIMPL*> nearby_agents( Mechanics_Agent* pCell );
	std::vector<Mechanics_Agent_PIMPL*> nearby_interacting_agents( Mechanics_Agent* pCell );

};

extern mechanics_environment mech_environment; // global instance of the mechanics environment
mechanics_environment& get_mechanics_environment(); // global accessor for the mechanics environment


} // namespace PhysiCell

#endif // __PhysiCell_mechanics_environment_h__
