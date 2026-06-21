#ifndef __PhysiCell_mechanics_agent_PIMPL_h__
#define __PhysiCell_mechanics_agent_PIMPL_h__

#include <vector>

#include "PhysiCell_mechanics_agent_interface.h"

namespace PhysiCell {

/**
 * @brief PIMPL base class for mechanics agents.
 *
 * Mirrors the BioFVM::Basic_Agent_PIMPL pattern for the mechanics subsystem.
 * Derived classes (e.g., PhysiCell::Cell) inherit from this and automatically
 * satisfy Mechanics_Agent_Interface by delegating every call to the wrapped
 * pImpl pointer.
 */
class Mechanics_Agent_PIMPL : public Mechanics_Agent_Interface
{
protected:
	/// Wrapped implementation — owned by this PIMPL; deleted in destructor.
	Mechanics_Agent_Interface* pImpl = nullptr;

public:
	Mechanics_Agent_PIMPL();
	explicit Mechanics_Agent_PIMPL(Mechanics_Agent_Interface* impl);
	virtual ~Mechanics_Agent_PIMPL();

	Mechanics_Agent_Interface* get_mechanics_implementation()
	{ return pImpl; }

	const Mechanics_Agent_Interface* get_mechanics_implementation() const
	{ return pImpl; }

	// ========================================================================
	// Mechanics_Agent_Interface — all delegate to pImpl
	// ========================================================================

	// ---- mech_agent_data — direct fields -----------------------------------

	double* get_velocity() override
	{ return pImpl->get_velocity(); }

	double& get_radius() override
	{ return pImpl->get_radius(); }

	bool& get_is_movable() override
	{ return pImpl->get_is_movable(); }

	std::vector<int>& get_neighbors() override
	{ return pImpl->get_neighbors(); }

	// ---- base_membrane_data ------------------------------------------------

	double& get_cell_BM_adhesion_strength() override
	{ return pImpl->get_cell_BM_adhesion_strength(); }

	double& get_cell_BM_repulsion_strength() override
	{ return pImpl->get_cell_BM_repulsion_strength(); }

	// ---- base_motility_data ------------------------------------------------

	bool& get_is_motile() override
	{ return pImpl->get_is_motile(); }

	double& get_persistence_time() override
	{ return pImpl->get_persistence_time(); }

	double& get_migration_speed() override
	{ return pImpl->get_migration_speed(); }

	double* get_migration_bias_direction() override
	{ return pImpl->get_migration_bias_direction(); }

	double& get_migration_bias() override
	{ return pImpl->get_migration_bias(); }

	double* get_motility_vector() override
	{ return pImpl->get_motility_vector(); }

	bool& get_restrict_to_2d() override
	{ return pImpl->get_restrict_to_2d(); }

	int& get_chemotaxis_index() override
	{ return pImpl->get_chemotaxis_index(); }

	int& get_chemotaxis_direction() override
	{ return pImpl->get_chemotaxis_direction(); }

	double* get_chemotactic_sensitivities() override
	{ return pImpl->get_chemotactic_sensitivities(); }

	// ---- base_potential_data -----------------------------------------------

	double& get_cell_cell_adhesion_strength() override
	{ return pImpl->get_cell_cell_adhesion_strength(); }

	double& get_cell_cell_repulsion_strength() override
	{ return pImpl->get_cell_cell_repulsion_strength(); }

	double* get_cell_adhesion_affinities() override
	{ return pImpl->get_cell_adhesion_affinities(); }

	double& get_relative_maximum_adhesion_distance() override
	{ return pImpl->get_relative_maximum_adhesion_distance(); }

	int& get_maximum_number_of_attachments() override
	{ return pImpl->get_maximum_number_of_attachments(); }

	double& get_attachment_elastic_constant() override
	{ return pImpl->get_attachment_elastic_constant(); }

	double& get_attachment_rate() override
	{ return pImpl->get_attachment_rate(); }

	double& get_detachment_rate() override
	{ return pImpl->get_detachment_rate(); }

	double& get_simple_pressure() override
	{ return pImpl->get_simple_pressure(); }

	double* get_previous_velocity() override
	{ return pImpl->get_previous_velocity(); }

	std::vector<int>& get_springs() override
	{ return pImpl->get_springs(); }
	
	
	double& get_relative_maximum_attachment_distance() override
	{ return pImpl->get_relative_maximum_attachment_distance(); }

	double& get_relative_detachment_distance() override
	{ return pImpl->get_relative_detachment_distance(); }

	double& get_maximum_attachment_rate() override
	{ return pImpl->get_maximum_attachment_rate(); }


	int get_current_mechanics_voxel_index() override
	{ return pImpl->get_current_mechanics_voxel_index(); }

	void remove_self_from_all_neighbors( void ) override
	{ pImpl->remove_self_from_all_neighbors(); }

	void remove_all_attached_cells( void ) override
	{ pImpl->remove_all_attached_cells(); }

	void remove_all_spring_attachments( void ) override
	{ pImpl->remove_all_spring_attachments(); }

	bool& get_is_out_of_domain() override
	{ return pImpl->get_is_out_of_domain(); }

	void update_voxel_in_container(void) override
	{ pImpl->update_voxel_in_container(); }

};

} // namespace PhysiCell

#endif // __PhysiCell_mechanics_agent_PIMPL_h__
