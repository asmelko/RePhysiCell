#ifndef __PhysiCell_core_entity_PIMPL_h__
#define __PhysiCell_core_entity_PIMPL_h__

#include <vector>

#include "PhysiCell_core_entity_interface.h"

namespace PhysiCell {

/**
 * @brief PIMPL base for classes (i.e. Cell) backed by a CoreEntity.
 *
 * Implements the CoreEntity interface by delegating every call to the wrapped
 * CoreEntity, which is chosen at runtime by the installed CoreImplementation.
 * The wrapped entity is owned and deleted by this class.
 */
class CoreEntityPIMPL : public CoreEntity
{
protected:
	CoreEntity* pImpl = nullptr;

public:
	CoreEntityPIMPL() = default;
	explicit CoreEntityPIMPL(CoreEntity* core_entity);
	virtual ~CoreEntityPIMPL();

	CoreEntityPIMPL(const CoreEntityPIMPL&) = delete;
	CoreEntityPIMPL& operator=(const CoreEntityPIMPL&) = delete;

	CoreEntity* get_core_entity()
	{ return pImpl; }

	const CoreEntity* get_core_entity() const
	{ return pImpl; }

	double* position() override
	{ return pImpl->position(); }

	Mechanics_Agent_Interface* get_mechanics_implementation() override
	{ return pImpl->get_mechanics_implementation(); }

	const Mechanics_Agent_Interface* get_mechanics_implementation() const
	{ return pImpl->get_mechanics_implementation(); }

	// ---- BioFVM::Basic_Agent_Interface ---------------------------------------

	double& get_total_volume() override
	{ return pImpl->get_total_volume(); }

	void set_total_volume(double a0) override
	{ pImpl->set_total_volume(a0); }

	void update_voxel_index() override
	{ pImpl->update_voxel_index(); }

	void release_internalized_substrates() override
	{ pImpl->release_internalized_substrates(); }

	void set_internal_uptake_constants(double dt) override
	{ pImpl->set_internal_uptake_constants(dt); }

	void register_microenvironment(BioFVM::Microenvironment_Interface* a0) override
	{ pImpl->register_microenvironment(a0); }

	BioFVM::Microenvironment_Interface* get_microenvironment_interface() override
	{ return pImpl->get_microenvironment_interface(); }

	int& ID() override
	{ return pImpl->ID(); }

	int& index() override
	{ return pImpl->index(); }

	bool& is_active() override
	{ return pImpl->is_active(); }

	double* secretion_rates() override
	{ return pImpl->secretion_rates(); }

	double* saturation_densities() override
	{ return pImpl->saturation_densities(); }

	double* uptake_rates() override
	{ return pImpl->uptake_rates(); }

	double* net_export_rates() override
	{ return pImpl->net_export_rates(); }

	double* internalized_total_substrates() override
	{ return pImpl->internalized_total_substrates(); }

	double* fraction_released_at_death() override
	{ return pImpl->fraction_released_at_death(); }

	double* fraction_transferred_when_ingested() override
	{ return pImpl->fraction_transferred_when_ingested(); }

	void simulate_secretion_and_uptake(double dt) override
	{ pImpl->simulate_secretion_and_uptake(dt); }

	int get_current_voxel_index() override
	{ return pImpl->get_current_voxel_index(); }

	double* nearest_density_vector() override
	{ return pImpl->nearest_density_vector(); }

	std::vector<double>& nearest_gradient(int substrate_index) override
	{ return pImpl->nearest_gradient(substrate_index); }

	std::vector<std::vector<double>>& nearest_gradient_vector() override
	{ return pImpl->nearest_gradient_vector(); }

	// ---- Mechanics_Agent_Interface -------------------------------------------

	int& type() override
	{ return pImpl->type(); }

	double* velocity() override
	{ return pImpl->velocity(); }

	double& radius() override
	{ return pImpl->radius(); }

	bool& is_movable() override
	{ return pImpl->is_movable(); }

	bool& is_out_of_domain() override
	{ return pImpl->is_out_of_domain(); }

	std::vector<CoreEntityPIMPL*>& neighbors() override
	{ return pImpl->neighbors(); }

	double& cell_BM_adhesion_strength() override
	{ return pImpl->cell_BM_adhesion_strength(); }

	double& cell_BM_repulsion_strength() override
	{ return pImpl->cell_BM_repulsion_strength(); }

	bool& is_motile() override
	{ return pImpl->is_motile(); }

	double& persistence_time() override
	{ return pImpl->persistence_time(); }

	double& migration_speed() override
	{ return pImpl->migration_speed(); }

	double* migration_bias_direction() override
	{ return pImpl->migration_bias_direction(); }

	double& migration_bias() override
	{ return pImpl->migration_bias(); }

	double* motility_vector() override
	{ return pImpl->motility_vector(); }

	bool& restrict_to_2d() override
	{ return pImpl->restrict_to_2d(); }

	int& chemotaxis_index() override
	{ return pImpl->chemotaxis_index(); }

	int& chemotaxis_direction() override
	{ return pImpl->chemotaxis_direction(); }

	double* chemotactic_sensitivities() override
	{ return pImpl->chemotactic_sensitivities(); }

	double& cell_cell_adhesion_strength() override
	{ return pImpl->cell_cell_adhesion_strength(); }

	double& cell_cell_repulsion_strength() override
	{ return pImpl->cell_cell_repulsion_strength(); }

	double* cell_adhesion_affinities() override
	{ return pImpl->cell_adhesion_affinities(); }

	double& relative_maximum_adhesion_distance() override
	{ return pImpl->relative_maximum_adhesion_distance(); }

	int& maximum_number_of_attachments() override
	{ return pImpl->maximum_number_of_attachments(); }

	double& attachment_elastic_constant() override
	{ return pImpl->attachment_elastic_constant(); }

	double& attachment_rate() override
	{ return pImpl->attachment_rate(); }

	double& detachment_rate() override
	{ return pImpl->detachment_rate(); }

	double& simple_pressure() override
	{ return pImpl->simple_pressure(); }

	double* previous_velocity() override
	{ return pImpl->previous_velocity(); }

	std::vector<std::pair<CoreEntityPIMPL*, bool>>& springs() override
	{ return pImpl->springs(); }

	double& relative_maximum_attachment_distance() override
	{ return pImpl->relative_maximum_attachment_distance(); }

	double& relative_detachment_distance() override
	{ return pImpl->relative_detachment_distance(); }

	double& maximum_attachment_rate() override
	{ return pImpl->maximum_attachment_rate(); }

	int current_mechanics_voxel_index() override
	{ return pImpl->current_mechanics_voxel_index(); }

	void remove_self_from_all_neighbors() override
	{ pImpl->remove_self_from_all_neighbors(); }

	void remove_all_spring_attachments() override
	{ pImpl->remove_all_spring_attachments(); }

	void update_voxel_in_container() override
	{ pImpl->update_voxel_in_container(); }

	bool assign_position(double x, double y, double z) override
	{ return pImpl->assign_position(x, y, z); }

	bool assign_position(const std::vector<double>& new_position) override
	{ return pImpl->assign_position(new_position); }
};

} // namespace PhysiCell

#endif // __PhysiCell_core_entity_PIMPL_h__
