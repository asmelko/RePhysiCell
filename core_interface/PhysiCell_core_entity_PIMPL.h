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

	int get_ID() const override
	{ return pImpl->get_ID(); }

	void set_ID(int new_ID) override
	{ pImpl->set_ID(new_ID); }

	int get_index() const override
	{ return pImpl->get_index(); }

	void set_index(int new_index) override
	{ pImpl->set_index(new_index); }

	bool get_is_active() const override
	{ return pImpl->get_is_active(); }

	void set_is_active(bool active) override
	{ pImpl->set_is_active(active); }

	double* get_secretion_rates() override
	{ return pImpl->get_secretion_rates(); }

	const double* get_secretion_rates() const override
	{ return pImpl->get_secretion_rates(); }

	double* get_saturation_densities() override
	{ return pImpl->get_saturation_densities(); }

	const double* get_saturation_densities() const override
	{ return pImpl->get_saturation_densities(); }

	double* get_uptake_rates() override
	{ return pImpl->get_uptake_rates(); }

	const double* get_uptake_rates() const override
	{ return pImpl->get_uptake_rates(); }

	double* get_net_export_rates() override
	{ return pImpl->get_net_export_rates(); }

	const double* get_net_export_rates() const override
	{ return pImpl->get_net_export_rates(); }

	double* get_internalized_total_substrates() override
	{ return pImpl->get_internalized_total_substrates(); }

	const double* get_internalized_total_substrates() const override
	{ return pImpl->get_internalized_total_substrates(); }

	double* get_fraction_released_at_death() override
	{ return pImpl->get_fraction_released_at_death(); }

	const double* get_fraction_released_at_death() const override
	{ return pImpl->get_fraction_released_at_death(); }

	double* get_fraction_transferred_when_ingested() override
	{ return pImpl->get_fraction_transferred_when_ingested(); }

	const double* get_fraction_transferred_when_ingested() const override
	{ return pImpl->get_fraction_transferred_when_ingested(); }

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

	int get_type() const override
	{ return pImpl->get_type(); }

	void set_type(int new_type) override
	{ pImpl->set_type(new_type); }

	double* get_velocity() override
	{ return pImpl->get_velocity(); }

	double& get_radius() override
	{ return pImpl->get_radius(); }

	bool& get_is_movable() override
	{ return pImpl->get_is_movable(); }

	bool& get_is_out_of_domain() override
	{ return pImpl->get_is_out_of_domain(); }

	std::vector<CoreEntityPIMPL*>& get_neighbors() override
	{ return pImpl->get_neighbors(); }

	double& get_cell_BM_adhesion_strength() override
	{ return pImpl->get_cell_BM_adhesion_strength(); }

	double& get_cell_BM_repulsion_strength() override
	{ return pImpl->get_cell_BM_repulsion_strength(); }

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

	std::vector<std::pair<CoreEntityPIMPL*, bool>>& get_springs() override
	{ return pImpl->get_springs(); }

	double& get_relative_maximum_attachment_distance() override
	{ return pImpl->get_relative_maximum_attachment_distance(); }

	double& get_relative_detachment_distance() override
	{ return pImpl->get_relative_detachment_distance(); }

	double& get_maximum_attachment_rate() override
	{ return pImpl->get_maximum_attachment_rate(); }

	int get_current_mechanics_voxel_index() override
	{ return pImpl->get_current_mechanics_voxel_index(); }

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
