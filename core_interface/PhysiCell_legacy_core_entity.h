#ifndef __PhysiCell_legacy_core_entity_h__
#define __PhysiCell_legacy_core_entity_h__

#include <array>
#include <memory>
#include <vector>

#include "PhysiCell_core_entity_interface.h"
#include "../BioFVM/BioFVM_basic_agent.h"
#include "../mechanics/PhysiCell_mechanics_agent.h"

namespace PhysiCell {

class Cell;

/**
 * @brief Core entity of the legacy implementation.
 *
 * Owns the position storage together with a legacy BioFVM::Basic_Agent and a legacy
 * Mechanics_Agent, both of which read and write the position held here.
 * Every interface method is forwarded to the respective agent.
 */
class LegacyCoreEntity : public CoreEntity
{
	std::array<double, 3> position_ = {0.0, 0.0, 0.0};
	std::unique_ptr<BioFVM::Basic_Agent> basic_agent_;
	std::unique_ptr<Mechanics_Agent> mechanics_agent_;

public:
	explicit LegacyCoreEntity(Cell* pCell);
	/** @brief Use a custom mechanics agent (e.g. an extension of Mechanics_Agent) instead of the default one. */
	LegacyCoreEntity(std::unique_ptr<Mechanics_Agent> mechanics_agent);
	~LegacyCoreEntity() override;

	double* position() override
	{ return position_.data(); }

	BioFVM::Basic_Agent* get_basic_agent()
	{ return basic_agent_.get(); }

	Mechanics_Agent* get_mechanics_agent()
	{ return mechanics_agent_.get(); }

	Mechanics_Agent_Interface* get_mechanics_implementation() override;

	/** @brief Replace the mechanics agent; the new one is bound to this entity's position. */
	void set_mechanics_agent(std::unique_ptr<Mechanics_Agent> mechanics_agent);

	// ---- BioFVM::Basic_Agent_Interface ---------------------------------------

	double& get_total_volume() override
	{ return basic_agent_->get_total_volume(); }

	void set_total_volume(double a0) override
	{ basic_agent_->set_total_volume(a0); }

	void update_voxel_index() override
	{ basic_agent_->update_voxel_index(); }

	void release_internalized_substrates() override
	{ basic_agent_->release_internalized_substrates(); }

	void set_internal_uptake_constants(double dt) override
	{ basic_agent_->set_internal_uptake_constants(dt); }

	void register_microenvironment(BioFVM::Microenvironment_Interface* a0) override
	{ basic_agent_->register_microenvironment(a0); }

	BioFVM::Microenvironment_Interface* get_microenvironment_interface() override
	{ return basic_agent_->get_microenvironment_interface(); }

	int get_ID() const override
	{ return basic_agent_->get_ID(); }

	void set_ID(int new_ID) override
	{ basic_agent_->set_ID(new_ID); }

	int get_index() const override
	{ return basic_agent_->get_index(); }

	void set_index(int new_index) override
	{ basic_agent_->set_index(new_index); }

	bool get_is_active() const override
	{ return basic_agent_->get_is_active(); }

	void set_is_active(bool active) override
	{ basic_agent_->set_is_active(active); }

	double* get_secretion_rates() override
	{ return basic_agent_->get_secretion_rates(); }

	const double* get_secretion_rates() const override
	{ return basic_agent_->get_secretion_rates(); }

	double* get_saturation_densities() override
	{ return basic_agent_->get_saturation_densities(); }

	const double* get_saturation_densities() const override
	{ return basic_agent_->get_saturation_densities(); }

	double* get_uptake_rates() override
	{ return basic_agent_->get_uptake_rates(); }

	const double* get_uptake_rates() const override
	{ return basic_agent_->get_uptake_rates(); }

	double* get_net_export_rates() override
	{ return basic_agent_->get_net_export_rates(); }

	const double* get_net_export_rates() const override
	{ return basic_agent_->get_net_export_rates(); }

	double* get_internalized_total_substrates() override
	{ return basic_agent_->get_internalized_total_substrates(); }

	const double* get_internalized_total_substrates() const override
	{ return basic_agent_->get_internalized_total_substrates(); }

	double* get_fraction_released_at_death() override
	{ return basic_agent_->get_fraction_released_at_death(); }

	const double* get_fraction_released_at_death() const override
	{ return basic_agent_->get_fraction_released_at_death(); }

	double* get_fraction_transferred_when_ingested() override
	{ return basic_agent_->get_fraction_transferred_when_ingested(); }

	const double* get_fraction_transferred_when_ingested() const override
	{ return basic_agent_->get_fraction_transferred_when_ingested(); }

	void simulate_secretion_and_uptake(double dt) override
	{ basic_agent_->simulate_secretion_and_uptake(dt); }

	int get_current_voxel_index() override
	{ return basic_agent_->get_current_voxel_index(); }

	double* nearest_density_vector() override
	{ return basic_agent_->nearest_density_vector(); }

	std::vector<double>& nearest_gradient(int substrate_index) override
	{ return basic_agent_->nearest_gradient(substrate_index); }

	std::vector<std::vector<double>>& nearest_gradient_vector() override
	{ return basic_agent_->nearest_gradient_vector(); }

	// ---- Mechanics_Agent_Interface -------------------------------------------

	int get_type() const override
	{ return mechanics_agent_->get_type(); }

	void set_type(int new_type) override
	{ mechanics_agent_->set_type(new_type); }

	double* get_velocity() override
	{ return mechanics_agent_->get_velocity(); }

	double& get_radius() override
	{ return mechanics_agent_->get_radius(); }

	bool& get_is_movable() override
	{ return mechanics_agent_->get_is_movable(); }

	bool& get_is_out_of_domain() override
	{ return mechanics_agent_->get_is_out_of_domain(); }

	std::vector<CoreEntityPIMPL*>& get_neighbors() override
	{ return mechanics_agent_->get_neighbors(); }

	double& get_cell_BM_adhesion_strength() override
	{ return mechanics_agent_->get_cell_BM_adhesion_strength(); }

	double& get_cell_BM_repulsion_strength() override
	{ return mechanics_agent_->get_cell_BM_repulsion_strength(); }

	bool& get_is_motile() override
	{ return mechanics_agent_->get_is_motile(); }

	double& get_persistence_time() override
	{ return mechanics_agent_->get_persistence_time(); }

	double& get_migration_speed() override
	{ return mechanics_agent_->get_migration_speed(); }

	double* get_migration_bias_direction() override
	{ return mechanics_agent_->get_migration_bias_direction(); }

	double& get_migration_bias() override
	{ return mechanics_agent_->get_migration_bias(); }

	double* get_motility_vector() override
	{ return mechanics_agent_->get_motility_vector(); }

	bool& get_restrict_to_2d() override
	{ return mechanics_agent_->get_restrict_to_2d(); }

	int& get_chemotaxis_index() override
	{ return mechanics_agent_->get_chemotaxis_index(); }

	int& get_chemotaxis_direction() override
	{ return mechanics_agent_->get_chemotaxis_direction(); }

	double* get_chemotactic_sensitivities() override
	{ return mechanics_agent_->get_chemotactic_sensitivities(); }

	double& get_cell_cell_adhesion_strength() override
	{ return mechanics_agent_->get_cell_cell_adhesion_strength(); }

	double& get_cell_cell_repulsion_strength() override
	{ return mechanics_agent_->get_cell_cell_repulsion_strength(); }

	double* get_cell_adhesion_affinities() override
	{ return mechanics_agent_->get_cell_adhesion_affinities(); }

	double& get_relative_maximum_adhesion_distance() override
	{ return mechanics_agent_->get_relative_maximum_adhesion_distance(); }

	int& get_maximum_number_of_attachments() override
	{ return mechanics_agent_->get_maximum_number_of_attachments(); }

	double& get_attachment_elastic_constant() override
	{ return mechanics_agent_->get_attachment_elastic_constant(); }

	double& get_attachment_rate() override
	{ return mechanics_agent_->get_attachment_rate(); }

	double& get_detachment_rate() override
	{ return mechanics_agent_->get_detachment_rate(); }

	double& get_simple_pressure() override
	{ return mechanics_agent_->get_simple_pressure(); }

	double* get_previous_velocity() override
	{ return mechanics_agent_->get_previous_velocity(); }

	std::vector<std::pair<CoreEntityPIMPL*, bool>>& get_springs() override
	{ return mechanics_agent_->get_springs(); }

	double& get_relative_maximum_attachment_distance() override
	{ return mechanics_agent_->get_relative_maximum_attachment_distance(); }

	double& get_relative_detachment_distance() override
	{ return mechanics_agent_->get_relative_detachment_distance(); }

	double& get_maximum_attachment_rate() override
	{ return mechanics_agent_->get_maximum_attachment_rate(); }

	int get_current_mechanics_voxel_index() override
	{ return mechanics_agent_->get_current_mechanics_voxel_index(); }

	void remove_self_from_all_neighbors() override
	{ mechanics_agent_->remove_self_from_all_neighbors(); }

	void remove_all_spring_attachments() override
	{ mechanics_agent_->remove_all_spring_attachments(); }

	void update_voxel_in_container() override
	{ mechanics_agent_->update_voxel_in_container(); }

	bool assign_position(double x, double y, double z) override
	{ return mechanics_agent_->assign_position(x, y, z); }

	bool assign_position(const std::vector<double>& new_position) override
	{ return mechanics_agent_->assign_position(new_position); }
};

} // namespace PhysiCell

#endif // __PhysiCell_legacy_core_entity_h__
