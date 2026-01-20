#ifndef __physicore_basic_agent_adapter_h__
#define __physicore_basic_agent_adapter_h__

#include <biofvm/agent.h>
#include <biofvm/microenvironment.h>
#include <cassert>
#include "../../BioFVM/BioFVM_basic_agent_interface.h"

namespace BioFVM {

class physicore_microenvironment_adapter;

/**
 * @brief Lightweight stack-based proxy adapter bridging BioFVM's Basic_Agent_Interface 
 *        with physicore's Structure-of-Arrays agent storage
 *
 * This adapter wraps an agent index and references to physicore's SoA data arrays,
 * providing BioFVM's mutable interface over physicore's immutable backend. It stores
 * BioFVM-specific state (velocity, ID, type) that doesn't exist in physicore.
 *
 * @note IMPORTANT USAGE CONSTRAINTS:
 * - simulate_secretion_and_uptake() NOT SUPPORTED - physicore uses centralized solver (throws exception)
 */
class physicore_basic_agent_adapter : public Basic_Agent_Interface
{
private:
	// Reference to physicore agent container
	physicore::biofvm::agent* agent;
	int index;
	
	physicore_microenvironment_adapter* microenvironment;

	static int ID_counter;
	
	// BioFVM-specific state not in physicore
	int ID;
	int type;
	int cached_voxel_index;
	
	// Cached position vector so we always return 3D vector
	// even in 2D simulations
	mutable std::vector<double> position_cache;

	std::vector<double> velocity;
	std::vector<double> previous_velocity;

	std::vector<std::vector<double>> gradient_cache;
	void compute_gradient_at_voxel(int voxel_index, int substrate_index, std::vector<double>& gradient) const;
	
protected:
	double* get_position_internal() override;

public:
	/**
	 * @brief Construct agent proxy wrapping physicore agent at given index
	 */
	physicore_basic_agent_adapter();
	
	virtual ~physicore_basic_agent_adapter();

	// ========================================================================
	// Volume methods
	// ========================================================================
	double& get_total_volume() override;
	void set_total_volume(double volume) override;
	void update_voxel_index() override;

	// ========================================================================
	// Microenvironment registration and access
	// ========================================================================
	void register_microenvironment(Microenvironment_Interface* microenv) override;
	Microenvironment_Interface* get_microenvironment_interface() override;

	// ========================================================================
	// ID and type accessors
	// ========================================================================
	int get_ID() const override;
	void set_ID(int new_ID) override;
	int get_index() const override;
	void set_index(int new_index) override;
	int get_type() const override;
	void set_type(int new_type) override;
	
	// ========================================================================
	// Position methods
	// ========================================================================
	bool assign_position(double x, double y, double z) override;
	bool assign_position(std::vector<double> new_position) override;
	const std::vector<double>& get_position() const override;
	void update_position(double dt) override;
	
	// ========================================================================
	// Velocity methods
	// ========================================================================
	std::vector<double>& get_velocity() override;
	const std::vector<double>& get_velocity() const override;
	std::vector<double>& get_previous_velocity() override;
	const std::vector<double>& get_previous_velocity() const override;
	
	// ========================================================================
	// Activity status
	// ========================================================================
	bool get_is_active() const override;
	void set_is_active(bool active) override;
	
	// ========================================================================
	// Substrate interaction rates (direct pointer access to SoA)
	// ========================================================================
	double* get_secretion_rates() override;
	const double* get_secretion_rates() const override;
	double* get_saturation_densities() override;
	const double* get_saturation_densities() const override;
	double* get_uptake_rates() override;
	const double* get_uptake_rates() const override;
	double* get_net_export_rates() override;
	const double* get_net_export_rates() const override;
	
	// ========================================================================
	// Internalized substrates
	// ========================================================================
	double* get_internalized_total_substrates() override;
	const double* get_internalized_total_substrates() const override;
	double* get_fraction_released_at_death() override;
	const double* get_fraction_released_at_death() const override;
	double* get_fraction_transferred_when_ingested() override;
	const double* get_fraction_transferred_when_ingested() const override;
	void release_internalized_substrates() override;
	void set_internal_uptake_constants(double dt) override;

	// ========================================================================
	// Secretion and uptake simulation (NOT SUPPORTED - physicore uses centralized solver)
	// ========================================================================
	/**
	 * @brief NOT SUPPORTED - physicore uses centralized solver architecture
	 * @throws std::runtime_error Always throws - use physicore solver instead
	 */
	void simulate_secretion_and_uptake(double dt) override;

	// ========================================================================
	// Voxel access
	// ========================================================================
	int get_current_voxel_index() override;
	
	// ========================================================================
	// Density and gradient access
	// ========================================================================
	double* nearest_density_vector() override;
	std::vector<double>& nearest_gradient(int substrate_index) override;
	std::vector<std::vector<double>>& nearest_gradient_vector() override;
};

} // namespace BioFVM

#endif // __physicore_basic_agent_adapter_h__
