#ifndef __physicore_microenvironment_adapter_h__
#define __physicore_microenvironment_adapter_h__

#include <biofvm/microenvironment.h>
#include <memory>

#include "../../BioFVM/BioFVM_microenvironment_interface.h"
#include "basic_agent_adapter.h"
#include "mesh_adapter.h"

namespace BioFVM {

/**
 * @brief Adapter implementing BioFVM::Microenvironment_Interface using physicore backend
 *
 * This adapter wraps physicore::biofvm::microenvironment to provide compatibility with
 * the BioFVM interface. The physicore backend uses an immutable, builder-based architecture,
 * which differs fundamentally from BioFVM's mutable runtime configuration.
 */
class physicore_microenvironment_adapter : public Microenvironment_Interface
{
	friend class physicore_basic_agent_adapter;
private:
	std::unique_ptr<physicore::biofvm::microenvironment> me;
	
	// Mesh wrapper for BioFVM compatibility
	std::unique_ptr<physicore_mesh_wrapper> mesh_wrapper;

	// Agent container wrapper
	std::unique_ptr<Agent_Container> agent_wrapper;
	
    std::vector<std::vector<std::vector<double>>> gradient_vectors;
	
public:
	virtual ~physicore_microenvironment_adapter() = default;

	physicore::biofvm::microenvironment* get_physicore_microenvironment();

	// ========================================================================
	// Units access
	// ========================================================================
	const std::string& get_time_units() const override;
	const std::string& get_spatial_units() const override;

	// ========================================================================
	// Query methods - Size information
	// ========================================================================
	unsigned int number_of_densities() const override;
	unsigned int number_of_voxels() const override;
	unsigned int number_of_voxel_faces() const override;

	// ========================================================================
	// Substrate/Density management
	// ========================================================================
	int find_density_index(const std::string& name) const override;

	// ========================================================================
	// Voxel/Position access
	// ========================================================================
	int voxel_index(int i, int j, int k) const override;
	std::vector<unsigned int> cartesian_indices(int n) const override;
	int nearest_voxel_index(const std::vector<double>& position) const override;
	std::vector<unsigned int> nearest_cartesian_indices(const std::vector<double>& position) const override;
	Voxel& voxels(int voxel_index) override;
	const Voxel& voxels(int voxel_index) const override;
	Voxel& nearest_voxel(const std::vector<double>& position) override;

	// ========================================================================
	// Density vector access
	// ========================================================================
	double* density_vector(int n) override;
	double* density_vector(int i, int j) override;
	double* density_vector(int i, int j, int k) override;
	double* nearest_density_vector(const std::vector<double>& position) override;
	double* nearest_density_vector(int voxel_index) override;

	// ========================================================================
	// Gradient computation and access
	// ========================================================================
	void compute_all_gradient_vectors() override;
	void reset_all_gradient_vectors() override;
	std::vector<std::vector<double>>& gradient_vector(int n) override;
	std::vector<std::vector<double>>& gradient_vector(int i, int j) override;
	std::vector<std::vector<double>>& gradient_vector(int i, int j, int k) override;
	std::vector<std::vector<double>>& nearest_gradient_vector(const std::vector<double>& position) override;

	// ========================================================================
	// Simulation methods
	// ========================================================================
	void simulate_time_step(double dt) override;
	void simulate_diffusion_decay(double dt) override;
	void simulate_bulk_sources_and_sinks(double dt) override;
	void simulate_cell_sources_and_sinks(double dt) override;

	// ========================================================================
	// Dirichlet boundary conditions
	// ========================================================================
	void set_substrate_dirichlet_activation(int substrate_index, bool new_value) override;
	void set_substrate_dirichlet_activation(int substrate_index, int index, bool new_value) override;
	void set_substrate_dirichlet_activation(int index, std::vector<bool>& new_value) override;
	bool get_substrate_dirichlet_activation(int substrate_index) const override;
	bool get_substrate_dirichlet_activation(int substrate_index, int index) const override;
	double get_substrate_dirichlet_value(int substrate_index, int index) const override;
	bool& is_dirichlet_node(int voxel_index) override;

	// ========================================================================
	// Mesh access
	// ========================================================================
	const Cartesian_Mesh& get_mesh() const override;

	// ========================================================================
	// Agent container access
	// ========================================================================
	Agent_Container* get_agent_container() override;
	const Agent_Container* get_agent_container() const override;
	void set_agent_container(Agent_Container* container) override;

	// ========================================================================
	// Metadata access
	// ========================================================================
	const std::vector<std::string>& get_density_names() const override;
	const std::vector<std::string>& get_density_units() const override;
	const double* get_diffusion_coefficients() const override;
	const double* get_decay_rates() const override;

	// ========================================================================
	// Name access
	// ========================================================================
	const std::string& get_name() const override;

	// ========================================================================
	// Display and I/O
	// ========================================================================
	void display_information(std::ostream& os) const override;
	void write_to_matlab(std::string filename) override;

	// ========================================================================
	// Configuration query methods
	// ========================================================================
	bool simulate_2D() const override;
	bool calculate_gradients() const override;

	bool setup_microenvironment_from_XML(const std::string& filename) override;
	void initialize() override;
};

} // namespace BioFVM

#endif // __physicore_microenvironment_adapter_h__
