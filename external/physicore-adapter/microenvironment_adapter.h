/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2025, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#ifndef __physicore_microenvironment_adapter_h__
#define __physicore_microenvironment_adapter_h__

#include <biofvm/microenvironment.h>
#include <map>
#include <memory>
#include <mutex>

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
 *
 * @note PERFORMANCE AND MUTABILITY CONSTRAINTS:
 * - The physicore microenvironment is IMMUTABLE after construction
 * - Methods that would mutate the environment (add_density, resize_space, etc.) will
 *   throw std::runtime_error
 * - To modify configuration, rebuild the microenvironment using create_from_config()
 * - Gradients are computed manually via finite differences and cached for performance
 * - Density vector access creates temporary views - avoid repeated calls in tight loops
 *
 * @note XML CONFIGURATION:
 * - Use setup_microenvironment_from_XML() or create_from_config() for initialization
 * - The adapter translates PhysiCell XML format to physicore's configuration format
 */
class physicore_microenvironment_adapter : public Microenvironment_Interface
{
	friend class physicore_basic_agent_adapter;
private:
	std::unique_ptr<physicore::biofvm::microenvironment> me;
	
	// Mesh wrapper for BioFVM compatibility
	std::unique_ptr<physicore_mesh_wrapper> mesh_wrapper;
	
	// Agent container wrapper
	std::unique_ptr<Agent_Container_Interface> agent_wrapper;
	
	// Cached data for interface compatibility
	std::mutex gradient_cache_mutex;
    std::map<int, std::vector<std::vector<double>>> gradient_cache;
	
	// Helper methods
	void invalidate_caches();
	std::vector<std::vector<double>> compute_gradient_vector(int n);
	void compute_gradient_at_voxel(int voxel_index, int substrate_index, std::vector<double>& gradient) const;

public:
	// /**
	//  * @brief Construct adapter from existing physicore microenvironment
	//  * @param environment Unique pointer to physicore microenvironment (ownership transferred)
	//  */
	// explicit physicore_microenvironment_adapter(std::unique_ptr<physicore::biofvm::microenvironment> environment);
	
	// /**
	//  * @brief Construct adapter from XML configuration file
	//  * @param config_path Path to physicore configuration file
	//  */
	// explicit physicore_microenvironment_adapter(const std::string& config_path);
	
	virtual ~physicore_microenvironment_adapter() = default;

	physicore::biofvm::microenvironment* get_physicore_microenvironment();

	// ========================================================================
	// Units access
	// ========================================================================
	std::string& get_time_units() override;
	const std::string& get_time_units() const override;
	std::string& get_spatial_units() override;
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
	// void add_density() override;
	// void add_density(const std::string& name, const std::string& units) override;
	// void add_density(const std::string& name, const std::string& units,
	//                  double diffusion_constant, double decay_rate) override;
	// void set_density(int index, const std::string& name, const std::string& units) override;
	// void set_density(int index, const std::string& name, const std::string& units,
	//                  double diffusion_constant, double decay_rate) override;
	// void resize_densities(int new_size) override;

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
	const double* density_vector(int n) const override;

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
	// void add_dirichlet_node(int voxel_index, std::vector<double>& value) override;
	// void update_dirichlet_node(int voxel_index, std::vector<double>& new_value) override;
	// void update_dirichlet_node(int voxel_index, int substrate_index, double new_value) override;
	// void remove_dirichlet_node(int voxel_index) override;
	// void apply_dirichlet_conditions() override;
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
	// Cartesian_Mesh& get_mesh() override;

	// ========================================================================
	// Agent container access
	// ========================================================================
	Agent_Container_Interface* get_agent_container() override;
	const Agent_Container_Interface* get_agent_container() const override;
	void set_agent_container(Agent_Container_Interface* container) override;

	// ========================================================================
	// Metadata access
	// ========================================================================
	// std::vector<std::string>& get_density_names() override;
	const std::vector<std::string>& get_density_names() const override;
	// std::vector<std::string>& get_density_units() override;
	const std::vector<std::string>& get_density_units() const override;
	// std::vector<double>& get_diffusion_coefficients() override;
	const double* get_diffusion_coefficients() const override;
	// std::vector<double>& get_decay_rates() override;
	const double* get_decay_rates() const override;

	// ========================================================================
	// Name access
	// ========================================================================
	// std::string& get_name() override;
	const std::string& get_name() const override;

	// ========================================================================
	// Display and I/O
	// ========================================================================
	void display_information(std::ostream& os) const override;
	void write_to_matlab(std::string filename) override;

	// ========================================================================
	// Spatial setup methods
	// ========================================================================
	// void resize_space(int x_nodes, int y_nodes, int z_nodes) override;
	// void resize_space(double x_start, double x_end, double y_start, double y_end,
	//                   double z_start, double z_end, int x_nodes, int y_nodes, int z_nodes) override;
	// void resize_space(double x_start, double x_end, double y_start, double y_end,
	//                   double z_start, double z_end, double dx_new, double dy_new, double dz_new) override;
	// void resize_space_uniform(double x_start, double x_end, double y_start, double y_end,
	//                           double z_start, double z_end, double dx_new) override;
	// void resize_voxels(int new_number_of_voxels) override;

	// ========================================================================
	// Update methods
	// ========================================================================
	void update_rates() override;

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
