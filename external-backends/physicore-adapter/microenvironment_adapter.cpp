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

#include "microenvironment_adapter.h"
#include "biofvm/microenvironment.h"
#include "mesh_adapter.h"
#include "agent_container_adapter.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace BioFVM {

// ============================================================================
// Adapter implementation - Constructors
// ============================================================================

physicore_microenvironment_adapter::physicore_microenvironment_adapter(
	std::unique_ptr<physicore::biofvm::microenvironment> environment)
	: me(std::move(environment))
{
	if (!me) {
		throw std::runtime_error("physicore_microenvironment_adapter: null microenvironment provided");
	}
	
	mesh_wrapper = std::make_unique<physicore_mesh_wrapper>(me->mesh);
	mesh_wrapper->update_dirichlet_flags(*me);
	
	if (me->agents) {
		agent_wrapper = std::make_unique<physicore_agent_container_wrapper>(me->agents);
	}
}

physicore_microenvironment_adapter::physicore_microenvironment_adapter(const std::string& config_path)
	: physicore_microenvironment_adapter(
		physicore::biofvm::microenvironment::create_from_config(config_path))
{
}

// ============================================================================
// Helper methods
// ============================================================================

void physicore_microenvironment_adapter::invalidate_caches() {
	gradient_cache.clear();
}

void physicore_microenvironment_adapter::compute_gradient_at_voxel(
	int voxel_index, int substrate_index, std::vector<double>& gradient) const 
{
	gradient.resize(3, 0.0);
	
	auto indices = cartesian_indices(voxel_index);
	int i = indices[0];
	int j = indices[1];
	int k = indices[2];
	
	const auto& grid = me->mesh.grid_shape;
	const auto& voxel_shape = me->mesh.voxel_shape;
	
	// Central finite difference for interior voxels, one-sided for boundaries
	
	// X-direction gradient
	if (i > 0 && i < static_cast<int>(grid[0]) - 1) {
		double rho_plus = me->solver->get_substrate_density(substrate_index, i + 1, j, k);
		double rho_minus = me->solver->get_substrate_density(substrate_index, i - 1, j, k);
		gradient[0] = (rho_plus - rho_minus) / (2.0 * voxel_shape[0]);
	} else if (i == 0 && grid[0] > 1) {
		double rho_0 = me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_1 = me->solver->get_substrate_density(substrate_index, i + 1, j, k);
		gradient[0] = (rho_1 - rho_0) / voxel_shape[0];
	} else if (i == static_cast<int>(grid[0]) - 1 && grid[0] > 1) {
		double rho_0 = me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_m1 = me->solver->get_substrate_density(substrate_index, i - 1, j, k);
		gradient[0] = (rho_0 - rho_m1) / voxel_shape[0];
	}
	
	// Y-direction gradient
	if (j > 0 && j < static_cast<int>(grid[1]) - 1) {
		double rho_plus = me->solver->get_substrate_density(substrate_index, i, j + 1, k);
		double rho_minus = me->solver->get_substrate_density(substrate_index, i, j - 1, k);
		gradient[1] = (rho_plus - rho_minus) / (2.0 * voxel_shape[1]);
	} else if (j == 0 && grid[1] > 1) {
		double rho_0 = me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_1 = me->solver->get_substrate_density(substrate_index, i, j + 1, k);
		gradient[1] = (rho_1 - rho_0) / voxel_shape[1];
	} else if (j == static_cast<int>(grid[1]) - 1 && grid[1] > 1) {
		double rho_0 = me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_m1 = me->solver->get_substrate_density(substrate_index, i, j - 1, k);
		gradient[1] = (rho_0 - rho_m1) / voxel_shape[1];
	}
	
	// Z-direction gradient
	if (k > 0 && k < static_cast<int>(grid[2]) - 1) {
		double rho_plus = me->solver->get_substrate_density(substrate_index, i, j, k + 1);
		double rho_minus = me->solver->get_substrate_density(substrate_index, i, j, k - 1);
		gradient[2] = (rho_plus - rho_minus) / (2.0 * voxel_shape[2]);
	} else if (k == 0 && grid[2] > 1) {
		double rho_0 = me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_1 = me->solver->get_substrate_density(substrate_index, i, j, k + 1);
		gradient[2] = (rho_1 - rho_0) / voxel_shape[2];
	} else if (k == static_cast<int>(grid[2]) - 1 && grid[2] > 1) {
		double rho_0 = me->solver->get_substrate_density(substrate_index, i, j, k);
		double rho_m1 = me->solver->get_substrate_density(substrate_index, i, j, k - 1);
		gradient[2] = (rho_0 - rho_m1) / voxel_shape[2];
	}
}

// ============================================================================
// Units access
// ============================================================================

std::string& physicore_microenvironment_adapter::get_time_units() {
	return me->time_units;
}

const std::string& physicore_microenvironment_adapter::get_time_units() const {
	return me->time_units;
}

std::string& physicore_microenvironment_adapter::get_spatial_units() {
	return me->space_units;
}

const std::string& physicore_microenvironment_adapter::get_spatial_units() const {
	return me->space_units;
}

// ============================================================================
// Query methods - Size information
// ============================================================================

unsigned int physicore_microenvironment_adapter::number_of_densities() const {
	return static_cast<unsigned int>(me->substrates_count);
}

unsigned int physicore_microenvironment_adapter::number_of_voxels() const {
	return static_cast<unsigned int>(me->mesh.voxel_count());
}

unsigned int physicore_microenvironment_adapter::number_of_voxel_faces() const {
	// Physicore doesn't use voxel faces
	return 0;
}

// ============================================================================
// Substrate/Density management
// ============================================================================

int physicore_microenvironment_adapter::find_density_index(const std::string& name) const {
	auto it = std::find(me->substrates_names.begin(), me->substrates_names.end(), name);
	if (it != me->substrates_names.end()) {
		return static_cast<int>(std::distance(me->substrates_names.begin(), it));
	}
	return -1;
}

// void physicore_microenvironment_adapter::add_density() {
// 	throw std::runtime_error("add_density: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify substrate configuration.");
// }

// void physicore_microenvironment_adapter::add_density(const std::string& name, const std::string& units) {
// 	throw std::runtime_error("add_density: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify substrate configuration.");
// }

// void physicore_microenvironment_adapter::add_density(const std::string& name, const std::string& units,
//                                                       double diffusion_constant, double decay_rate) {
// 	throw std::runtime_error("add_density: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify substrate configuration.");
// }

// void physicore_microenvironment_adapter::set_density(int index, const std::string& name, const std::string& units) {
// 	throw std::runtime_error("set_density: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify substrate configuration.");
// }

// void physicore_microenvironment_adapter::set_density(int index, const std::string& name, const std::string& units,
//                                                       double diffusion_constant, double decay_rate) {
// 	throw std::runtime_error("set_density: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify substrate configuration.");
// }

// void physicore_microenvironment_adapter::resize_densities(int new_size) {
// 	throw std::runtime_error("resize_densities: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify substrate configuration.");
// }

// ============================================================================
// Voxel/Position access
// ============================================================================

int physicore_microenvironment_adapter::voxel_index(int i, int j, int k) const {
	return static_cast<int>(me->mesh.linearize(
		static_cast<physicore::index_t>(i),
		static_cast<physicore::index_t>(j),
		static_cast<physicore::index_t>(k)
	));
}

std::vector<unsigned int> physicore_microenvironment_adapter::cartesian_indices(int n) const {
	const auto& grid = me->mesh.grid_shape;
	std::vector<unsigned int> indices(3);
	
	size_t idx = static_cast<size_t>(n);
	indices[2] = static_cast<unsigned int>(idx / (grid[0] * grid[1]));
	idx %= (grid[0] * grid[1]);
	indices[1] = static_cast<unsigned int>(idx / grid[0]);
	indices[0] = static_cast<unsigned int>(idx % grid[0]);
	
	return indices;
}

int physicore_microenvironment_adapter::nearest_voxel_index(const std::vector<double>& position) const {
	std::vector<physicore::real_t> pos(3);
	pos[0] = position.size() > 0 ? position[0] : 0.0;
	pos[1] = position.size() > 1 ? position[1] : 0.0;
	pos[2] = position.size() > 2 ? position[2] : 0.0;
	
	auto voxel_pos = me->mesh.voxel_position(pos);
	return static_cast<int>(me->mesh.linearize(voxel_pos[0], voxel_pos[1], voxel_pos[2]));
}

std::vector<unsigned int> physicore_microenvironment_adapter::nearest_cartesian_indices(
	const std::vector<double>& position) const 
{
	int idx = nearest_voxel_index(position);
	return cartesian_indices(idx);
}

Voxel& physicore_microenvironment_adapter::voxels(int voxel_index) {
	return mesh_wrapper->voxels[voxel_index];
}

const Voxel& physicore_microenvironment_adapter::voxels(int voxel_index) const {
	return mesh_wrapper->voxels[voxel_index];
}

Voxel& physicore_microenvironment_adapter::nearest_voxel(const std::vector<double>& position) {
	int idx = nearest_voxel_index(position);
	return voxels(idx);
}

// ============================================================================
// Density vector access
// ============================================================================

double* physicore_microenvironment_adapter::density_vector(int n) {
	auto indices = cartesian_indices(n);
    return &me->get_substrate_density(0, indices[0], indices[1], indices[2]);
}

double* physicore_microenvironment_adapter::density_vector(int i, int j) {
    return &me->get_substrate_density(0, i, j, 0);
}

double* physicore_microenvironment_adapter::density_vector(int i, int j, int k) {
    return &me->get_substrate_density(0, i, j, k);
}

double* physicore_microenvironment_adapter::nearest_density_vector(const std::vector<double>& position) {
    auto indices = me->mesh.voxel_position(position);
    return density_vector(indices[0], indices[1], indices[2]);
}

double* physicore_microenvironment_adapter::nearest_density_vector(int voxel_index) {
	return density_vector(voxel_index);
}

const double* physicore_microenvironment_adapter::density_vector(int n) const {
	return density_vector(n);
}

// ============================================================================
// Gradient computation and access
// ============================================================================

std::vector<std::vector<double>> physicore_microenvironment_adapter::compute_gradient_vector(int n) {
	std::vector<std::vector<double>> gradients(me->substrates_count);
	
	for (size_t s = 0; s < me->substrates_count; ++s) {
		compute_gradient_at_voxel(n, s, gradients[s]);
	}
	
    return gradients;
}

void physicore_microenvironment_adapter::compute_all_gradient_vectors() {
	// We are using lazy evaluation for gradients
}

void physicore_microenvironment_adapter::reset_all_gradient_vectors() {
    std::lock_guard<std::mutex> lock(gradient_cache_mutex);
	gradient_cache.clear();
}

std::vector<std::vector<double>>& physicore_microenvironment_adapter::gradient_vector(int n) {
    {
        std::lock_guard<std::mutex> lock(gradient_cache_mutex);
	    auto it = gradient_cache.find(n);

        if (it != gradient_cache.end()) {
            return it->second;
        }
    }

	auto gradiens = compute_gradient_vector(n);

    {
        std::lock_guard<std::mutex> lock(gradient_cache_mutex);
        gradient_cache[n] = std::move(gradiens);
        return gradient_cache[n];
    }
}

std::vector<std::vector<double>>& physicore_microenvironment_adapter::gradient_vector(int i, int j) {
	int idx = voxel_index(i, j, 0);
	return gradient_vector(idx);
}

std::vector<std::vector<double>>& physicore_microenvironment_adapter::gradient_vector(int i, int j, int k) {
	int idx = voxel_index(i, j, k);
	return gradient_vector(idx);
}

std::vector<std::vector<double>>& physicore_microenvironment_adapter::nearest_gradient_vector(
	const std::vector<double>& position) 
{
	int idx = nearest_voxel_index(position);
	return gradient_vector(idx);
}

// ============================================================================
// Simulation methods
// ============================================================================

void physicore_microenvironment_adapter::simulate_diffusion_decay(double dt) {
	// Physicore handles diffusion-decay in run_single_timestep
	// Calculate how many iterations needed for dt
	if (me->diffusion_timestep > 0) {
		physicore::index_t iterations = static_cast<physicore::index_t>(
			std::ceil(dt / me->diffusion_timestep)
		);
		me->solver->solve(*me, iterations);
	}
	
	invalidate_caches();
}

void physicore_microenvironment_adapter::simulate_bulk_sources_and_sinks(double dt) {
	// Physicore handles bulk sources via bulk_fnc during solve
	// This is a no-op since it's integrated into the solver
}

void physicore_microenvironment_adapter::simulate_cell_sources_and_sinks(double dt) {
	// Physicore handles cell sources via agent container during solve
	// This is a no-op since it's integrated into the solver
}

// ============================================================================
// Dirichlet boundary conditions
// ============================================================================

// void physicore_microenvironment_adapter::add_dirichlet_node(int voxel_index, std::vector<double>& value) {
// 	throw std::runtime_error("add_dirichlet_node: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify Dirichlet conditions.");
// }

// void physicore_microenvironment_adapter::update_dirichlet_node(int voxel_index, std::vector<double>& new_value) {
// 	throw std::runtime_error("update_dirichlet_node: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify Dirichlet conditions.");
// }

// void physicore_microenvironment_adapter::update_dirichlet_node(int voxel_index, int substrate_index, double new_value) {
// 	throw std::runtime_error("update_dirichlet_node: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify Dirichlet conditions.");
// }

// void physicore_microenvironment_adapter::remove_dirichlet_node(int voxel_index) {
// 	throw std::runtime_error("remove_dirichlet_node: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify Dirichlet conditions.");
// }

// void physicore_microenvironment_adapter::apply_dirichlet_conditions() {
// 	// Physicore applies Dirichlet conditions automatically during solve
// 	// This is a no-op
// }

void physicore_microenvironment_adapter::set_substrate_dirichlet_activation(int substrate_index, bool new_value) {
	throw std::runtime_error("Deprecated function. Use set_substrate_dirichlet_activation with voxel index.");
}

void physicore_microenvironment_adapter::set_substrate_dirichlet_activation(int substrate_index, int index, bool new_value) {
	auto indices = cartesian_indices(index);
    auto value = get_substrate_dirichlet_value(substrate_index, index);
    me->update_dirichlet_interior_voxel({indices[0], indices[1], indices[2]}, substrate_index, value, new_value);
    mesh_wrapper->update_dirichlet_flags(*me);
    me->update_dirichlet_conditions();
}

void physicore_microenvironment_adapter::set_substrate_dirichlet_activation(int index, std::vector<bool>& new_value) {
	throw std::runtime_error("set_substrate_dirichlet_activation: physicore microenvironment is immutable after construction. "
	                         "Rebuild using create_from_config() to modify Dirichlet conditions.");
}

bool physicore_microenvironment_adapter::get_substrate_dirichlet_activation(int substrate_index) const {
	throw std::runtime_error("Deprecated function. Use get_substrate_dirichlet_activation with voxel index.");
}

bool physicore_microenvironment_adapter::get_substrate_dirichlet_activation(int substrate_index, int index) const {
	// Check if specific voxel has Dirichlet condition for this substrate
    auto indices = cartesian_indices(index);
	// Search interior Dirichlet nodes
	for (size_t i = 0; i < me->dirichlet_interior_voxels_count; ++i) {
        if (me->mesh.dims >= 1 && me->dirichlet_interior_voxels[i * me->mesh.dims + 1] != indices[0]) continue;
        if (me->mesh.dims >= 2 && me->dirichlet_interior_voxels[i * me->mesh.dims + 2] != indices[1]) continue;
        if (me->mesh.dims >= 3 && me->dirichlet_interior_voxels[i * me->mesh.dims + 3] != indices[2]) continue;
		
        size_t offset = i * me->substrates_count + substrate_index;
        return me->dirichlet_interior_conditions && me->dirichlet_interior_conditions[offset];
	}
	
	return false;
}

double physicore_microenvironment_adapter::get_substrate_dirichlet_value(int substrate_index, int index) const {
	// Search interior Dirichlet nodes
    auto indices = cartesian_indices(index);

	for (size_t i = 0; i < me->dirichlet_interior_voxels_count; ++i) {
        if (me->mesh.dims >= 1 && me->dirichlet_interior_voxels[i * me->mesh.dims + 1] != indices[0]) continue;
        if (me->mesh.dims >= 2 && me->dirichlet_interior_voxels[i * me->mesh.dims + 2] != indices[1]) continue;
        if (me->mesh.dims >= 3 && me->dirichlet_interior_voxels[i * me->mesh.dims + 3] != indices[2]) continue;

        size_t offset = i * me->substrates_count + substrate_index;
        if (me->dirichlet_interior_values) {
            return me->dirichlet_interior_values[offset];
        }
	}
	
	return 0.0;
}

bool& physicore_microenvironment_adapter::is_dirichlet_node(int voxel_index) {
	return mesh_wrapper->voxels[voxel_index].is_Dirichlet;
}

// ============================================================================
// Mesh access
// ============================================================================

const Cartesian_Mesh& physicore_microenvironment_adapter::get_mesh() const {
	return *mesh_wrapper;
}

// Cartesian_Mesh& physicore_microenvironment_adapter::get_mesh() {
// 	return *mesh_wrapper;
// }

// ============================================================================
// Agent container access
// ============================================================================

Agent_Container_Interface* physicore_microenvironment_adapter::get_agent_container() {
	return agent_wrapper.get();
}

const Agent_Container_Interface* physicore_microenvironment_adapter::get_agent_container() const {
	return agent_wrapper.get();
}

void physicore_microenvironment_adapter::set_agent_container(Agent_Container_Interface* container) {
	throw std::runtime_error("set_agent_container: Cannot replace physicore agent container. "
	                         "Agent container is managed by physicore microenvironment.");
}

// ============================================================================
// Metadata access
// ============================================================================

// std::vector<std::string>& physicore_microenvironment_adapter::get_density_names() {
// 	return me->substrates_names;
// }

const std::vector<std::string>& physicore_microenvironment_adapter::get_density_names() const {
	return me->substrates_names;
}

// std::vector<std::string>& physicore_microenvironment_adapter::get_density_units() {
// 	return me->substrates_units;
// }

const std::vector<std::string>& physicore_microenvironment_adapter::get_density_units() const {
	return me->substrates_units;
}

// std::vector<double>& physicore_microenvironment_adapter::get_diffusion_coefficients() {
// 	ensure_metadata_cache();
// 	return diffusion_coeffs_cache;
// }

const double* physicore_microenvironment_adapter::get_diffusion_coefficients() const {
	return me->diffusion_coefficients.get();
}

// std::vector<double>& physicore_microenvironment_adapter::get_decay_rates() {
// 	ensure_metadata_cache();
// 	return decay_rates_cache;
// }

const double* physicore_microenvironment_adapter::get_decay_rates() const {
	return me->decay_rates.get();
}

// ============================================================================
// Name access
// ============================================================================

// std::string& physicore_microenvironment_adapter::get_name() {
// 	return me->name;
// }

const std::string& physicore_microenvironment_adapter::get_name() const {
	return me->name;
}

// ============================================================================
// Display and I/O
// ============================================================================

void physicore_microenvironment_adapter::display_information(std::ostream& os) const {
	me->print_info(os);
}

void physicore_microenvironment_adapter::write_to_matlab(std::string filename) {
	
    int number_of_data_entries = mesh_wrapper->voxels.size();
	int size_of_each_datum = 3 + 1 + number_of_densities() ; // x,y,z, volume, densities

	FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "multiscale_microenvironment" );  

	// storing data as cols 
	for( int i=0; i < number_of_data_entries ; i++ )
	{
		fwrite( (char*) &( mesh_wrapper->voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
		fwrite( (char*) &( mesh_wrapper->voxels[i].center[1] ) , sizeof(double) , 1 , fp );   
		fwrite( (char*) &( mesh_wrapper->voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
		fwrite( (char*) &( mesh_wrapper->voxels[i].volume ) , sizeof(double) , 1 , fp ); 

		// densities  

		for( unsigned int j=0 ; j < number_of_densities() ; j++)
		{ fwrite( (char*) &( density_vector(i)[j] ) , sizeof(double) , 1 , fp ); }
	}

	fclose( fp ); 
	return;
}

// ============================================================================
// Spatial setup methods
// ============================================================================

// void physicore_microenvironment_adapter::resize_space(int x_nodes, int y_nodes, int z_nodes) {
// 	throw std::runtime_error("resize_space: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify spatial configuration.");
// }

// void physicore_microenvironment_adapter::resize_space(double x_start, double x_end, double y_start, double y_end,
//                                                        double z_start, double z_end, int x_nodes, int y_nodes, int z_nodes) {
// 	throw std::runtime_error("resize_space: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify spatial configuration.");
// }

// void physicore_microenvironment_adapter::resize_space(double x_start, double x_end, double y_start, double y_end,
//                                                        double z_start, double z_end, double dx_new, double dy_new, double dz_new) {
// 	throw std::runtime_error("resize_space: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify spatial configuration.");
// }

// void physicore_microenvironment_adapter::resize_space_uniform(double x_start, double x_end, double y_start, double y_end,
//                                                                double z_start, double z_end, double dx_new) {
// 	throw std::runtime_error("resize_space_uniform: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify spatial configuration.");
// }

// void physicore_microenvironment_adapter::resize_voxels(int new_number_of_voxels) {
// 	throw std::runtime_error("resize_voxels: physicore microenvironment is immutable after construction. "
// 	                         "Rebuild using create_from_config() to modify spatial configuration.");
// }

// ============================================================================
// Update methods
// ============================================================================

void physicore_microenvironment_adapter::update_rates() {
	// Physicore doesn't have a separate update_rates step
	// Rates are configured at build time
}

// ============================================================================
// Configuration query methods
// ============================================================================

bool physicore_microenvironment_adapter::simulate_2D() const {
	return me->mesh.dims == 2;
}

bool physicore_microenvironment_adapter::calculate_gradients() const {
	return false;
}

bool physicore_microenvironment_adapter::setup_microenvironment_from_XML(const std::string& filename) {
    try {
        me = physicore::biofvm::microenvironment::create_from_config(filename);
    } catch (const std::exception& e) {
        return false;
    }
    return true;
}

} // namespace BioFVM