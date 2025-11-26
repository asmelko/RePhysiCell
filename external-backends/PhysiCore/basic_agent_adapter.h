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

#ifndef __physicore_basic_agent_adapter_h__
#define __physicore_basic_agent_adapter_h__

#include <biofvm/agent.h>
#include <cassert>
#include "../../BioFVM/BioFVM_basic_agent_interface.h"

namespace BioFVM {

/**
 * @brief Lightweight stack-based proxy adapter bridging BioFVM's Basic_Agent_Interface 
 *        with physicore's Structure-of-Arrays agent storage
 *
 * This adapter wraps an agent index and references to physicore's SoA data arrays,
 * providing BioFVM's mutable interface over physicore's immutable backend. It stores
 * BioFVM-specific state (velocity, ID, type) that doesn't exist in physicore.
 *
 * @note IMPORTANT USAGE CONSTRAINTS:
 * - This is a STACK-BASED PROXY - create on-demand, use immediately, don't store long-term
 * - Proxy lifetime must not exceed physicore agent container's lifetime
 * - SINGLE-THREADED USE ONLY - not thread-safe due to shared SoA access
 * - Position/velocity caches use lazy synchronization - call sync methods before physicore solver steps
 * - simulate_secretion_and_uptake() NOT SUPPORTED - physicore uses centralized solver (throws exception)
 *
 * @note SYNCHRONIZATION:
 * - Position modifications via assign_position() write directly to physicore SoA
 * - Position reads via get_position() return cached vector (synced lazily on first access)
 * - Velocity is adapter-only state (not in physicore) - use update_position() with caution
 * - Call sync_position_to_physicore() before physicore solver if position cache was modified
 *
 * Example usage:
 * @code
 * // Create proxy on stack for temporary use
 * physicore_basic_agent_adapter proxy(&agent_container, index);
 * proxy.register_microenvironment(&microenv);
 * 
 * // Use BioFVM interface
 * proxy.assign_position(10.0, 5.0, 0.0);
 * double vol = proxy.get_total_volume();
 * double* rates = proxy.get_secretion_rates();
 * 
 * // Proxy destroyed automatically when out of scope
 * @endcode
 */
class physicore_basic_agent_adapter : public Basic_Agent_Interface
{
private:
	// Reference to physicore agent container
	physicore::biofvm::agent* agent;
	physicore::index_t index;
	
	// BioFVM microenvironment link
	Microenvironment_Interface* microenvironment;
	
	// BioFVM-specific state not in physicore
	int ID;
	int type;
	bool is_active;
	int cached_voxel_index;
	
	// Cached vectors with lazy synchronization
	mutable std::vector<double> position_cache;
	mutable bool position_cache_dirty;
	std::vector<double> velocity;
	std::vector<double> previous_velocity;

public:
	/**
	 * @brief Construct agent proxy wrapping physicore agent at given index
	 * @param agent_ptr Pointer to physicore agent container
	 * @param index Agent index in physicore SoA arrays
	 * @throws std::out_of_range if index >= agent count
	 */
	physicore_basic_agent_adapter(
		physicore::biofvm::agent* agent_ptr,
		physicore::index_t index);
	
	virtual ~physicore_basic_agent_adapter() = default;

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
	Microenvironment_Interface* get_microenvironment() override;

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
	double* get_position() override;
	const double* get_position() const override;
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
	std::vector<double>& nearest_density_vector() override;
	std::vector<double>& nearest_gradient(int substrate_index) override;
	std::vector<std::vector<double>>& nearest_gradient_vector() override;

	// ========================================================================
	// Synchronization helpers
	// ========================================================================
	/**
	 * @brief Write cached position back to physicore SoA if modified
	 * 
	 * Call this before physicore solver steps if position was accessed via
	 * get_position() and potentially modified through the returned reference.
	 */
	void sync_position_to_physicore();
};

} // namespace BioFVM

#endif // __physicore_basic_agent_adapter_h__
