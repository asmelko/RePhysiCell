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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

#ifndef __PhysiCell_constants_h__
#define __PhysiCell_constants_h__

#include <string>
#include <unordered_map>

namespace PhysiCell
{
	
namespace Constants
{
	constexpr double pi=3.1415926535897932384626433832795;
	
	constexpr double cell_removal_threshold_volume = 20; // 20 cubic microns -- about 1% of typical cell 
	constexpr int keep_pushed_out_cells_in_outer_voxel=1;
	constexpr int solid_boundary = 2;
	constexpr int default_boundary_condition_for_pushed_out_agents = keep_pushed_out_cells_in_outer_voxel;		
	
	constexpr int deterministic_necrosis = 0;
	constexpr int stochastic_necrosis = 1;
	
	constexpr int mesh_min_x_index=0;
	constexpr int mesh_min_y_index=1;
	constexpr int mesh_min_z_index=2;
	constexpr int mesh_max_x_index=3;
	constexpr int mesh_max_y_index=4;
	constexpr int mesh_max_z_index=5;			
	
	constexpr int mesh_lx_face_index=0;
	constexpr int mesh_ly_face_index=1;
	constexpr int mesh_lz_face_index=2;
	constexpr int mesh_ux_face_index=3;
	constexpr int mesh_uy_face_index=4;
	constexpr int mesh_uz_face_index=5;
	
	// currently recognized cell cycle models 
	constexpr int advanced_Ki67_cycle_model= 0;
	constexpr int basic_Ki67_cycle_model=1;
	constexpr int flow_cytometry_cycle_model=2;
	constexpr int live_apoptotic_cycle_model=3;
	constexpr int total_cells_cycle_model=4;
	constexpr int live_cells_cycle_model = 5; 
	constexpr int flow_cytometry_separated_cycle_model = 6; 
	constexpr int cycling_quiescent_model = 7; 
	
	// currently recognized death models 
	constexpr int apoptosis_death_model = 100; 
	constexpr int necrosis_death_model = 101; 
	constexpr int autophagy_death_model = 102; 
	
	constexpr int custom_cycle_model=9999; 
	
	// currently recognized cell cycle and death phases 
	// cycle phases
	constexpr int Ki67_positive_premitotic=0; 
	constexpr int Ki67_positive_postmitotic=1; 
	constexpr int Ki67_positive=2; 
	constexpr int Ki67_negative=3; 
	constexpr int G0G1_phase=4;
	constexpr int G0_phase=5;
	constexpr int G1_phase=6; 
	constexpr int G1a_phase=7; 
	constexpr int G1b_phase=8;
	constexpr int G1c_phase=9;
	constexpr int S_phase=10;
	constexpr int G2M_phase=11;
	constexpr int G2_phase=12;
	constexpr int M_phase=13;
	constexpr int live=14;
	
	constexpr int G1pm_phase = 15;
	constexpr int G1ps_phase = 16; 
	
	constexpr int cycling = 17; 
	constexpr int quiescent = 18; 
	
	
	constexpr int custom_phase = 9999;
	// death phases
	constexpr int apoptotic=100;
	constexpr int necrotic_swelling=101;
	constexpr int necrotic_lysed=102;
	constexpr int necrotic=103; 
	constexpr int debris=104; 
};

extern std::string time_units;
extern std::string space_units;
extern double diffusion_dt; 
extern double mechanics_dt;
extern double phenotype_dt;
extern double intracellular_dt;


extern std::unordered_map<std::string,int> cycle_model_codes;
int find_cycle_model_code( std::string model_name ); 

};

#endif
