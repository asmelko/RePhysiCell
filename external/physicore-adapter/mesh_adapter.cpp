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

#include "mesh_adapter.h"

namespace BioFVM {

physicore_mesh_wrapper::physicore_mesh_wrapper(const physicore::biofvm::cartesian_mesh& mesh)
	: physicore_mesh(mesh)
{
	// Initialize BioFVM Cartesian_Mesh members based on physicore mesh
	
	// Set coordinates
	const auto& grid = physicore_mesh.grid_shape;
	const auto& bbox_min = physicore_mesh.bounding_box_mins;
	const auto& bbox_max = physicore_mesh.bounding_box_maxs;
	const auto& voxel_shape = physicore_mesh.voxel_shape;
	
	// Calculate spacing
	dx = static_cast<double>(voxel_shape[0]);
	dy = static_cast<double>(voxel_shape[1]);
	dz = static_cast<double>(voxel_shape[2]);
	
	// Set volumes and surface areas
	dV = dx * dy * dz;
	dS_xy = dx * dy;
	dS_yz = dy * dz;
	dS_xz = dx * dz;
	dS = 2.0 * (dS_xy + dS_yz + dS_xz); // Total surface area
	
	// Set bounding box
	bounding_box.resize(6);
	bounding_box[0] = static_cast<double>(bbox_min[0]);
	bounding_box[1] = static_cast<double>(bbox_min[1]);
	bounding_box[2] = static_cast<double>(bbox_min[2]);
	bounding_box[3] = static_cast<double>(bbox_max[0]);
	bounding_box[4] = static_cast<double>(bbox_max[1]);
	bounding_box[5] = static_cast<double>(bbox_max[2]);
	
	// Build coordinate vectors
	x_coordinates.clear();
	y_coordinates.clear();
	z_coordinates.clear();
	
	for (size_t i = 0; i < grid[0]; ++i) {
		x_coordinates.push_back(bbox_min[0] + (i + 0.5) * dx);
	}
	for (size_t j = 0; j < grid[1]; ++j) {
		y_coordinates.push_back(bbox_min[1] + (j + 0.5) * dy);
	}
	for (size_t k = 0; k < grid[2]; ++k) {
		z_coordinates.push_back(bbox_min[2] + (k + 0.5) * dz);
	}
	
	// Set mesh flags
	Cartesian_mesh = true;
	uniform_mesh = (dx == dy && dy == dz);
	regular_mesh = true;
	use_voxel_faces = false;
	
	// Initialize voxels
	size_t num_voxels = physicore_mesh.voxel_count();
	voxels.resize(num_voxels);
	
	for (size_t idx = 0; idx < num_voxels; ++idx) {
		voxels[idx].mesh_index = static_cast<int>(idx);
		voxels[idx].volume = dV;
		voxels[idx].is_Dirichlet = false; // Will be set based on Dirichlet conditions
		
		// Calculate voxel center
		auto indices = cartesian_indices(static_cast<unsigned int>(idx));
		voxels[idx].center.resize(3);
		voxels[idx].center[0] = x_coordinates[indices[0]];
		voxels[idx].center[1] = y_coordinates[indices[1]];
		voxels[idx].center[2] = z_coordinates[indices[2]];
	}
}

void physicore_mesh_wrapper::update_dirichlet_flags(const physicore::biofvm::microenvironment& me) {
	// Mark Dirichlet interior voxels
	for (size_t i = 0; i < me.dirichlet_interior_voxels_count; ++i) {
		auto indices =  std::span(me.dirichlet_interior_voxels.get() + i * me.mesh.dims, me.mesh.dims);
		auto voxel_idx = physicore_mesh.linearize(indices[0], indices[1], indices[2]);
		if (voxel_idx < voxels.size()) {
			voxels[voxel_idx].is_Dirichlet = true;
		}
	}
}

} // namespace BioFVM
