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
