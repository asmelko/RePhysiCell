#ifndef __physicore_mesh_adapter_h__
#define __physicore_mesh_adapter_h__

#include <biofvm/microenvironment.h>
#include "../../BioFVM/BioFVM_mesh.h"

namespace BioFVM {

/**
 * @brief Adapter wrapping physicore's cartesian_mesh as BioFVM's Cartesian_Mesh
 *
 * This wrapper creates a BioFVM-compatible Cartesian_Mesh interface around
 * physicore's cartesian_mesh. It initializes coordinate vectors, voxels,
 * and mesh properties based on the physicore mesh configuration.
 */
class physicore_mesh_wrapper : public Cartesian_Mesh
{
private:
	const physicore::biofvm::cartesian_mesh& physicore_mesh;
	
public:
	/**
	 * @brief Construct mesh wrapper from physicore cartesian_mesh
	 * @param mesh Reference to physicore mesh (must outlive this wrapper)
	 */
	explicit physicore_mesh_wrapper(const physicore::biofvm::cartesian_mesh& mesh);
	
	/**
	 * @brief Update Dirichlet flags on voxels based on microenvironment configuration
	 * @param me Reference to physicore microenvironment with Dirichlet node data
	 */
	void update_dirichlet_flags(const physicore::biofvm::microenvironment& me);
};

} // namespace BioFVM

#endif // __physicore_mesh_adapter_h__
