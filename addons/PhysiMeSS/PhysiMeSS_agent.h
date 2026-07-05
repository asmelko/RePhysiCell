#ifndef __PhysiMeSS_agent_h__
#define __PhysiMeSS_agent_h__

#include <list>
#include <vector>

#include "../../mechanics/PhysiCell_mechanics_agent.h"

namespace PhysiCell { class Cell; }

/**
 * @brief Base mechanics-implementation class for all PhysiMeSS agents.
 *
 * Inherits from Mechanics_Agent (the pImpl side of the new mechanics interface)
 * and adds the PhysiMeSS extended neighbor list and multi-voxel fibre tracking.
 * The owning Cell* is accessible via get_cell() through the inherited pOwner pointer.
 */
class PhysiMeSS_Agent : public PhysiCell::Mechanics_Agent
{
    private:
    public:
    
    std::vector<PhysiMeSS_Agent*> physimess_neighbors;
    std::list<int> physimess_voxels;

    explicit PhysiMeSS_Agent(PhysiCell::Cell* pCell);
    virtual ~PhysiMeSS_Agent(){};

    virtual void register_fibre_voxels() {};
    virtual void deregister_fibre_voxels() {};

    std::list<int> find_agent_voxels();
    void find_agent_neighbors();

    /// Returns the owning Cell (safe: Cell inherits from Mechanics_Agent_PIMPL = pOwner).
    PhysiCell::Cell* get_cell() const;

    virtual bool is_fibre() const { return false; }

};

#endif