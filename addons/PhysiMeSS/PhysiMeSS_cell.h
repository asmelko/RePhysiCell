#ifndef __PhysiMeSS_cell_h__
#define __PhysiMeSS_cell_h__

#include "../../core/PhysiCell_cell.h"
#include "PhysiMeSS_agent.h"

class PhysiMeSS_FibreAgent;

/**
 * @brief Mechanics-implementation class for PhysiMeSS cells.
 *
 * Specializes PhysiMeSS_Agent for cell-type agents: tracks stuck/unstuck state
 * and implements cell-fibre interaction forces.
 */
class PhysiMeSS_CellAgent : public PhysiMeSS_Agent
{
    private:
    public:
    
    int stuck_counter;
    int unstuck_counter;

    explicit PhysiMeSS_CellAgent(PhysiCell::Cell* pCell);
    virtual ~PhysiMeSS_CellAgent() = default;

    bool is_fibre() const override { return false; }

    void add_potentials_from_fibre(PhysiMeSS_FibreAgent* fibre);

    void register_fibre_voxels() override;
    void deregister_fibre_voxels() override;

    void force_update_motility_vector(double dt_);

    virtual void degrade_fibre(PhysiMeSS_FibreAgent* pFibre);
};

/**
 * @brief Thin wrapper for PhysiMeSS cells.
 *
 * A PhysiCell::Cell subclass that swaps its mechanics pImpl to PhysiMeSS_CellAgent,
 * providing access to extended PhysiMeSS neighbor lists and multi-voxel tracking.
 */
class PhysiMeSS_Cell : public PhysiCell::Cell
{
public:
    explicit PhysiMeSS_Cell();
    virtual ~PhysiMeSS_Cell() = default;

    /// Returns the PhysiMeSS-specific mechanics implementation.
    PhysiMeSS_CellAgent* get_physimess_agent() const;
};

#endif