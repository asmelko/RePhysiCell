#ifndef __PhysiMeSS_fibre_h__
#define __PhysiMeSS_fibre_h__

#include "PhysiMeSS_agent.h"

#include "../../core/PhysiCell_cell.h"

class PhysiMeSS_CellAgent;

bool isFibre(PhysiCell::Cell* pCell);

bool isFibre(PhysiCell::Cell_Definition* cellDef);

std::vector<PhysiCell::Cell_Definition*>* getFibreCellDefinitions();

/**
 * @brief Mechanics-implementation class for PhysiMeSS fibres.
 *
 * Specializes PhysiMeSS_Agent for fibre-type agents: tracks fibre geometry,
 * crosslinks, and fibre-cell/fibre-fibre interaction forces.
 */
class PhysiMeSS_FibreAgent : public PhysiMeSS_Agent
{
public:
    std::vector<PhysiMeSS_FibreAgent*> fibres_crosslinkers;
    std::vector<double> fibres_crosslink_point;

    double mLength;
    double mRadius;
    int X_crosslink_count;
    int fail_count;

    explicit PhysiMeSS_FibreAgent(PhysiCell::Cell* pCell);
    virtual ~PhysiMeSS_FibreAgent() {};
    void assign_fibre_orientation() ;

    void check_out_of_bounds(std::vector<double>& position);
    void add_potentials_from_fibre(PhysiMeSS_FibreAgent* other_fibre);
    void add_potentials_from_cell(PhysiMeSS_CellAgent* cell);

    void register_fibre_voxels() override;
    void deregister_fibre_voxels() override;

    std::vector<double> nearest_point_on_fibre(std::vector<double> point, std::vector<double>& displacement);

    void check_fibre_crosslinks(PhysiMeSS_FibreAgent* fibre_neighbor);
    void add_crosslinks();

    bool is_fibre() const override { return true; }
};

/**
 * @brief Thin wrapper for PhysiMeSS fibres.
 *
 * A PhysiCell::Cell subclass that swaps its mechanics pImpl to PhysiMeSS_FibreAgent,
 * providing access to extended PhysiMeSS neighbor lists and multi-voxel fibre tracking.
 */
class PhysiMeSS_Fibre : public PhysiCell::Cell
{
public:
    explicit PhysiMeSS_Fibre();
    virtual ~PhysiMeSS_Fibre() = default;

    /// Returns the PhysiMeSS-specific mechanics implementation.
    PhysiMeSS_FibreAgent* get_physimess_agent() const;
};

#endif