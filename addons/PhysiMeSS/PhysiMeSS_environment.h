#ifndef __PhysiMeSS_environment_h__
#define __PhysiMeSS_environment_h__

#include "../../mechanics/PhysiCell_mechanics_environment.h"

/**
 * @brief Specialized mechanics environment for PhysiMeSS simulations.
 *
 * Overrides compute_velocities() to integrate PhysiMeSS-specific cell-fibre
 * and fibre-fibre interaction forces.
 */
class PhysiMeSS_Environment : public PhysiCell::mechanics_environment
{
public:
    PhysiMeSS_Environment() = default;
    virtual ~PhysiMeSS_Environment() = default;

    /// Compute forces for all agents including PhysiMeSS cell-fibre interactions.
    void compute_velocities(double dt) override;
};

/// Global PhysiMeSS mechanics environment instance.
extern PhysiMeSS_Environment physimess_environment;

#endif // __PhysiMeSS_environment_h__
