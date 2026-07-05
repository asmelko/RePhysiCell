#ifndef __PhysiMeSS_implementation_h__
#define __PhysiMeSS_implementation_h__

#include "../../mechanics/PhysiCell_mechanics_implementation.h"

class PhysiMeSS_Environment;

/**
 * @brief Mechanics implementation factory for PhysiMeSS.
 *
 * Instantiates PhysiMeSS-specific environment and agent types.
 */
class PhysiMeSS_mechanics_implementation : public PhysiCell::Mechanics_implementation
{
public:
    static PhysiMeSS_mechanics_implementation* get_instance();

    PhysiCell::Mechanics_Environment_Interface* get_mechanics_environment() override;
    PhysiCell::Mechanics_Agent_Interface* create_mechanics_agent(PhysiCell::Cell* pCell) override;

private:
    PhysiMeSS_mechanics_implementation() = default;
    static PhysiMeSS_mechanics_implementation* instance;
};

#endif // __PhysiMeSS_implementation_h__
