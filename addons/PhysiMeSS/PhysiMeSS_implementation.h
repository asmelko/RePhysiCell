#ifndef __PhysiMeSS_implementation_h__
#define __PhysiMeSS_implementation_h__

#include "../../core_interface/PhysiCell_core_implementation.h"

class PhysiMeSS_Environment;

/**
 * @brief Core implementation factory for PhysiMeSS.
 *
 * Reuses the legacy BioFVM parts and instantiates PhysiMeSS-specific
 * mechanics environment and agent types.
 */
class PhysiMeSS_core_implementation : public PhysiCell::LegacyCoreImplementation
{
public:
    static PhysiMeSS_core_implementation* get_instance();

    PhysiCell::Mechanics_Environment_Interface* get_mechanics_environment() override;
    PhysiCell::CoreEntity* create_core_entity(PhysiCell::Cell* pCell) override;

private:
    PhysiMeSS_core_implementation() = default;
    static PhysiMeSS_core_implementation* instance;
};

#endif // __PhysiMeSS_implementation_h__
