#include "PhysiMeSS_implementation.h"
#include "PhysiMeSS_environment.h"
#include "PhysiMeSS_cell.h"
#include "PhysiMeSS_fibre.h"
#include "../../core_interface/PhysiCell_legacy_core_entity.h"

using namespace PhysiCell;

PhysiMeSS_core_implementation* PhysiMeSS_core_implementation::instance = nullptr;

PhysiMeSS_core_implementation* PhysiMeSS_core_implementation::get_instance()
{
    if (!instance)
    {
        instance = new PhysiMeSS_core_implementation();
    }
    return instance;
}

Mechanics_Environment_Interface* PhysiMeSS_core_implementation::get_mechanics_environment()
{
    return &physimess_environment;
}

CoreEntity* PhysiMeSS_core_implementation::create_core_entity(Cell* pCell)
{
    // the actual agent type is allocated in the Cell subclass (PhysiMeSS_Cell or PhysiMeSS_Fibre) constructor
    // lets just allocate a generic PhysiMeSS_Agent here to satisfy the interface
    return new LegacyCoreEntity(std::make_unique<PhysiMeSS_Agent>(pCell));
}
