#include "PhysiMeSS_implementation.h"
#include "PhysiMeSS_environment.h"
#include "PhysiMeSS_cell.h"
#include "PhysiMeSS_fibre.h"

using namespace PhysiCell;

PhysiMeSS_mechanics_implementation* PhysiMeSS_mechanics_implementation::instance = nullptr;

PhysiMeSS_mechanics_implementation* PhysiMeSS_mechanics_implementation::get_instance()
{
    if (!instance)
    {
        instance = new PhysiMeSS_mechanics_implementation();
    }
    return instance;
}

Mechanics_Environment_Interface* PhysiMeSS_mechanics_implementation::get_mechanics_environment()
{
    return &physimess_environment;
}

Mechanics_Agent_Interface* PhysiMeSS_mechanics_implementation::create_mechanics_agent(Cell* pCell)
{
    // the actual agent type is allocated in the Cell subclass (PhysiMeSS_Cell or PhysiMeSS_Fibre) constructor
    // lets just allocate a generic PhysiMeSS_Agent here to satisfy the interface
    return new PhysiMeSS_Agent(pCell);
}
