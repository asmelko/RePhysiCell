#ifndef __PhysiCell_core_entity_interface_h__
#define __PhysiCell_core_entity_interface_h__

#include <vector>

#include "../BioFVM/BioFVM_basic_agent_interface.h"
#include "../mechanics/PhysiCell_mechanics_agent_interface.h"

namespace PhysiCell
{

/**
 * @brief Everything a simulation core provides for a single cell: the BioFVM agent,
 * the mechanics agent and the single source of truth for the cell's position.
 */
class CoreEntity :  public BioFVM::Basic_Agent_Interface, public Mechanics_Agent_Interface
{
public:
    /// The single source of truth for the position of the entity, see Basic_Agent_Interface::position().
    double* position() override = 0;

    /**
     * @brief The object implementing the mechanics interface, for implementation-specific downcasts.
     * By default it is the entity itself.
     */
    virtual Mechanics_Agent_Interface* get_mechanics_implementation()
    { return this; }
};

}

#endif
