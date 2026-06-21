#include "PhysiCell_mechanics_legacy_implementation.h"

#include "PhysiCell_mechanics_agent.h"
#include "PhysiCell_mechanics_environment.h"

namespace PhysiCell {

Mechanics_Environment_Interface* PhysiCell_mechanics_legacy_implementation::get_mechanics_environment() {
    return &mech_environment;
}

Mechanics_Agent_Interface* PhysiCell_mechanics_legacy_implementation::create_mechanics_agent(BioFVM::Basic_Agent_Interface* pBasicAgent, Cell* pCell) {
    return new Mechanics_Agent(pBasicAgent, pCell);
}

} // namespace PhysiCell
