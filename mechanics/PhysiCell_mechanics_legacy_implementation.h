#ifndef PHYSICELL_MECHANICS_LEGACY_IMPLEMENTATION_H
#define PHYSICELL_MECHANICS_LEGACY_IMPLEMENTATION_H

#include "PhysiCell_mechanics_implementation.h"

namespace PhysiCell {

class PhysiCell_mechanics_legacy_implementation : public PhysiCell::Mechanics_implementation
{
public:
	Mechanics_Environment_Interface* get_mechanics_environment() override;

    Mechanics_Agent_Interface* create_mechanics_agent(BioFVM::Basic_Agent_Interface* pBasicAgent, Cell* pCell) override;
};

}

#endif // PHYSICELL_MECHANICS_LEGACY_IMPLEMENTATION_H
