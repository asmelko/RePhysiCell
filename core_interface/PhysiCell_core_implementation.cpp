#include "../core/PhysiCell_cell.h"
#include "PhysiCell_core_implementation.h"

#include "../BioFVM/BioFVM_agent_container.h"
#include "../BioFVM/BioFVM_basic_agent.h"
#include "../BioFVM/BioFVM_microenvironment.h"
#include "PhysiCell_legacy_core_entity.h"
#include "../mechanics/PhysiCell_mechanics_environment.h"

namespace PhysiCell {

CoreImplementation* CoreImplementation::instance_ = nullptr;

CoreImplementation* CoreImplementation::get_instance()
{
	return instance_;
}

void CoreImplementation::set_instance(CoreImplementation* implementation)
{
	instance_ = implementation;
}

BioFVM::Microenvironment_Interface* LegacyCoreImplementation::get_microenvironment()
{
	return &BioFVM::microenvironment;
}

Mechanics_Environment_Interface* LegacyCoreImplementation::get_mechanics_environment()
{
	return &mech_environment;
}

CoreEntity* LegacyCoreImplementation::create_core_entity(Cell* pCell)
{
	return new LegacyCoreEntity(pCell);
}

std::vector<Cell*>* LegacyCoreImplementation::get_all_cells()
{
	return &cells_;
}

Mechanics_Environment_Interface* get_mechanics_environment_i()
{
	if (CoreImplementation::get_instance() == nullptr)
	{
		return nullptr;
	}
	return CoreImplementation::get_instance()->get_mechanics_environment();
}

} // namespace PhysiCell

namespace BioFVM {

Microenvironment_Interface* get_microenvironment_i()
{
	if (PhysiCell::CoreImplementation::get_instance() == nullptr)
	{
		return nullptr;
	}
	return PhysiCell::CoreImplementation::get_instance()->get_microenvironment();
}

} // namespace BioFVM
