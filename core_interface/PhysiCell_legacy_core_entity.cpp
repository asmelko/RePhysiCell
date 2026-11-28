#include "PhysiCell_legacy_core_entity.h"

#include "../BioFVM/BioFVM_basic_agent.h"
#include "../mechanics/PhysiCell_mechanics_agent.h"

#include <stdexcept>

namespace PhysiCell {

LegacyCoreEntity::LegacyCoreEntity(Cell* pCell)
	: basic_agent_(std::make_unique<BioFVM::Basic_Agent>(position_)),
	  mechanics_agent_(std::make_unique<Mechanics_Agent>(pCell, position_))
{}

LegacyCoreEntity::LegacyCoreEntity(std::function<std::unique_ptr<Mechanics_Agent>(std::vector<double>&)> mechanics_factory)
	: basic_agent_(std::make_unique<BioFVM::Basic_Agent>(position_))
{
	set_mechanics_agent(mechanics_factory(position_));
}

LegacyCoreEntity::~LegacyCoreEntity() = default;

Mechanics_Agent_Interface* LegacyCoreEntity::get_mechanics_implementation()
{
	return mechanics_agent_.get();
}

void LegacyCoreEntity::set_mechanics_agent(std::unique_ptr<Mechanics_Agent> mechanics_agent)
{
	if (!mechanics_agent || !mechanics_agent->uses_position(position_))
	{
		throw std::invalid_argument("Mechanics agent must use the core entity's position storage");
	}
	mechanics_agent_ = std::move(mechanics_agent);
}

} // namespace PhysiCell
