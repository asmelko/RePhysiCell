#include "PhysiCell_legacy_core_entity.h"

#include "../BioFVM/BioFVM_basic_agent.h"
#include "../mechanics/PhysiCell_mechanics_agent.h"

namespace PhysiCell {

LegacyCoreEntity::LegacyCoreEntity(Cell* pCell)
	: LegacyCoreEntity(std::make_unique<Mechanics_Agent>(pCell))
{}

LegacyCoreEntity::LegacyCoreEntity(std::unique_ptr<Mechanics_Agent> mechanics_agent)
	: basic_agent_(std::make_unique<BioFVM::Basic_Agent>())
{
	basic_agent_->bind_position(position_.data());
	set_mechanics_agent(std::move(mechanics_agent));
}

LegacyCoreEntity::~LegacyCoreEntity() = default;

Mechanics_Agent_Interface* LegacyCoreEntity::get_mechanics_implementation()
{
	return mechanics_agent_.get();
}

void LegacyCoreEntity::set_mechanics_agent(std::unique_ptr<Mechanics_Agent> mechanics_agent)
{
	mechanics_agent_ = std::move(mechanics_agent);
	mechanics_agent_->bind_position(position_.data());
}

} // namespace PhysiCell
