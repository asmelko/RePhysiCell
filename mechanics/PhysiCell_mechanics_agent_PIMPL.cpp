#include "PhysiCell_mechanics_agent_PIMPL.h"

namespace PhysiCell {

Mechanics_Agent_PIMPL::Mechanics_Agent_PIMPL()
	: pImpl(nullptr)
{
}

Mechanics_Agent_PIMPL::Mechanics_Agent_PIMPL(Mechanics_Agent_Interface* impl)
	: pImpl(impl)
{
}

Mechanics_Agent_PIMPL::~Mechanics_Agent_PIMPL()
{
	delete pImpl;
	pImpl = nullptr;
}

} // namespace PhysiCell
