#include "PhysiCell_core_entity_PIMPL.h"

namespace PhysiCell {

CoreEntityPIMPL::CoreEntityPIMPL(CoreEntity* core_entity)
	: pImpl(core_entity)
{}

CoreEntityPIMPL::~CoreEntityPIMPL()
{
	delete pImpl;
	pImpl = nullptr;
}

} // namespace PhysiCell
