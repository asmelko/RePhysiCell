#include "PhysiCell_mechanics_implementation.h"

namespace PhysiCell {

Mechanics_implementation* Mechanics_implementation::instance_ = nullptr;

Mechanics_implementation* Mechanics_implementation::get_instance()
{
	return instance_;
}

void Mechanics_implementation::set_instance(Mechanics_implementation* implementation)
{
	instance_ = implementation;
}

Mechanics_Environment_Interface* get_mechanics_environment_i()
{
	if (Mechanics_implementation::get_instance() == nullptr)
	{
		return nullptr;
	}
	return Mechanics_implementation::get_instance()->get_mechanics_environment();
}

} // namespace PhysiCell
