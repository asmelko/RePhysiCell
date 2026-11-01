#ifndef __PhysiCell_core_implementation_h__
#define __PhysiCell_core_implementation_h__

#include <vector>

#include "PhysiCell_core_entity_interface.h"
#include "../BioFVM/BioFVM_microenvironment_interface.h"
#include "../mechanics/PhysiCell_mechanics_environment_interface.h"

namespace PhysiCell {

class Cell;

/**
 * @brief Abstract factory for the whole simulation core (BioFVM + mechanics).
 *
 * A concrete subclass is installed once at simulation startup via set_instance().
 */
class CoreImplementation
{
	static CoreImplementation* instance_;

public:
	virtual ~CoreImplementation() = default;

	/** @brief Return the single active microenvironment. */
	virtual BioFVM::Microenvironment_Interface* get_microenvironment() = 0;

	/** @brief Return the single active mechanics environment. */
	virtual Mechanics_Environment_Interface* get_mechanics_environment() = 0;

	/** @brief Create the core entity backing a new cell; the cell takes ownership. */
	virtual CoreEntity* create_core_entity(Cell* pCell) = 0;

	/** @brief Return the implementation-owned list of all cells. */
	virtual std::vector<Cell*>* get_all_cells() = 0;

	/** @brief Return the globally installed implementation, or nullptr. */
	static CoreImplementation* get_instance();

	/** @brief Install the implementation used for this simulation. */
	static void set_instance(CoreImplementation* implementation);
};

/** @brief The original (legacy) BioFVM and mechanics implementation. */
class LegacyCoreImplementation : public CoreImplementation
{
public:
	BioFVM::Microenvironment_Interface* get_microenvironment() override;
	Mechanics_Environment_Interface* get_mechanics_environment() override;
	CoreEntity* create_core_entity(Cell* pCell) override;
	std::vector<Cell*>* get_all_cells() override;

private:
	std::vector<Cell*> cells_;
};

/** @brief Convenience accessor for the mechanics environment. */
Mechanics_Environment_Interface* get_mechanics_environment_i();

} // namespace PhysiCell

#endif // __PhysiCell_core_implementation_h__
