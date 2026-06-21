#ifndef __PhysiCell_mechanics_implementation_h__
#define __PhysiCell_mechanics_implementation_h__

#include "PhysiCell_mechanics_environment_interface.h"
#include "PhysiCell_mechanics_agent_interface.h"

namespace BioFVM {
class Basic_Agent_Interface;
}

namespace PhysiCell {

class Cell;

/**
 * @brief Abstract factory for the mechanics subsystem.
 *
 * Mirrors BioFVM::BioFVM_implementation. A concrete subclass is installed
 * once at simulation startup via set_instance(). After that point,
 * get_mechanics_environment_i() may be used to obtain the active
 * Mechanics_Environment_Interface from anywhere in the code.
 */
class Mechanics_implementation
{
	static Mechanics_implementation* instance_;

public:
	virtual ~Mechanics_implementation() = default;

	/** @brief Return the single active mechanics environment. */
	virtual Mechanics_Environment_Interface* get_mechanics_environment() = 0;

	/** @brief Create and register a new mechanics agent. */
	virtual Mechanics_Agent_Interface* create_mechanics_agent(BioFVM::Basic_Agent_Interface* pBasicAgent, Cell* pCell) = 0;

	/** @brief Return the list of all tracked mechanics agents. */
	// virtual std::vector<Mechanics_Agent_Interface*>* get_all_mechanics_agents() = 0;

	// ---- Singleton accessors -----------------------------------------------

	/** @brief Return the globally installed Mechanics_implementation, or nullptr. */
	static Mechanics_implementation* get_instance();

	/**
	 * @brief Install the concrete implementation to be used for this simulation.
	 *
	 * Must be called before any mechanics operation is performed.
	 */
	static void set_instance(Mechanics_implementation* implementation);
};

/** @brief Convenience free function mirroring BioFVM::get_microenvironment_i(). */
Mechanics_Environment_Interface* get_mechanics_environment_i();

} // namespace PhysiCell

#endif // __PhysiCell_mechanics_implementation_h__
