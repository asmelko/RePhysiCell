#ifndef __BioFVM_position_entity_h__
#define __BioFVM_position_entity_h__

#include <vector>

namespace BioFVM {

/**
 * @brief Minimal spatial entity that owns a 3D position vector.
 *
 * Position_Entity is a plain, non-polymorphic struct intended to be the
 * single source of truth for an agent's spatial location.  It carries no
 * virtual methods so it can be embedded by value with zero overhead.
 *
 * Ownership semantics
 * -------------------
 *   - Basic_Agent and Mechanics_Agent hold a *non-owning* Position_Entity*
 *     that normally points to their own embedded `default_position` member.
 *   - Cell inherits from Position_Entity.  Its constructor rebinds both
 *     agents to point at the Cell's own Position_Entity subobject, so every
 *     read/write goes to the same memory with no extra indirection.
 */
struct Position_Entity
{
	std::vector<double> position; // always size 3: {x, y, z}

	Position_Entity() : position(3, 0.0) {}
};

} // namespace BioFVM

#endif
