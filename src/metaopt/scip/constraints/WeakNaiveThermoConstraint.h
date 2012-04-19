/*
 * WeakNaiveThermoConstraint.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef WEAKNAIVETHERMOCONSTRAINT_H_
#define WEAKNAIVETHERMOCONSTRAINT_H_

#include "metaopt/model/scip/addon/PotentialDifferences.h"
#include "metaopt/model/scip/addon/ReactionDirections.h"

namespace metaopt {

/**
 * Creates the naive implementation of the weak formulation of thermodynamic constraints.
 * The weak formulation states, that a flux through a reaction is possible, if the reaction has negative potential difference.
 *
 * In particular, a reaction can have non-zero potential difference but no flow.
 * This property is the difference to the strong formulation of thermodynamic constraints.
 *
 * This constraint requires ReactionDirections as addon.
 * It basically simply initializes the constraints of the reaction direction vars of all internal reactions
 */
void createWeakNaiveThermoConstraint(ScipModelPtr model, ReactionDirectionsPtr dirs);

} /* namespace metaopt */
#endif /* WEAKNAIVETHERMOCONSTRAINT_H_ */
