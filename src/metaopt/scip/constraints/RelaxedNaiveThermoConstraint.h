/*
 * WeakNaiveThermoConstraint.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef RELAXEDNAIVETHERMOCONSTRAINT_H_
#define RELAXEDNAIVETHERMOCONSTRAINT_H_

#include "metaopt/model/scip/addon/PotentialDifferences.h"
#include "metaopt/model/scip/addon/ReactionDirections.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Creates a naive implementation of thermodynamic constraints where all strict inequalities are relaxed to weak inequalities.
 *
 * The relaxed constraint simply states that flux and potential difference must have different signs, or one of them has to be zero.
 *
 * Note, that this formulation is significantly weaker than even the weak formulation of thermodynamic constraints.
 *
 * This constraint requires ReactionDirections as addon.
 * It basically simply initializes the constraints of the reaction direction vars of all internal reactions
 */
void createRelaxedNaiveThermoConstraint(ScipModelPtr model, ReactionDirectionsPtr dirs);

} /* namespace metaopt */
#endif /* RELAXEDNAIVETHERMOCONSTRAINT_H_ */
