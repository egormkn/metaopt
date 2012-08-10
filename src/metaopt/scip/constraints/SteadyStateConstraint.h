/*
 * SteadyStateConstraint.h
 *
 *  Created on: 11.04.2012
 *      Author: arnem
 */

#ifndef STEADYSTATECONSTRAINT_H_
#define STEADYSTATECONSTRAINT_H_

#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Creates linear constraints that enforce the steady state assumption on the flux variables.
 */
void createSteadyStateConstraint(ScipModelPtr smodel);

} /* namespace metaopt */
#endif /* STEADYSTATECONSTRAINT_H_ */
