/*
 * SteadyStateConstraint.h
 *
 *  Created on: 11.04.2012
 *      Author: arnem
 */

#ifndef STEADYSTATECONSTRAINT_H_
#define STEADYSTATECONSTRAINT_H_

#include "metaopt/model/scip/ScipModel.h"

namespace metaopt {

void createSteadyStateConstraint(ScipModelPtr model);

} /* namespace metaopt */
#endif /* STEADYSTATECONSTRAINT_H_ */
