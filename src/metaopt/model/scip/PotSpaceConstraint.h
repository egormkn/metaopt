/*
 * PotSpaceConstraint.h
 *
 *  Created on: 16.07.2012
 *      Author: arnem
 */

#ifndef POTSPACECONSTRAINT_H_
#define POTSPACECONSTRAINT_H_

#include <boost/unordered_map.hpp>
#include "scip/scip.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Use this constraint to formulate side constraints on the potential space variables, even if we do not create the potential variables explicitly.
 * Formulates _coef * mu >= 0
 */
struct PotSpaceConstraint {
	boost::unordered_map<MetabolitePtr, double> _coef;
};

typedef boost::shared_ptr<PotSpaceConstraint> PotSpaceConstraintPtr;


} /* namespace metaopt */
#endif /* POTSPACECONSTRAINT_H_ */
