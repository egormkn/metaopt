/*
 * ModelAddOn.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef MODELADDON_H_
#define MODELADDON_H_

#include "metaopt/model/scip/ScipModel.h"
#include "Solution.h"

namespace metaopt {

/**
 * ModelAddOnS are used to add variables that are needed for certain formulation to the basic ScipModel.
 *
 * Since implementations of Heuristics have to set values for all variables,
 * each add on has to supply a method computeSolutionVals that tries to compute solution values
 * for the additional variables from the basic variables of the ScipModel or other add ons.
 */
class ModelAddOn {
public:
	ModelAddOn(ScipModelPtr model);
	virtual ~ModelAddOn();

	/**
	 * Tries to compute solution values for the given solution.
	 *
	 * @returns true, if successfully computes solution values.
	 */
	virtual bool computeSolutionVals(SolutionPtr sol);

protected:
	ScipModelPtr _model;
};

} /* namespace metaopt */
#endif /* MODELADDON_H_ */
