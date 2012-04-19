/*
 * PotentialDifferences.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef POTENTIALDIFFERENCES_H_
#define POTENTIALDIFFERENCES_H_

#include "metaopt/model/Model.h"
#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/model/scip/ModelAddOn.h"
#include <boost/unordered_map.hpp>

namespace metaopt {

/**
 * This class stores the variables of potential differences.
 *
 * Potential differences are not needed for all formulations and is only used to model the problem.
 * Hence, these variables are separated from the others in ScipModel.
 */
class PotentialDifferences : public ModelAddOn {
public:
	PotentialDifferences(ScipModelPtr model);
	virtual ~PotentialDifferences();

	/**
	 * Fetch the variable for the potential difference of the specified reaction.
	 *
	 * In principle potential difference variables can be created for all reactions,
	 * however it does not make much sense to compute potential differences for exchange reactions.
	 *
	 * The corresponding SCIP_VAR is only created if it is needed (i.e. by a call to this method).
	 *
	 * The returned SCIP_VAR lives as long as these PotentialDifferences live.
	 */
	SCIP_VAR* getPotDiff(ReactionPtr rxn);

	/**
	 * Returns the potential difference of the specified reaction.
	 * If the problem has been solved, the value of an optimal solution is returned.
	 * Else, the value of the current LP-relaxation is returned. The LP-relaxation, however has to be solved to optimality.
	 *
	 * Precondition: ScipModel::hasCurrentFlux
	 *
	 * If the model has no current flux, a PreconditionViolatedException is thrown
	 */

	double getCurrentPotDiff(ReactionPtr rxn);

	/**
	 * Tests if the corresponding variable has already been created.
	 */
	bool hasPotDiff(ReactionPtr rxn);

	bool computeSolutionVals(SolutionPtr sol) const;

private:
	boost::unordered_map<ReactionPtr, SCIP_VAR*> _potDiff;
};

typedef boost::shared_ptr<PotentialDifferences> PotentialDifferencesPtr;

inline PotentialDifferencesPtr createPotentialDifferences(ScipModelPtr model) {
	return PotentialDifferencesPtr(new PotentialDifferences(model));
}

} /* namespace metaopt */
#endif /* POTENTIALDIFFERENCES_H_ */
