/*
 * ReactionDirections.h
 *
 *  Created on: 19.04.2012
 *      Author: arnem
 */

#ifndef REACTIONDIRECTIONS_H_
#define REACTIONDIRECTIONS_H_

#include <boost/unordered_map.hpp>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "metaopt/scip/ScipError.h"
#include "PotentialDifferences.h"
#include "metaopt/model/scip/ModelAddOn.h"

namespace metaopt {

class ReactionDirections : public ModelAddOn {
public:
	ReactionDirections(ScipModelPtr model, PotentialDifferencesPtr potDiff);
	virtual ~ReactionDirections();

	/**
	 * Returns the direction of the specified reaction.
	 * 1 is forward, 0 is backward. It may also return fractional values.
	 * If the problem has been solved, the value of an optimal solution is returned.
	 * Else, the value of the current LP-relaxation is returned. The LP-relaxation, however has to be solved to optimality.
	 *
	 * Precondition: ScipModel::hasCurrentFlux
	 *
	 * If the model has no current flux, a PreconditionViolatedException is thrown
	 */

	double getCurrentDirection(ReactionPtr rxn);

	/**
	 * Tests if the corresponding variable has already been created.
	 */
	bool hasDirection(ReactionPtr rxn);

	/**
	 * Fetch the variable for the reaction direction of the specified reaction.
	 *
	 * In principle direction variables can be created for all reactions,
	 * however it does not make much sense to compute directions for exchange reactions.
	 *
	 * The corresponding SCIP_VAR is only created if it is needed (i.e. by a call to this method).
	 *
	 * The returned SCIP_VAR lives as long as these ReactionDirections live.
	 */
	SCIP_VAR* getDirection(ReactionPtr rxn);

	bool computeSolutionVals(SolutionPtr sol) const;

private:
	PotentialDifferencesPtr _potDiff;
	boost::unordered_map<ReactionPtr, SCIP_VAR*> _dirs;
};

typedef boost::shared_ptr<ReactionDirections> ReactionDirectionsPtr;

inline ReactionDirectionsPtr createReactionDirections( ScipModelPtr model, PotentialDifferencesPtr potDiff) {
	return ReactionDirectionsPtr(new ReactionDirections(model, potDiff));
}

} /* namespace metaopt */
#endif /* REACTIONDIRECTIONS_H_ */
