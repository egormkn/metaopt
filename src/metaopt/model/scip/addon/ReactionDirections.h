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

struct ReactionDirVars {
	SCIP_VAR* dir;
	SCIP_VAR* slack_flux_fwd;
	SCIP_VAR* slack_flux_bwd;
	SCIP_VAR* slack_pot_fwd;
	SCIP_VAR* slack_pot_bwd;

	ReactionDirVars() {
		dir = NULL;
		slack_flux_fwd = NULL;
		slack_flux_bwd = NULL;
		slack_pot_fwd = NULL;
		slack_pot_bwd = NULL;
	}
};

/**
 * ReactionDirections couple flux directions to potential differences.
 *
 * However, this link is realized by relaxing the strict inequalities that are used usually.
 * In particular, if a potential difference is zero, nonnegative flux in any direction will still be allowed.
 *
 * To deal with this problem, please create a WeakThermoConstraint. (TODO: implement!)
 */
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

	/**
	 * fixes the decision variable.
	 */
	virtual void setDirection(SCIP_NODE* node, ReactionPtr rxn, bool fwd);

	virtual boost::shared_ptr<const boost::unordered_set<ReactionPtr> > retainFixedDirections(boost::shared_ptr<const boost::unordered_set<ReactionPtr> > candidates);

	/**
	 * Helper method for ScipModel destruction process.
	 * This method should only be called from the destructor of ScipModel.
	 * It processes the destruction of this ModelAddOn that still requires data of ScipModel.
	 * A reference to ScipModel is supplied, since weak_ptrs are already expired during destruction.
	 * After this method has been called, the indicator destroyed is set.
	 *
	 * Important: If this addon has dependencies on other addons, these dependencies must be destroyed, too.
	 */
	virtual void destroy(ScipModel* model);

private:
	PotentialDifferencesPtr _potDiff;
	boost::unordered_map<ReactionPtr, ReactionDirVars> _dirs;
};

typedef boost::shared_ptr<ReactionDirections> ReactionDirectionsPtr;

inline ReactionDirectionsPtr createReactionDirections( ScipModelPtr model, PotentialDifferencesPtr potDiff) {
	return ReactionDirectionsPtr(new ReactionDirections(model, potDiff));
}

} /* namespace metaopt */
#endif /* REACTIONDIRECTIONS_H_ */
