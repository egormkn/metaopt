/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    fast-tfva - efficient thermodynamic constrained flux variability analysis.
    Copyright (C) 2012  Arne Müller, arne.mueller@fu-berlin.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * PotentialDifferences.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef POTENTIALDIFFERENCES_H_
#define POTENTIALDIFFERENCES_H_

#include "model/Model.h"
#include "model/scip/ScipModel.h"
#include "model/scip/ModelAddOn.h"
#include <boost/unordered_map.hpp>
#include "Properties.h"

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

	/**
	 * Helper method for ScipModel destruction process.
	 * This method should only be called from the destructor of ScipModel.
	 * It processes the destruction of this ModelAddOn that still requires data of ScipModel.
	 * A reference to ScipModel is supplied, since boost::weak_ptrs are already expired during destruction.
	 * After this method has been called, the indicator destroyed is set.
	 *
	 * Important: If this addon has dependencies on other addons, these dependencies must be destroyed, too.
	 */
	virtual void destroy(ScipModel* model);

private:
	boost::unordered_map<ReactionPtr, SCIP_VAR*> _potDiff;
};

typedef boost::shared_ptr<PotentialDifferences> PotentialDifferencesPtr;

inline PotentialDifferencesPtr createPotentialDifferences(ScipModelPtr model) {
	return PotentialDifferencesPtr(new PotentialDifferences(model));
}

} /* namespace metaopt */
#endif /* POTENTIALDIFFERENCES_H_ */
