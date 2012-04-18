/*
 * PotentialDifferences.cpp
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#include "PotentialDifferences.h"

namespace metaopt {

PotentialDifferences::PotentialDifferences(ScipModelPtr model) : ModelAddOn(model) {
	// nothing to do
}

PotentialDifferences::~PotentialDifferences() {
	// nothing to do
}

SCIP_VAR* PotentialDifferences::getPotDiff(ReactionPtr rxn) {

}

bool PotentialDifferences::computeSolutionVals(SolutionPtr sol) {

}

} /* namespace metaopt */
