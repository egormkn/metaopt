/*
 * WeakNaiveThermoConstraint.cpp
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#include "WeakNaiveThermoConstraint.h"

namespace metaopt {

void createWeakNaiveThermoConstraint(ScipModelPtr model, ReactionDirectionsPtr dirs) {
	ModelPtr m = model->getModel();
	foreach(ReactionPtr rxn, m->getReactions()) {
		if(!rxn->isExchange()) {
			// create the var, even if we don't do anything with it now.
			// this also initializes the corresponding constraints
			dirs->getDirection(rxn);
		}
	}
}


} /* namespace metaopt */
