/*
 * FVA.cpp
 *
 *  Created on: 14.07.2012
 *      Author: arne
 */

#include "FVA.h"
#include "metaopt/Properties.h"

namespace metaopt {

void fva(ModelPtr model, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max ) {
	LPFluxPtr flux(new LPFlux(model, true));
	foreach (ReactionPtr r, model->getReactions()) {
		flux->setObj(r, 0);
	}
	fva(flux, min, max);
}

void fva(LPFluxPtr flux, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max ) {
	ModelPtr model = flux->getModel();

	flux->setObjSense(true);
	foreach(ReactionPtr a, model->getReactions()) {
		flux->setObj(a,1);
		flux->solvePrimal();
		max[a] = flux->getObjVal();
		flux->setObj(a,0);
	}
	foreach(ReactionPtr a, model->getReactions()) {
		flux->setObj(a,-1);
		flux->solvePrimal();
		min[a] = -flux->getObjVal();
		flux->setObj(a,0);
	}
}

} /* namespace metaopt */
