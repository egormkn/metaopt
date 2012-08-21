/*
 * FCA.cpp
 *
 *  Created on: 14.07.2012
 *      Author: arne
 */

#include "FCA.h"

#include <iostream>

#include "metaopt/algorithms/FVA.h"

// the condition of metabolic networks can be very bad. In particular here, we need some extra factor.
#define SCALE_FACTOR 10000

using namespace boost;
using namespace std;

namespace metaopt {

void fca(ModelPtr model, CouplingPtr coupling) {
	/**
	 * this is the naive, stupid implementation.
	 * TODO: a smart implementation (e.g. like in F2C2 by Larhlimi et al. 2012)
	 */

	unordered_map<ReactionPtr,double > min, max;

	LPFluxPtr test(new LPFlux(model, true));
	foreach (ReactionPtr a, model->getReactions()) {
		test->setObj(a, 0);
		if(a->getLb() < -EPSILON) test->setLb(a,-SCALE_FACTOR);
		else test->setLb(a,0);
		if(a->getUb() > EPSILON) test->setUb(a,SCALE_FACTOR);
		else test->setUb(a,0);

	}
	// compute blocked fluxes to save computation time
	fva(test, min, max);

	int count = 0;
	foreach(ReactionPtr a, model->getReactions()) {
		if(min[a] < -EPSILON) {
			test->setLb(a,0);
			foreach(ReactionPtr b, model->getReactions()) {
				if(max[b] > EPSILON) {
					test->setObj(b,1);
					test->solvePrimal();
					if(test->getObjVal() < EPSILON) {
						// b^+ -> a^-
						coupling->addCoupled(DirectedReaction(b,true), DirectedReaction(a,false));
					}
					test->setObj(b,0);
				}
			}
			foreach(ReactionPtr b, model->getReactions()) {
				if(min[b] < -EPSILON) {
					test->setObj(b,-1);
					test->solvePrimal();
					if(test->getObjVal() < EPSILON) {
						// b^- -> a^-
						coupling->addCoupled(DirectedReaction(b,false), DirectedReaction(a,false));
					}
					test->setObj(b,0);
				}
			}
			test->setLb(a,-SCALE_FACTOR);
		}
		if(max[a] > EPSILON) {
			test->setUb(a,0);
			foreach(ReactionPtr b, model->getReactions()) {
				if(max[b] > EPSILON) {
					test->setObj(b,1);
					test->solvePrimal();
					if(test->getObjVal() < EPSILON) {
						// b^+ -> a^+
						if(test->getObjVal() > EPSILON*EPSILON) {
							cout << b->getName() + "_fwd -> " + a->getName() + "_fwd" << " " << test->getObjVal() << endl;
						}
						coupling->addCoupled(DirectedReaction(b,true), DirectedReaction(a,true));
					}
					test->setObj(b,0);
				}
			}
			foreach(ReactionPtr b, model->getReactions()) {
				if(min[b] < -EPSILON) {
					test->setObj(b,-1);
					test->solvePrimal();
					if(test->getObjVal() < EPSILON) {
						// b^- -> a^+
						coupling->addCoupled(DirectedReaction(b,false), DirectedReaction(a,true));
					}
					test->setObj(b,0);
				}
			}
			test->setUb(a,SCALE_FACTOR);
		}
#ifndef SILENT
		count++;
		cout << count << " ";
		cout.flush();
		if(count % 20 == 0) cout << endl;
#endif
	}
#ifndef SILENT
	cout << coupling->getStat() << endl;
#endif
}

} /* namespace metaopt */
