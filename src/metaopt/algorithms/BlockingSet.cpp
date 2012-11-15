/*
 * BlockingSet.cpp
 *
 *  Created on: 15.11.2012
 *      Author: arnem
 */

#include "BlockingSet.h"

#include <vector>
#include <boost/unordered_map.hpp>
#include "metaopt/model/scip/LPFlux.h"
#include "metaopt/Properties.h"

using namespace std;
using namespace boost;

namespace metaopt {

struct BlockingSetScore {
	int score;
	vector<int> cycles;

	BlockingSetScore() : score(0) {}
};

shared_ptr<vector<DirectedReaction> > findBlockingSet(LPFluxPtr flux, LPFluxPtr model) {
	ModelPtr m = flux->getModel();
	assert(m == model->getModel());

	// by default test is set to maximize.
	// we will adjust test in such way that flow through objective and flux forcing reaction is maximized.
	LPFlux test(flux->getModel(), false);
	test.setZeroObj();
	double sense = model->isMaximize()? 1 : -1;

	foreach(ReactionPtr r, m->getInternalReactions()) {
		if(model->getLb(r) > EPSILON) {
			test.setObj(r, 1);
		}
		else if(model->getUb(r) < -EPSILON) {
			test.setObj(r, -1);
		}
		else if(model->getObj(r) > EPSILON) {
			test.setObj(r, sense);
		}
		else if(model->getObj(r) < -EPSILON) {
			test.setObj(r, -sense);
		}
	}

	unordered_map<DirectedReaction, BlockingSetScore> scores;
	vector<unordered_set<DirectedReaction> > cycles;
	shared_ptr<vector<DirectedReaction> > blocked(new vector<DirectedReaction>);

	// block reactions until no cycle is found
	bool foundCycle = true;
	// use resolve to mark, if we did not solve the LPs to optimality and we hence have to resolve
	// it also counts how often we ran resolve after each other
	int resolve = 0;

	while(foundCycle) {
		if(resolve > 0) {
			model->solveDual();
			test.solveDual();
		}
		else {
			model->solveDual();
			if(model->isFeasible()) {
				test.setDirectionBounds(model); // only allow cycle-flows on reactions that carry flux
				test.solveDual();
			}
		}
		// we do not need optimality to run this
		if(model->isFeasible() && test.isFeasible() && test.getObjVal() > EPSILON) {
			int cycleId = cycles.size(); // index of the cycle that we will create
			cycles.push_back(unordered_set<DirectedReaction>());
			unordered_set<DirectedReaction>& cycle = cycles.back();
			ReactionPtr best;
			int best_score = 0;
			// collect reactions in the cycle
			// and update scores
			// we also find the reaction with the best score, because that one, we will block
			foreach(ReactionPtr r, m->getInternalReactions()) {
				if(test.getFlux(r) > EPSILON && flux->getFlux(r) < EPSILON) {
					DirectedReaction d(r, true);
					cycle.insert(d);
					BlockingSetScore& s = scores[d];
					s.score++;
					s.cycles.push_back(cycleId);
					if(s.score > best_score) {
						best_score = s.score;
						best = r;
					}
				}
				if(test.getFlux(r) < -EPSILON && flux->getFlux(r) > -EPSILON) {
					DirectedReaction d(r, false);
					cycle.insert(d);
					BlockingSetScore& s = scores[d];
					s.score++;
					s.cycles.push_back(cycleId);
					if(s.score > best_score) {
						best_score = s.score;
						best = r;
					}
				}
			}
			assert(best_score > 0); // we must find at least one reaction

			// block flow through best in its current direction
			if(test.getFlux(best) > EPSILON) {
				assert(model->getFlux(best) > EPSILON);
				model->setUb(best, 0);
				blocked->push_back(DirectedReaction(best, true));
			}
			else {
				assert(test.getFlux(best) < -EPSILON);
				assert(model->getFlux(best) < -EPSILON);
				model->setLb(best, 0);
				blocked->push_back(DirectedReaction(best, false));
			}
			resolve = 0; // we changed something on the bounds, so we cannot only resolve the old LPs
		}
		else {
			// to decide that we are done, we need optimality
			if(model->isOptimal() && test.isOptimal()) {
				foundCycle = false;
			}
			else {
				// LPs are not solved to optimality, so just resolve
				resolve++;
				// this may end in an endless loop, hence we abort if resolve is increased too much
				if(resolve > 10) {
					assert(false);
					BOOST_THROW_EXCEPTION(EndlessLoopError() << iteration_count(10));
				}
			}
		}
	}

	/**
	 * Possible Improvement:
	 * Some reactions are contained in several cycles.
	 * We could now try to improve the algorithm by first blocking the reactions with highest score.
	 * This may allow us to use blocking reactions with more effect. And thus, we may be able to use less and thus restrict the solution space less.
	 *
	 * But, it is unclear how useful it really is, so I'm not doing it now.
	 * This also has the effect, that the cycles map is effectively unused.
	 */

	return blocked;
}

} /* namespace metaopt */
