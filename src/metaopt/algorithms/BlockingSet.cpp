/*
 * BlockingSet.cpp
 *
 *  Created on: 15.11.2012
 *      Author: arnem
 */

#include "BlockingSet.h"

#include <vector>
#include <iostream>
#include <boost/unordered_map.hpp>
#include "metaopt/model/scip/LPFlux.h"
#include "metaopt/Properties.h"

using namespace std;
using namespace boost;

namespace metaopt {


#if 0
// first implementation, I think we can do better
// it also contains a bug

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
#else

// new, hopefully more clever implementation

shared_ptr<vector<DirectedReaction> > findBlockingSet(LPFluxPtr test, LPFluxPtr flux) {
	// preconditions:
	// test is already solved to optimality
	assert(test->isOptimal());

	shared_ptr<vector<DirectedReaction> > block(new vector<DirectedReaction>);

	foreach(ReactionPtr rxn, flux->getModel()->getInternalReactions()) {
		int stat = test->getColumnStatus(rxn);
		if(stat == SCIP_BASESTAT_LOWER && test->getLb(rxn) < -EPSILON && flux->getFlux(rxn) > -EPSILON) {
			cout << rxn->getName() << " redcost=" << test->getReducedCost(rxn) << endl;
			// we can block negative flux through this reaction
			double redcost = test->getReducedCost(rxn);
			if(redcost > EPSILON || redcost < -EPSILON) { // is this check valid?
				block->push_back(DirectedReaction(rxn, false));
			}
		}
		else if(stat == SCIP_BASESTAT_UPPER && test->getUb(rxn) > EPSILON && flux->getFlux(rxn) < EPSILON) {
			cout << rxn->getName() << " redcost=" << test->getReducedCost(rxn) << endl;
			// we can block positive flux through this reaction
			double redcost = test->getReducedCost(rxn);
			if(redcost > EPSILON || redcost < -EPSILON) { // is this check valid?
				block->push_back(DirectedReaction(rxn, true));
			}
		}
	}

	// since we now blocked the maximal flow through the network, there is no flow possible anymore and we have found the desired blocking set.

	return block;
}

shared_ptr<vector<DirectedReaction> > findBlockingSet(LPFluxPtr flux) {
	ModelPtr model = flux->getModel();
	double sense = flux->isMaximize() ? 1 : -1;

	LPFluxPtr test(new LPFlux(flux->getModel(), false));

	// first set the bounds of the test flux apropriately
	foreach(ReactionPtr rxn, model->getInternalReactions()) {
		if(rxn->getUb() > EPSILON) {
			// if positive flux is possible but not used, create a bound of 1 to restrict positive flow
			if(flux->getFlux(rxn) < EPSILON) test->setUb(rxn, 1);
			// otherwise we cannot block this reaction and hence, do not create a bound
			else test->setUb(rxn, INFINITY);
		}
		else {
			test->setUb(rxn, 0);
		}
		if(rxn->getLb() < -EPSILON) {
			// if negative flux is possible but not used, create a bound of 1 to restrict positive flow
			if(flux->getFlux(rxn) > -EPSILON) test->setLb(rxn, -1);
			// otherwise we cannot block this reaction and hence, do not create a bound
			else test->setLb(rxn, -INFINITY);
		}
		else {
			test->setLb(rxn, 0);
		}
	}
	// reset the objective
	test->setZeroObj();

	shared_ptr<vector<DirectedReaction> > block(new vector<DirectedReaction>);

	/*
	 * Strategy:
	 * We start blocking flux through unimportant reactions first.
	 * This makes sure, that for the important reactions we only block as few reactions as needed.
	 * So, we start with the problematic reactions, then proceed to the flux forcing reactions and finally deal with objective reactions.
	 */

	foreach(ReactionPtr rxn, model->getProblematicReactions()) {
		// with flux forcing reactions we deal later
		// for objective reactions, we only check the reaction unfavoured by the objective function
		if(!rxn->isFluxForcing()) {
			// positive flux case
			if(rxn->getUb() > EPSILON && rxn->getObj() < EPSILON) { // if the objective is also positive, we deal with it later
				test->setObj(rxn, 1); // ask for flow through the problematic reaction
				test->solvePrimal();

				double f = flux->getFlux(rxn);
				if(test->getObjVal() > EPSILON) {
					if(f < EPSILON) {
						// the reaction really needs to be blocked,
						// but it can be blocked itself.
						// so just block the beast itself
						block->push_back(DirectedReaction(rxn, true));
						test->setUb(rxn, 0);
					}
					else {
						// we are not allowed to block the beast itself
						// so compute a blocking set
						shared_ptr<vector<DirectedReaction> > b = findBlockingSet(test, flux);
						foreach(DirectedReaction& d, *b) {
							block->push_back(d);
							if(d._fwd) test->setUb(d._rxn, 0);
							else test->setLb(d._rxn, 0);
						}
					}
				}
			}
			// negative flux case
			if(rxn->getLb() < -EPSILON && rxn->getObj() > -EPSILON) { // if the objective is also negative, we deal with it later
				test->setObj(rxn, -1); // ask for flow through the problematic reaction
				test->solvePrimal();

				double f = flux->getFlux(rxn);
				if(test->getObjVal() > EPSILON) {
					// the reaction really needs to be blocked,
					if(f > -EPSILON) {
						// but it can be blocked itself.
						// so just block the beast itself
						block->push_back(DirectedReaction(rxn, false));
						test->setLb(rxn, 0);
					}
					else {
						// we are not allowed to block the beast itself
						// so compute a blocking set
						shared_ptr<vector<DirectedReaction> > b = findBlockingSet(test, flux);
						foreach(DirectedReaction& d, *b) {
							block->push_back(d);
							if(d._fwd) test->setUb(d._rxn, 0);
							else test->setLb(d._rxn, 0);
						}
					}
				}
			}
			// we blocked the problematic reaction, if it needed to be blocked

			// undo changes to objective function
			test->setObj(rxn, 0);
		}
	}
	// we dealt with the problematic reactions
	// now we process the flux forcing reactions
	foreach(ReactionPtr rxn, model->getFluxForcingReactions()) {
		if(rxn->getLb() > EPSILON) {
			test->setObj(rxn, 1); // ask for flow through the flux forcing reaction
		}
		else {
			assert(rxn->getUb() < -EPSILON);
			test->setObj(rxn, -1);
		}
		test->solvePrimal();
		// we always have to compute the blocking set
		shared_ptr<vector<DirectedReaction> > b = findBlockingSet(test, flux);
		// and enforce it
		foreach(DirectedReaction& d, *b) {
			block->push_back(d);
			if(d._fwd) test->setUb(d._rxn, 0);
			else test->setLb(d._rxn, 0);
		}
		// undo the objective change
		test->setObj(rxn, 0);
	}
	// we have now dealt with flux-forcing and problematic reactions
	// now deal with the objective reactions
	// theoretically we can deal with these reactions all together
	// however, some of the objective reactions may also be problematic and then we have to deal with them separately.
	// so, we deal with one objective reaction after another.
	// this also makes sure, that we don't get problems with cycles consisting only of objective reactions
	foreach(ReactionPtr rxn, model->getObjectiveReactions()) {
		if(rxn->getObj() > EPSILON) {
			test->setObj(rxn, sense); // ask for flow through the objective reaction according to the sense of the optimization problem
		}
		else {
			test->setObj(rxn, -sense);
		}
		test->solvePrimal();
		// we always have to compute the blocking set
		shared_ptr<vector<DirectedReaction> > b = findBlockingSet(test, flux);
		// and enforce it
		foreach(DirectedReaction& d, *b) {
			block->push_back(d);
			if(d._fwd) test->setUb(d._rxn, 0);
			else test->setLb(d._rxn, 0);
		}
		// undo the objective change
		test->setObj(rxn, 0);
	}

	// we are done and we can return the blocking set!

	return block;
}

#endif



} /* namespace metaopt */
