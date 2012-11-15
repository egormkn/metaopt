/*
 * BlockingSet.h
 *
 *  Created on: 15.11.2012
 *      Author: arnem
 */

#ifndef BLOCKINGSET_H_
#define BLOCKINGSET_H_

#include <vector>
#include "metaopt/model/scip/LPFlux.h"
#include "metaopt/model/DirectedReaction.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Given a thermodynamically flux vector, computes the set of reaction directions that need to be blocked
 * such that all LP solutions can be turned thermodynamically feasible by the cycle deletion heuristic.
 * This means, that the only remaining cycles must not contain objective nor flux forcing reactions.
 *
 * The algorithm will only block reaction directions that are not used by the given solution.
 *
 * This method only enforce the looplaw. If you have potentials this method will not work.
 *
 * @param flux the current solution of this LPFlux must be a thermodynamically feasible flux. Only the flux values are used. If you have a Solution from a ScipModel at hand, you can wrap this in an LPFlux.
 * @param model the LP that should be modified such that only thermodynamically feasible flux is possible. After the function terminates, this LP will contain exactly the computated solution.
 */
boost::shared_ptr<std::vector<DirectedReaction> > findBlockingSet(LPFluxPtr flux, LPFluxPtr model);

} /* namespace metaopt */
#endif /* BLOCKINGSET_H_ */
