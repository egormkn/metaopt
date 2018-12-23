/*
 * BlockingSet.h
 *
 *  Created on: 15.11.2012
 *      Author: arnem
 */

#ifndef BLOCKINGSET_H_
#define BLOCKINGSET_H_

#include <vector>
#include "model/scip/LPFlux.h"
#include "model/DirectedReaction.h"
#include "Properties.h"

namespace metaopt {

/**
 * Given a thermodynamically flux vector, computes the set of reaction directions that need to be blocked
 * such that all LP solutions can be turned thermodynamically feasible by the cycle deletion heuristic.
 * This means, that the only remaining cycles must not contain objective nor flux forcing nor problematic reactions.
 *
 * The algorithm will only block reaction directions that are not used by the given solution.
 *
 * This method only enforce the looplaw. If you have potentials this method will not work.
 *
 * @param flux the current solution of this LPFlux must be a thermodynamically feasible flux.
 * 		  Only the flux values are used. If you have a Solution from a ScipModel at hand, you can wrap this in an LPFlux.
 * 		  flux also supplies the underlying model (the model it was constructed for).
 * 		  It also supplies the objective sense.
 */
boost::shared_ptr<std::vector<DirectedReaction> > findBlockingSet(LPFluxPtr flux);

} /* namespace metaopt */
#endif /* BLOCKINGSET_H_ */
