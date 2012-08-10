/*
 * FluxForcing.h
 *
 *  Created on: 23.07.2012
 *      Author: arnem
 */

#ifndef FLUXFORCING_H_
#define FLUXFORCING_H_

#include <boost/unordered_map.hpp>
#include "metaopt/model/DirectedReaction.h"
#include "metaopt/model/Model.h"
#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/algorithms/ModelFactory.h"
#include "metaopt/Properties.h"


namespace metaopt {

/**
 * By making the reaction <code>fluxForcing</code> flux forcing how much flux can we reduce through the target reaction?
 *
 * This function assumes that no variable has objective coefficients initially.
 * After the execution of the function, the model will be back in its old state.
 *
 * @param model the model to operate on
 * @param factory gives the ScipModel instance
 * @param source the source reaction (e.g. light)
 * @param target the target reaction (e.g. biomass)
 * @param fluxForcing the reaction to make fluxForcing
 * @param y_1 the level for which to reduce maximal target production rate
 * @param y_2 the level for which the target production rate may not be minimized
 * @return maximal flux through target as possible when using fluxForcing as much as possible. If no flux reduction is detected to be possible, INFINITY is returned.
 */
double computeReducedFlux(ModelPtr model, ModelFactory& factory, DirectedReaction source, DirectedReaction target, DirectedReaction fluxForcing, double y_1, double y_2);

/**
 * computes reduced fluxes for all reaction in the model.
 * Only those reactions that allow a reduction are stored in the output.
 *
 * This function clears objective coefficients of reactions and metabolite potentials.
 */
void computeReducedFluxes(ModelPtr model, ModelFactory& factory, DirectedReaction source, DirectedReaction target, double y_1, double y_2, boost::unordered_map<DirectedReaction, double>& reduced, double& full);


} /* namespace metaopt */
#endif /* FLUXFORCING_H_ */
