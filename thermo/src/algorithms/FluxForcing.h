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
 * FluxForcing.h
 *
 *  Created on: 23.07.2012
 *      Author: arnem
 */

#ifndef FLUXFORCING_H_
#define FLUXFORCING_H_

#include <boost/unordered_map.hpp>
#include "model/DirectedReaction.h"
#include "model/Model.h"
#include "model/scip/ScipModel.h"
#include "algorithms/ModelFactory.h"
#include "Properties.h"


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
