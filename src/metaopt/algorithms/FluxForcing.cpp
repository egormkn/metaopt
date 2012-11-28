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
 * FluxForcing.cpp
 *
 *  Created on: 23.07.2012
 *      Author: arnem
 */

#include <iostream>
#include "FluxForcing.h"
#include "metaopt/model/DirectedReaction.h"
#include "metaopt/model/Model.h"
#include "metaopt/model/scip/ScipModel.h"

using namespace std;
using namespace boost;

namespace metaopt {

/**
 * sets upper bound for directed reaction.
 * @returns old bound value.
 */
double setUb(DirectedReaction d, double ub) {
	if(d._fwd) {
		double oldub = d._rxn->getUb();
		if(d._rxn->getLb() > ub) ub = d._rxn->getLb();
		d._rxn->setUb(ub);
		return oldub;
	}
	else {
		double oldub = -d._rxn->getLb();
		if(d._rxn->getUb() < -ub) ub = -d._rxn->getUb();
		d._rxn->setLb(-ub);
		return oldub;
	}
}

/**
 * sets lower bound for directed reaction.
 * @return old bound value
 */
double setLb(DirectedReaction d, double lb) {
	if(d._fwd) {
		double oldlb = d._rxn->getLb();
		if(d._rxn->getUb() < lb) lb = d._rxn->getUb();
		d._rxn->setLb(lb);
		return oldlb;
	}
	else {
		double oldlb = -d._rxn->getUb();
		if(d._rxn->getLb() > -lb) lb = -d._rxn->getLb();
		d._rxn->setUb(-lb);
		return oldlb;
	}
}

double setObj(DirectedReaction d, double obj) {
	if(d._fwd) {
		double oldobj = d._rxn->getObj();
		d._rxn->setObj(obj);
		return oldobj;
	}
	else {
		double oldobj = d._rxn->getObj();
		d._rxn->setObj(-obj);
		return oldobj;
	}
}

double computeReducedFlux(ModelPtr model, ModelFactory& factory, DirectedReaction source, DirectedReaction target, DirectedReaction fluxForcing, double y_1, double y_2) {
#ifndef NDEBUG
	foreach(ReactionPtr rxn, model->getReactions()) {
		assert(!rxn->isObjective());
	}
	foreach(MetabolitePtr met, model->getMetabolites()) {
		assert(!met->isObjective());
	}
#endif

	// we are dealing with fluxes, and will add constraints to fluxes
	// so we should start with more precision than actually desired.

	/* we run three consecutive optimization steps
	 *
	 * 1. Computation of max_t_y_2         - rangePrecision
	 * 2. Computation of lambda_t, lambda  - lambdaPrecision
	 * 3. Computation of maxFlux           - targetPrecision
	 */

	PrecisionPtr targetPrecision = model->getFluxPrecision();
	PrecisionPtr lambdaPrecision = targetPrecision->getSlavePrecision();
	model->setFluxPrecision(lambdaPrecision->getSlavePrecision()); // set it to rangePrecision

	// compute \lambda^t(y_2)
	setObj(target,1);
	double old_source_ub = setUb(source, y_2);
	ScipModelPtr scip = factory.build(model);
	scip->solve();
	double max_t_y_2 = scip->getObjectiveValue();
	double foo = scip->getCurrentFlux(target._rxn);
	setObj(target,0);

	// compute \lambda^r_t(y_2)
	// for this, we use a constraint that we computed in the previous optimization step,
	// so we now have to reduce the precision to lambda precision
	model->setFluxPrecision(lambdaPrecision);
	setObj(fluxForcing,1);
	//double old_target_lb = setLb(target, max_t_y_2-EPSILON/100);
	double old_target_lb = setLb(target, max_t_y_2);
	ScipModelPtr scip2 = factory.build(model);
	scip2->solve();
	double lambda_t;
	if(scip2->isOptimal()) {
		lambda_t = max((double) 0,scip2->getObjectiveValue());
	}
	else {
		/*SCIPprintBestSol(scip->getScip(), stdout, true);
		{
			setLb(target, max_t_y_2-EPSILON);
			ScipModelPtr scip3 = factory.build(model);
			scip3->solve();
			assert(false);
		}*/
		lambda_t = -1;
	}
	setLb(target, old_target_lb);

	double lambda = 0;
	if(lambda_t > lambdaPrecision->getCheckTol()) {
		// compute \lambda^r(y_1)
		setUb(source, y_1);
		scip = factory.build(model);
		scip->solve();
		if(scip->isOptimal()) {
			lambda = max(0.0, scip->getObjectiveValue());
		}
		else {
			lambda = -1;
		}
	}
	setObj(fluxForcing,0);

	double k = min(lambda_t, lambda);

	double maxFlux = INFINITY;
	if(k > lambdaPrecision->getCheckTol()) {
		// compute actual flux reduction
		setObj(target,1);
		double old_fluxFocing_lb = setLb(fluxForcing, k);
		model->setFluxPrecision(targetPrecision);
		scip = factory.build(model);
		scip->solve();
		if(scip->isOptimal()) {
			maxFlux = scip->getObjectiveValue();
		}
		else {
			maxFlux = -INFINITY;  // mark with impossibly good value -> indicates that something went wrong
		}
		setLb(fluxForcing, old_fluxFocing_lb);
		setObj(target,0);
	}
	else if(k < -0.5) {
		maxFlux = -INFINITY; // mark with impossibly good value -> indicates that something went wrong
	}
	// cleanup last bound
	setUb(source, old_source_ub);

	return maxFlux;
}

void computeReducedFluxes(ModelPtr model, ModelFactory& factory, DirectedReaction source, DirectedReaction target, double y_1, double y_2, boost::unordered_map<DirectedReaction, double>& reduced, double& full) {
	// clear objectives
	foreach(ReactionPtr rxn, model->getReactions()) {
		rxn->setObj(0);
	}
	foreach(MetabolitePtr met, model->getMetabolites()) {
		met->setPotObj(0);
	}

	/*
	 * Compute theoretical maximal flux.
	 * This is only used to find the reactions which should be outputted.
	 * But it is not used to constrain the feasible space, so we don't have to solve it with additional precision.
	 */
	double oldub = setUb(source, y_1);
	setObj(target,1);
	ScipModelPtr scip = factory.build(model);
	scip->solve();
	full = scip->getObjectiveValue();
	setUb(source, oldub);
	setObj(target, 0);

	/**
	 * Compute the reduced fluxes for the reactions.
	 */
	const PrecisionPtr& modelPrec = model->getFluxPrecision();
	// do the computation
	int i = 0;
	foreach(ReactionPtr rxn, model->getReactions()) {
		DirectedReaction d(rxn, true);
		double red = computeReducedFlux(model, factory, source, target, d, y_1, y_2);
		if(red + modelPrec->getCheckTol() <= full) {
			reduced[d] = red;
		}
		d = DirectedReaction(rxn, false);
		red = computeReducedFlux(model, factory, source, target, d, y_1, y_2);
		if(red + modelPrec->getCheckTol() <= full) {
			reduced[d] = red;
		}
		cout << endl << endl;
		i++;
		cout << "Computed " << i << " " << rxn->getName() << endl << endl;
	}
}


} /* namespace metaopt */
