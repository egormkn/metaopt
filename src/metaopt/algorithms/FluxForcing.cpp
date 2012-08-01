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
		assert(-EPSILON < rxn->getObj() && rxn->getObj() < EPSILON);
	}
	foreach(MetabolitePtr met, model->getMetabolites()) {
		assert(-EPSILON < met->getPotObj() && met->getPotObj() < EPSILON);
	}
#endif

	// compute \lambda^t(y_2)
	setObj(target,1);
	double old_source_ub = setUb(source, y_2);
	ScipModelPtr scip = factory.build(model);
	scip->solve();
	double max_t_y_2 = scip->getObjectiveValue();
	double foo = scip->getCurrentFlux(target._rxn);
	setObj(target,0);

	// compute \lambda^r_t(y_2)
	setObj(fluxForcing,1);
	double old_target_lb = setLb(target, max_t_y_2-EPSILON/100);
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
	if(lambda_t > EPSILON) {
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
	if(k > EPSILON) {
		// compute actual flux reduction
		setObj(target,1);
		double old_fluxFocing_lb = setLb(fluxForcing, k);
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

	double oldub = setUb(source, y_1);
	setObj(target,1);
	ScipModelPtr scip = factory.build(model);
	scip->solve();
	full = scip->getObjectiveValue();
	setUb(source, oldub);
	setObj(target, 0);

	// do the computation
	int i = 0;
	foreach(ReactionPtr rxn, model->getReactions()) {
		DirectedReaction d(rxn, true);
		double red = computeReducedFlux(model, factory, source, target, d, y_1, y_2);
		if(red + EPSILON <= full) {
			reduced[d] = red;
		}
		d = DirectedReaction(rxn, false);
		red = computeReducedFlux(model, factory, source, target, d, y_1, y_2);
		if(red + EPSILON <= full) {
			reduced[d] = red;
		}
		cout << endl << endl;
		i++;
		cout << "Computed " << i << " " << rxn->getName() << endl << endl;
	}
}


} /* namespace metaopt */
