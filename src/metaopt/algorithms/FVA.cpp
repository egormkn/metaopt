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
 * FVA.cpp
 *
 *  Created on: 14.07.2012
 *      Author: arne
 */

#include <iostream>
#include <time.h>
#include "scip/scip.h"

#include "FVA.h"
#include "metaopt/Properties.h"
#include "metaopt/model/scip/LPPotentials.h"
#include "metaopt/scip/constraints/SteadyStateConstraint.h"
#include "metaopt/scip/constraints/ThermoConstraintHandler.h"
#include "metaopt/scip/heur/CycleDeletionHeur.h"

using namespace boost;

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
		if(!flux->isOptimal()) {
			flux->solvePrimal();
		}
		assert(flux->isOptimal());
		max[a] = flux->getObjVal();
		flux->setObj(a,0);
	}
	foreach(ReactionPtr a, model->getReactions()) {
		flux->setObj(a,-1);
		flux->solvePrimal();
		if(!flux->isOptimal()) {
			flux->solvePrimal();
		}
		assert(flux->isOptimal());

		min[a] = -flux->getObjVal();
		flux->setObj(a,0);
	}
}

void fva(ModelPtr model, ModelFactory& factory, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max ) {
	foreach(ReactionPtr a, model->getReactions()) {
		a->setObj(0);
	}
	foreach(MetabolitePtr a, model->getMetabolites()) {
		a->setPotObj(0);
	}

	int i = 0;
	foreach(ReactionPtr a, model->getReactions()) {
		a->setObj(1);
		ScipModelPtr scip = factory.build(model);
		scip->setObjectiveSense(true);
		scip->solve();
		max[a] = scip->getObjectiveValue();
		scip = factory.build(model);
		scip->setObjectiveSense(false);
		scip->solve();
		min[a] = scip->getObjectiveValue();
		a->setObj(0);

		std::cout << std::endl;
		std::cout << "solved " << ++i << " : " << a->getName() << std::endl;
		std::cout << std::endl;
	}
}

class FVAThermoModelFactory : public ModelFactory {
public:
	ScipModelPtr build(ModelPtr m) {
		ScipModelPtr scip(new ScipModel(m));
		createSteadyStateConstraint(scip);
		createThermoConstraint(scip);
		createCycleDeletionHeur(scip);
		//registerExitEventHandler(scip);

		return scip;
	}
};

/**
 * checks if a loopless free flux with the same objective value can be attained
 */
bool isLooplessFluxAttainable(LPFluxPtr sol, LPFluxPtr helper) {
	ModelPtr model = sol->getModel();
	helper->setDirectionBounds(sol);
	helper->setZeroObj();
	foreach(ReactionPtr r, model->getObjectiveReactions()) {
		if(!r->isExchange()) {
			double val = sol->getFlux(r);
			if(val > EPSILON) {
				helper->setObj(r, 1);
			}
			else if(val < -EPSILON) {
				helper->setObj(r, -1);
			}
		}
	}
	foreach(ReactionPtr r, model->getFluxForcingReactions()) {
		if(!r->isExchange()) {
			double val = sol->getFlux(r);
			if(val > EPSILON) {
				helper->setObj(r, 1);
			}
			else if(val < -EPSILON) {
				helper->setObj(r, -1);
			}
		}
	}
	foreach(ReactionPtr r, model->getProblematicReactions()) {
		if(!r->isExchange()) {
			double val = sol->getFlux(r);
			if(val > EPSILON) {
				helper->setObj(r, 1);
			}
			else if(val < -EPSILON) {
				helper->setObj(r, -1);
			}
		}
	}
	helper->solve();
	if(!helper->isOptimal()) { // if we weren't able to compute the optimum, we better be pessimistic
		return false;
	}
#ifndef SILENT
	std::cout << "test: " << helper->getObjVal() << std::endl;
#endif
	return helper->getObjVal() < EPSILON;
}

/**
 * checks if a thermodynamically feasible flux with the same objective value can be attained
 */
bool isThermoFluxAttainable(LPFluxPtr sol, LPFluxPtr helper, LPPotentialsPtr potTest) {
	// make sure that subtracting cycles does not violate flux bounds or objective value
	assert(isLooplessFluxAttainable(sol, helper));
	ModelPtr model = sol->getModel();
#ifndef NDEBUG
	int debugi = 0;
#endif
	do {
		helper->setDirectionBounds(sol);
		helper->setZeroObj();
		foreach(ReactionPtr r, model->getReactions()) {
			double val = sol->getFlux(r);
			if(val > EPSILON) {
				helper->setObj(r, 1);
			}
			else if(val < -EPSILON) {
				helper->setObj(r, -1);
			}
		}
		helper->solve();
		if(!helper->isFeasible()) { // something strange, abort
			return false;
		}
		if(helper->getObjVal() > EPSILON) {
			double scale = sol->getSubScale(helper);
			if(scale < -0.5) return false;
			sol->subtract(helper, scale);
		}
#ifndef NDEBUG
		debugi++;
		if(debugi % 1000 == 0) std::cout << "FVA.cpp " << __LINE__ << " caught endless loop" << std::endl;
#endif
	} while(helper->getObjVal() > EPSILON);
	// now we computed a guess of a thermodynamically feasible flow, now we have to check
	potTest->setDirections(sol);
	bool result;
	if(!potTest->testStrictFeasible(result)) {
		return false;
	}
	else return result;
}

void tfva(ModelPtr model, FVASettingsPtr settings, unordered_map<ReactionPtr,double >& min , unordered_map<ReactionPtr,double >& max ) {
	/*
	 * Depending on the kind of thermodynamic information given there are different kinds of speedups possible.
	 * If no potential bounds are given, we just have loopless FVA and we can check if FBA = tFBA by solving two LPs
	 * In the other case, we basically have to run the cycle elimination heur
	 */

	clock_t start = clock();
	double runningTime = 0;

	bool simple = true;

	foreach(MetabolitePtr met, model->getMetabolites()) {
		if(isinf(met->getPotLb()) != -1) {
			simple = false;
		}
		if(isinf(met->getPotUb()) != 1) {
			simple = false;
		}
	}

	if(simple) std::cout << "tfva problem has simple structure " << std::endl;

	FVAThermoModelFactory factory;

	/**
	 * reset objective functions
	 */
	foreach(ReactionPtr a, model->getReactions()) {
		a->setObj(0);
	}
	foreach(MetabolitePtr a, model->getMetabolites()) {
		a->setPotObj(0);
	}

	/**
	 * for the cases where it is sufficient to run an LP, we just run an LP
	 */
	LPFluxPtr flux(new LPFlux(model, true));

	// the helper flux is used to test if the LP solution is already optimal
	LPFluxPtr helper(new LPFlux(model, false));
	helper->setObjSense(true);

	LPPotentialsPtr potTest;
	if(!simple) {
		potTest = LPPotentialsPtr(new LPPotentials(model));
	}

	flux->setObjSense(true);
	foreach(ReactionPtr a, settings->reactions) {
		a->setObj(1);
		cout << "max " << a->getName() << endl;
		flux->setObj(a,1);
		flux->solvePrimal();
#ifndef NDEBUG
		if(flux->isOptimal()) {
			cout << "opt-flux = " << flux->getObjVal() << endl;
		}
#endif
		// for shortcut looplessflux must always be attainable
		// if it is simple, it is sufficient, else we have to do more
		if(flux->isOptimal() && isLooplessFluxAttainable(flux, helper) && (simple || isThermoFluxAttainable(flux, helper, potTest))) {
			max[a] = flux->getObjVal();
		}
		else {
			ScipModelPtr scip = factory.build(model);
			if(settings->timeout > EPSILON) {
				BOOST_SCIP_CALL( SCIPsetRealParam(scip->getScip(), "limits/time", settings->timeout) );
			}
			scip->setObjectiveSense(true);
			scip->solve();
			assert(scip->isOptimal());
			max[a] = scip->getObjectiveValue();
		}
		flux->setObj(a,0);
		a->setObj(0);

		runningTime = (double) (clock() - start) / CLOCKS_PER_SEC;

		if(settings->timeout > EPSILON && runningTime > settings->timeout) {
			cout << endl;
			cout << "aborted by timeout of " << settings->timeout << " seconds" << endl;
			BOOST_THROW_EXCEPTION( TimeoutError() );
		}
	}

	flux->setObjSense(false);
	foreach(ReactionPtr a, settings->reactions) {
		a->setObj(1);
		cout << "min " << a->getName() << endl;
		flux->setObj(a,1);
		flux->solvePrimal();
		// for shortcut looplessflux must always be attainable
		// if it is simple, it is sufficient, else we have to do more
		if(flux->isOptimal() && isLooplessFluxAttainable(flux, helper) && (simple || isThermoFluxAttainable(flux, helper, potTest))) {
			min[a] = flux->getObjVal();
		}
		else {
			ScipModelPtr scip = factory.build(model);
			if(settings->timeout > EPSILON) {
				BOOST_SCIP_CALL( SCIPsetRealParam(scip->getScip(), "limits/time", settings->timeout) );
			}
			scip->setObjectiveSense(false);
			scip->solve();
			assert(scip->isOptimal());
			min[a] = scip->getObjectiveValue();
		}
		flux->setObj(a,0);
		a->setObj(0);

		runningTime = (double) (clock() - start) / CLOCKS_PER_SEC;

		if(settings->timeout > EPSILON && runningTime > settings->timeout) {
			cout << endl;
			cout << "aborted by timeout of " << settings->timeout << " seconds" << endl;
			BOOST_THROW_EXCEPTION( TimeoutError() );
		}

	}
}


} /* namespace metaopt */
