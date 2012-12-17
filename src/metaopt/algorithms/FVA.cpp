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
#include <fstream>
#include <time.h>
#include <math.h>
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
	CouplingPtr coupling;

	ScipModelPtr build(ModelPtr m) {
		ScipModelPtr scip(new ScipModel(m));
		createSteadyStateConstraint(scip);
		if(coupling.use_count() >= 1) {
			createThermoConstraint(scip, coupling);
		}
		else {
			createThermoConstraint(scip);
		}
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
	const PrecisionPtr& solPrec = sol->getPrecision();
	helper->setDirectionBounds(sol);
	helper->setZeroObj();
	foreach(ReactionPtr r, model->getObjectiveReactions()) {
		if(!r->isExchange()) {
			double val = sol->getFlux(r);
			if(val > solPrec->getCheckTol()) {
				helper->setObj(r, 1);
			}
			else if(val < -solPrec->getCheckTol()) {
				helper->setObj(r, -1);
			}
		}
	}
	foreach(ReactionPtr r, model->getFluxForcingReactions()) {
		if(!r->isExchange()) {
			double val = sol->getFlux(r);
			if(val > solPrec->getCheckTol()) {
				helper->setObj(r, 1);
			}
			else if(val < -solPrec->getCheckTol()) {
				helper->setObj(r, -1);
			}
		}
	}
	foreach(ReactionPtr r, model->getProblematicReactions()) {
		if(!r->isExchange()) {
			double val = sol->getFlux(r);
			if(val > solPrec->getCheckTol()) {
				helper->setObj(r, 1);
			}
			else if(val < -solPrec->getCheckTol()) {
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
	return helper->getObjVal() < helper->getPrecision()->getCheckTol();
}

/**
 * checks if a thermodynamically feasible flux with the same objective value can be attained
 */
bool isThermoFluxAttainable(LPFluxPtr sol, LPFluxPtr helper, LPPotentialsPtr potTest) {
	// make sure that subtracting cycles does not violate flux bounds or objective value
	assert(isLooplessFluxAttainable(sol, helper));
	ModelPtr model = sol->getModel();
	const PrecisionPtr& solPrec = sol->getPrecision();
	const PrecisionPtr& helperPrec = helper->getPrecision();
#ifndef NDEBUG
	int debugi = 0;
#endif
	do {
		helper->setDirectionBounds(sol);
		helper->setZeroObj();
		foreach(ReactionPtr r, model->getInternalReactions()) {
			double val = sol->getFlux(r);
			if(val > solPrec->getCheckTol()) {
				helper->setObj(r, 1);
			}
			else if(val < -solPrec->getCheckTol()) {
				helper->setObj(r, -1);
			}
		}
		helper->solve();
		if(!helper->isFeasible()) { // something strange, abort
			return false;
		}
		if(helper->getObjVal() > helperPrec->getCheckTol()) {
			double scale = sol->getSubScale(helper);
			if(scale < -0.5) return false;
			sol->subtract(helper, scale);
		}
#ifndef NDEBUG
		debugi++;
		if(debugi % 1000 == 0) std::cout << "FVA.cpp " << __LINE__ << " caught endless loop" << std::endl;
#endif
	} while(helper->getObjVal() > helperPrec->getCheckTol());
	// now we computed a guess of a thermodynamically feasible flow, now we have to check
	potTest->setDirections(sol);
	bool result;
	if(!potTest->testStrictFeasible(result)) {
		return false;
	}
	else return result;
}

int foo = 0;

void tfva(ModelPtr model, FVASettingsPtr settings, unordered_map<ReactionPtr,double >& min , unordered_map<ReactionPtr,double >& max ) {
	/*
	 * Depending on the kind of thermodynamic information given there are different kinds of speedups possible.
	 * If no potential bounds are given, we just have loopless FVA and we can check if FBA = tFBA by solving two LPs
	 * In the other case, we basically have to run the cycle elimination heur
	 */

#define REDUCE_DOMAIN 0

	clock_t start = clock();
	double runningTime = 0;

	bool simple = true;

	foreach(MetabolitePtr met, model->getMetabolites()) {
		if(isinf(met->getPotLb()) == 0) {
			simple = false;
		}
		if(isinf(met->getPotUb()) == 0) {
			simple = false;
		}
	}

	if(simple) std::cout << "tfva problem has simple structure " << std::endl;

	// TODO: Sometimes it is important that the results computed by FVA are not only valid bounds but also feasible.
	// in those cases we should check the computed solution for feasibility and if necessary make it feasible.

	// increase dual model flux precision, since we use computation results to add constraints
	PrecisionPtr orig_precision = model->getFluxPrecision();
	model->setFluxPrecision(orig_precision->getDualSlavePrecision());

	FVAThermoModelFactory factory;
	factory.coupling = settings->coupling;

	/**
	 * reset objective functions
	 */
	foreach(ReactionPtr a, model->getReactions()) {
		a->setObj(0);
	}
	foreach(MetabolitePtr a, model->getMetabolites()) {
		a->setPotObj(0);
	}

	/*
	 * for the cases where it is sufficient to run an LP, we just run an LP
	 *
	 * We will maximize fluxes and minimizes fluxes in parallel.
	 * This way we can exploit that we do not need to solve CPs for reactions which are already fixed by LP
	 * We use two different LPs for it, so that we don't have to change the objective too much.
	 */
	LPFluxPtr max_flux(new LPFlux(model, true));
	LPFluxPtr min_flux(new LPFlux(model, true));

	// the helper flux is used to test if the LP solution is already optimal
	LPFluxPtr helper(new LPFlux(model, false));
	helper->setObjSense(true);

	// helper needs more precision, since it is called iteratively to remove loops
	helper->setPrecision(model->getFluxPrecision()->getPrimalSlavePrecision());
	// for all the other cases we just use model precision, i.e.
	PrecisionPtr prec = model->getFluxPrecision();

	LPPotentialsPtr potTest;
	if(!simple) {
		potTest = LPPotentialsPtr(new LPPotentials(model));
	}

	max_flux->setObjSense(true);
	min_flux->setObjSense(false);

#if 0
	stringstream ss;
	ss << "debug" << ++foo <<".lp";

	SCIPlpiWriteLP(min_flux->getLPI(), ss.str().c_str());

	ss.clear();
	ss << "debug" << foo << ".map";

	ofstream mapf(ss.str().c_str());
	mapf << "Reactions" << endl;
	foreach(ReactionPtr a, model->getReactions()) {
		mapf << a->getName() << " " << min_flux->getIndex(a) << endl;
	}
	mapf << "Metabolites" << endl;
	foreach(MetabolitePtr a, model->getMetabolites()) {
		if(!a->hasBoundaryCondition()) {
			mapf << a->getName() << " " << min_flux->getIndex(a) << endl;
		}
	}
	mapf.close();
#endif

	int i = 1;
	int num_rxns = settings->reactions.size();

	foreach(ReactionPtr a, settings->reactions) {
		/*
		 * First solve the LPs
		 */
		a->setObj(1);
		cout << "max " << a->getName() << endl;
		max_flux->setObj(a,1);
		max_flux->solvePrimal();
#ifndef NDEBUG
		if(max_flux->isOptimal()) {
			cout << "opt-flux = " << max_flux->getObjVal() << endl;
		}
#endif

		cout << "min " << a->getName() << endl;
		min_flux->setObj(a,1);
		min_flux->solvePrimal();
#ifndef NDEBUG
		if(min_flux->isOptimal()) {
			cout << "opt-flux = " << min_flux->getObjVal() << endl;
		}
#endif

		if(max_flux->isFeasible() && min_flux->isFeasible() && max_flux->getObjVal() - min_flux->getObjVal() < prec->getCheckTol()) {
			// both fluxes are feasible and the same, hence they also must be optimal
			max[a] = max_flux->getObjVal();
			min[a] = min_flux->getObjVal();
		}
		else {
			// we were not able to derive that we already found the optimum, so deal with min and max separately

			/*
			 * Maximization
			 */

			// for shortcut looplessflux must always be attainable
			// if it is simple, it is sufficient, else we have to do more
			if(max_flux->isOptimal() && isLooplessFluxAttainable(max_flux, helper) && (simple || isThermoFluxAttainable(max_flux, helper, potTest))) {
				max[a] = max_flux->getObjVal();
			}
			else {
				ScipModelPtr scip = factory.build(model);
				if(settings->timeout > 1) { // a timeout of less than a second makes no sense
					BOOST_SCIP_CALL( SCIPsetRealParam(scip->getScip(), "limits/time", settings->timeout) );
				}
				scip->setObjectiveSense(true);
				scip->solve();

#ifndef NDEBUG
				if(!scip->isOptimal()) {
					SCIPprintStatistics(scip->getScip(), NULL);
					cout << "LP primal infeasible = " << SCIPlpiIsPrimalInfeasible(max_flux->getLPI()) << endl;
					if( SCIPlpiHasDualRay(max_flux->getLPI())) {
						double dualfarkas[model->getMetabolites().size()];
						SCIPlpiGetDualfarkas(max_flux->getLPI(), dualfarkas);
						foreach(MetabolitePtr met, model->getMetabolites()) {
							double d = dualfarkas[max_flux->getIndex(met)];
							if(d < -prec->getDualFeasTol() || d > prec->getDualFeasTol()) {
								cout << met->getName() << " = " << d << endl;
							}
						}
					}
					else {
						cout << "no dual ray available" << endl;
					}
					SCIPwriteOrigProblem(scip->getScip(), "debug.lp", NULL, TRUE);
				}
#endif
				assert(scip->isOptimal());
				double opt = scip->getObjectiveValue();
				max[a] = opt;
#if REDUCE_DOMAIN
				a->setUb(opt); // we computed an upper bound, so use it for future computations
#endif
			}

			/*
			 * Minimization
			 */

			// for shortcut looplessflux must always be attainable
			// if it is simple, it is sufficient, else we have to do more
			if(min_flux->isOptimal() && isLooplessFluxAttainable(min_flux, helper) && (simple || isThermoFluxAttainable(min_flux, helper, potTest))) {
				min[a] = min_flux->getObjVal();
			}
			else {
				ScipModelPtr scip = factory.build(model);
				if(settings->timeout > 1) { // a timeout of less than a second makes no sense
					BOOST_SCIP_CALL( SCIPsetRealParam(scip->getScip(), "limits/time", settings->timeout) );
				}
				scip->setObjectiveSense(false);
				scip->solve();
#ifndef NDEBUG
				if(!scip->isOptimal()) {
					SCIPprintStatistics(scip->getScip(), NULL);
					cout << "LP primal infeasible = " << SCIPlpiIsPrimalInfeasible(min_flux->getLPI()) << endl;
					if( SCIPlpiHasDualRay(min_flux->getLPI())) {
						double dualfarkas[model->getMetabolites().size()];
						SCIPlpiGetDualfarkas(min_flux->getLPI(), dualfarkas);
						foreach(MetabolitePtr met, model->getMetabolites()) {
							double d = dualfarkas[min_flux->getIndex(met)];
							if(d < -prec->getDualFeasTol() || d > prec->getDualFeasTol()) {
								cout << met->getName() << " = " << d << endl;
							}
						}
					}
					else {
						cout << "no dual ray available" << endl;
					}
					SCIPwriteOrigProblem(scip->getScip(), "debug.lp", NULL, TRUE);
				}
#endif
				assert(scip->isOptimal());
				double opt = scip->getObjectiveValue();
				min[a] = opt;
#if REDUCE_DOMAIN
				a->setLb(opt); // we computed a lower bound, so use it for future computations
#endif
			}

			/*
			 * Actually, we are now adding constraints and not decrease the precision.
			 * This ok, since we solve with high dual precision.
			 * This means, the error that we do is rather in primal infeasibilities than in dual infeasibilities,
			 * i.e. the bounds that we compute tend to be a bit weaker.
			 * In particular it is unlikely that they are overtight and produce infeasibilities.
			 */
			double maxflux = max[a];
			double minflux = min[a];
			//double maxflux = maxf > minf ? maxf : minf;
			//double minflux = maxf > minf ? minf : maxf;
			max_flux->setUb(a, maxflux);
			min_flux->setUb(a, maxflux);
			max_flux->setLb(a, minflux);
			min_flux->setLb(a, minflux);
			max_flux->solveDual();
			min_flux->solveDual();
		}
		/*
		 * reset objective function
		 */
		max_flux->setObj(a,0);
		min_flux->setObj(a,0);

		a->setObj(0);

		/*
		 * Check, if we are still in the run time limit
		 */
		runningTime = (double) (clock() - start) / CLOCKS_PER_SEC;

		if(settings->timeout > 1 && runningTime > settings->timeout) {  // a timeout of less than a second makes no sense
			cout << endl;
			cout << "aborted by timeout of " << settings->timeout << " seconds" << endl;
			BOOST_THROW_EXCEPTION( TimeoutError() );
		}

		cout << "finished iteration " << i << " of " << num_rxns << endl;
		i++;
	}

	// reset precision
	model->setFluxPrecision(orig_precision);
}


} /* namespace metaopt */
