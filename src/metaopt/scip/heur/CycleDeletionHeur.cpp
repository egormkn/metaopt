/*
 * CycleDeletionHeur.cpp
 *
 *  Created on: 23.04.2012
 *      Author: arnem
 */

#include <iostream>
#include "CycleDeletionHeur.h"
#include "metaopt/model/scip/Solution.h"
#include "metaopt/Properties.h"
#include "metaopt/scip/ScipError.h"

using namespace scip;
using namespace boost;

namespace metaopt {

CycleDeletionHeur::CycleDeletionHeur(ScipModelPtr scip) :
		ObjHeur(scip->getScip(), "CycleDeletionHeur",
				"compute solution by subtracting internal cycles",
				'C',
				-10000, 1, 0, -1, SCIP_HEURTIMING_AFTERLPLOOP, false) {
	_scip = scip;
	// initialize with original objective
	_difficultyTestFlux = LPFluxPtr(new LPFlux(scip->getModel(), false));
	_difficultyTestFlux->setObjSense(scip->isMaximize());

	_tflux = LPFluxPtr(new LPFlux(scip->getModel(), true));
	// _tflux doesn't need an Objsense, since we do not use it to solve optimization problems.
	_cycle = LPFluxPtr(new LPFlux(scip->getModel(), false));
	_cycle->setObjSense(true);

	_potentials = LPPotentialsPtr(new LPPotentials(scip->getModel()));
}

CycleDeletionHeur::~CycleDeletionHeur() {
	// nothing to do
}

SCIP_RETCODE CycleDeletionHeur::scip_exec(SCIP* scip, SCIP_HEUR* heur, SCIP_HEURTIMING timing, SCIP_RESULT* result) {
	if(!isDifficult()) {
		std::cout << "running heur" << std::endl;
		SCIP_SOL* raw_sol;
		SCIP_CALL( SCIPcreateOrigSol(scip, &raw_sol, heur));
		SolutionPtr sol = wrap(raw_sol, getScip());

		// use standard method to compute feasible solution
		bool success = computeFluxSolution(sol);

		if(!success) {
			//cout << "failed to find feasible flux in heuristic" << endl;
			*result=SCIP_DIDNOTFIND;
			return SCIP_OKAY;
		}

		// if we have potential variables, also compute potentials
		success = computePotentials(sol);

		if(!success) {
			*result=SCIP_DIDNOTFIND;
			return SCIP_OKAY;
		}

		// compute values for remaining addons
		success = getScip()->computeAddOnValues(sol);

		if(!success) {
			*result=SCIP_DIDNOTFIND;
			return SCIP_OKAY;
		}

		// set solution values
		unsigned int stored;

		// usually, we could also free the solution, but the wrapping does not allow this
		// TODO: is this a performance deficit?
		SCIP_CALL(SCIPtrySol(scip,sol.get(),
			TRUE,
			TRUE,
			TRUE,
			TRUE,
			&stored ));

		if(stored) {
			*result=SCIP_FOUNDSOL;
		}
		else {
			*result = SCIP_DIDNOTFIND;
		}
		return SCIP_OKAY;
	}
	else {
		*result = SCIP_DIDNOTRUN;
		return SCIP_OKAY;
	}
}

bool CycleDeletionHeur::isDifficult() {
	ScipModelPtr scip = getScip();
	// first check if the LP was really solved to optimality
	// only if this is the case, run the heuristic
	if( !scip->hasCurrentFlux() ) {
		return true;
	}

	_difficultyTestFlux->setDirectionBounds(scip); // use current LP-sol to set bounds

	_difficultyTestFlux->solveDual(); // we only change bounds, so solve using dual Simplex

	// zero flux is always possible, so only check feasibility for sanity
	assert( _difficultyTestFlux->isFeasible() );
	return _difficultyTestFlux->getObjVal() > EPSILON || _difficultyTestFlux->getObjVal() < -EPSILON;
}

bool CycleDeletionHeur::computeFluxSolution(SolutionPtr sol) {
	ScipModelPtr scip = getScip();
	// use LP sol
	_tflux->set(scip);

	bool hasFlux = true;

	// create new flux containing no exchange reactions, to find internal cycles
	_cycle->setDirectionObj(_tflux);
	// iteratively search for a cycle and subtract that cycle
	do {
		// only allow flux through reactions that still carry flux
		_cycle->setDirectionBounds(_tflux);
		_cycle->solve();
		assert(_cycle->isFeasible());
		hasFlux = _cycle->getObjVal() > EPSILON;
		if(hasFlux) {
			// now subtract computed flux from solution flux
			double scale = _tflux->getSubScale(_cycle);
			assert(scale > 0);
			_tflux->subtract(_cycle, scale);
		}
	}
	while(hasFlux);
	// copy computed solution into sol
	foreach(ReactionPtr rxn, scip->getModel()->getReactions()) {
		double val = _tflux->getFlux(rxn);
		SCIP_VAR* var = scip->getFlux(rxn);
		BOOST_SCIP_CALL( SCIPsetSolVal(scip->getScip(), sol.get(), var, val) );
	}

	// ignore flux forcing reactions for now
	// TODO: check to save computation on further processing
	return true;
}

inline double computePotDiff(LPPotentialsPtr& _potentials, ReactionPtr& rxn) {
	double val = 0;
	foreach(Stoichiometry s, rxn->getStoichiometries()) {
		val += s.second * _potentials->getPotential(s.first);
	}
	return val;
}

bool CycleDeletionHeur::perturb() {
	ScipModelPtr scip = getScip();
	foreach(ReactionPtr rxn, scip->getModel()->getReactions()) {
		if(!rxn->isExchange()) {
			double val = computePotDiff(_potentials, rxn);
			if(val >= 0) _potentials->setDirection(rxn, true);
			else _potentials->setDirection(rxn, false);
			if(val > -1 && val < 1) { //TODO: ugly
				// resolve
				_potentials->solveDual();
				if(!_potentials->isFeasible()) {
					return false;
				}
#ifndef NDEBUG
				val = computePotDiff(_potentials, rxn);
				assert(val < -1+EPSILON || val > 1-EPSILON);
#endif
			}
		}
	}
	return true;
}

bool CycleDeletionHeur::computePotentials(SolutionPtr sol) {
	ScipModelPtr scip = getScip();
	if(scip->hasPotentials()) {
		_potentials->setDirections(_tflux);
		_potentials->solveDual(); // we only changed bounds, so use dual simplex
		if(_potentials->isFeasible()) {
			perturb(); // make sure all potential differences are either <= -1 or >= 1
		}
		if(_potentials->isFeasible()) {
			// copy solution into sol
			foreach(MetabolitePtr met, scip->getModel()->getMetabolites()) {
				double val = _potentials->getPotential(met);
				if(scip->hasPotentialVar(met)) {
					SCIP_VAR* var = scip->getPotential(met);
					BOOST_SCIP_CALL( SCIPsetSolVal(scip->getScip(), sol.get(), var, val) );
				}
			}
			return true;
		}
		else {
			return false;
		}
	}
	else {
		// nothing to do
		return true;
	}
}

void createCycleDeletionHeur(ScipModelPtr scip) {
	// create insecure pointer, but thats ok since Scip will do all the allocation handling.
	// ofcourse, it gets more complicated, if we want to do more with the heur.
	CycleDeletionHeur* heur = new CycleDeletionHeur(scip);
	BOOST_SCIP_CALL( SCIPincludeObjHeur(scip->getScip(), heur, true) );
}

} /* namespace metaopt */
