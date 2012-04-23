/*
 * CycleDeletionHeur.cpp
 *
 *  Created on: 23.04.2012
 *      Author: arnem
 */

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
	_difficultyTestFlux = LPFluxPtr(new LPFlux(_scip->getModel(), false));
	_difficultyTestFlux->setObjSense(_scip->isMaximize());

	_tflux = LPFluxPtr(new LPFlux(_scip->getModel(), true));
	// _tflux doesn't need an Objsense, since we do not use it to solve optimization problems.
	_cycle = LPFluxPtr(new LPFlux(_scip->getModel(), false));
	_cycle->setObjSense(true);
}

CycleDeletionHeur::~CycleDeletionHeur() {
	// nothing to do
}

SCIP_RETCODE CycleDeletionHeur::scip_exec(SCIP* scip, SCIP_HEUR* heur, SCIP_HEURTIMING timing, SCIP_RESULT* result) {
	if(!isDifficult()) {
		SCIP_SOL* raw_sol;
		SCIP_CALL( SCIPcreateOrigSol(scip, &raw_sol, heur));
		SolutionPtr sol = wrap(raw_sol, _scip);

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
		success = _scip->computeAddOnValues(sol);

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
	// first check if the LP was really solved to optimality
	// only if this is the case, run the heuristic
	if( _scip->hasCurrentFlux() ) {
		return true;
	}

	_difficultyTestFlux->setDirectionBounds(_scip); // use current LP-sol to set bounds

	_difficultyTestFlux->solveDual(); // we only change bounds, so solve using dual Simplex

	// zero flux is always possible, so only check feasibility for sanity
	assert( _difficultyTestFlux->isFeasible() );
	return _difficultyTestFlux->getObjVal() > EPSILON || _difficultyTestFlux->getObjVal() < -EPSILON;
}

bool CycleDeletionHeur::computeFluxSolution(SolutionPtr sol) {
	// use LP sol
	_tflux->set(_scip);

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
	foreach(ReactionPtr rxn, _scip->getModel()->getReactions()) {
		double val = _tflux->getFlux(rxn);
		SCIP_VAR* var = _scip->getFlux(rxn);
		BOOST_SCIP_CALL( SCIPsetSolVal(_scip->getScip(), sol.get(), var, val) );
	}

	// ignore flux forcing reactions for now
	// TODO: check to save computation on further processing
	return true;
}

bool CycleDeletionHeur::computePotentials(SolutionPtr sol) {

}

} /* namespace metaopt */
