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
	// use default precision for _difficultyTestFlux, because it uses only -1/0/1 bounds.
	_difficultyTestFlux->setObjSense(scip->isMaximize());

	_tflux = LPFluxPtr(new LPFlux(scip->getModel(), true));
	_tflux->setPrecision(scip->getPrecision()); // _tflux eventually produces a solution for the original problem, so it should have the same precision
	// _tflux doesn't need an Objsense, since we do not use it to solve optimization problems.
	_cycle = LPFluxPtr(new LPFlux(scip->getModel(), false));
	_cycle->setPrecision(scip->getPrecision()->getSlavePrecision()); // solve _cycle with slave precision, because we will subtract it several times and errors may accumulate
	_cycle->setObjSense(true);

	_potentials = LPPotentialsPtr(new LPPotentials(scip->getModel()));
	_potentials->setPrecision(scip->getPotPrecision());
}

CycleDeletionHeur::~CycleDeletionHeur() {
	// nothing to do
}

SCIP_RETCODE CycleDeletionHeur::scip_exec(SCIP* scip, SCIP_HEUR* heur, SCIP_HEURTIMING timing, SCIP_RESULT* result) {
	if(!isDifficult()) {
		//std::cout << "running heur" << std::endl;
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
	if( !_difficultyTestFlux->isFeasible() ) {
		// numerical issues may cause this - we then don't want to fail completely
		return true; // if we have numerical issues, the problem is difficult!
	}

	const PrecisionPtr& prec = _difficultyTestFlux->getPrecision();

	return _difficultyTestFlux->getObjVal() > prec->getCheckTol() || _difficultyTestFlux->getObjVal() < -prec->getCheckTol();
}

bool CycleDeletionHeur::computeFluxSolution(SolutionPtr sol) {
	ScipModelPtr scip = getScip();
	// use LP sol
	_tflux->set(scip);

	bool hasFlux = true;

	const PrecisionPtr& cyclePrec = _cycle->getPrecision();

	// create new flux containing no exchange reactions, to find internal cycles
	_cycle->setDirectionObj(_tflux);
	// iteratively search for a cycle and subtract that cycle
#ifndef NDEBUG
	int debugi = 0;
#endif
	do {
		// only allow flux through reactions that still carry flux
		_cycle->setDirectionBounds(_tflux);
		_cycle->solve();
		if(!_cycle->isFeasible()) {
			// usually, this should never happen, but we may run into numerical issues
			// these numerical issues should not kill the program
			return false;
		}
		hasFlux = _cycle->getObjVal() > cyclePrec->getCheckTol();
		if(hasFlux) {
			// now subtract computed flux from solution flux
			double scale = _tflux->getSubScale(_cycle);
			assert(scale > 0);
			_tflux->subtract(_cycle, scale);
		}
#ifndef NDEBUG
		debugi++;
		if(debugi % 100 == 0) std::cout << "CycleDeletionHeur.cpp " << __LINE__ << " caught endless loop" << std::endl;
#endif
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

/*bool CycleDeletionHeur::perturb() {
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
}*/

bool CycleDeletionHeur::computePotentials(SolutionPtr sol) {
	ScipModelPtr scip = getScip();
	if(scip->hasPotentials()) {
		_potentials->setDirections(_tflux);
		// first check if we can find an optimal solution
		bool feasible;
		if(!_potentials->testStrictFeasible(feasible)) {
			// we failed to determine feasibility status, so be pessimistic and return false
			return false;
		}
		if(!feasible) return false;
		_potentials->optimize(); // we don't care if we failed to compute an optimal solution.
		// it suffices if the computed solution is feasible
		if(!_potentials->isCurrentSolutionFeasible()) {
			return false;
		}
		// copy solution into sol
		foreach(MetabolitePtr met, scip->getModel()->getMetabolites()) {
			double val = _potentials->getPotential(met);
			if(scip->hasPotentialVar(met)) {
				SCIP_VAR* var = scip->getPotential(met);
				BOOST_SCIP_CALL( SCIPsetSolVal(scip->getScip(), sol.get(), var, val) );
			}
		}
#ifndef NDEBUG
		const PrecisionPtr& fluxPrec = _tflux->getPrecision();
		const PrecisionPtr& potPrec = _potentials->getPrecision();
		// check if potential differences match to flux directions
		foreach(ReactionPtr rxn, scip->getModel()->getInternalReactions()) {
			double potDiff = 0;
			foreach(Stoichiometry s, rxn->getStoichiometries()) {
				potDiff += s.second*_potentials->getPotential(s.first);
			}
			double flux = _tflux->getFlux(rxn);
			if(!(flux*potDiff < 0 || (flux < fluxPrec->getCheckTol() && flux > -fluxPrec->getCheckTol()) || (potDiff < potPrec->getCheckTol() && potDiff > -potPrec->getCheckTol()))) { // also allow potDiff = 0, because of closure of the domain
				_potentials->save();
				foreach(Stoichiometry s, rxn->getStoichiometries()) {
					std::cout << s.first->getName() << " ("<< _potentials->getVar(s.first) << ") : " << s.second << "*" << _potentials->getPotential(s.first) << std::endl;
				}
				std::cout << rxn->getName() << ": " << _potentials->getCon(rxn) << std::endl;
				assert(false);
			}
		}
#endif
		return true;
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
