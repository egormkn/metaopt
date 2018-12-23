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
 * ReducedScipFluxModel.cpp
 *
 *  Created on: 04.07.2012
 *      Author: arnem
 */

#include "ReducedScipFluxModel.h"
#include "Properties.h"

namespace metaopt {

ReducedScipFluxModel::ReducedScipFluxModel(ScipModelPtr origin)
	: _origin(origin)
{
	//nothing
}

ReducedScipFluxModel::~ReducedScipFluxModel() {
	// nothing
}

bool ReducedScipFluxModel::isMaximize() const {
	//ScipModelPtr scip = getScip();
	return false; // the transformed problem is always a minimization problem.
}

SCIP_VAR* ReducedScipFluxModel::getFlux(ReactionPtr rxn) {
	assert(hasFluxVar(rxn));
	return _reactions.at(rxn);
}


void ReducedScipFluxModel::setFlux(ReactionPtr rxn, SCIP_VAR* var) {
	_reactions[rxn] = var;
}

double ReducedScipFluxModel::getCurrentFluxUb(ReactionPtr rxn) {
	assert(hasFluxVar(rxn));
	return SCIPvarGetUbLocal(_reactions[rxn]);
}

double ReducedScipFluxModel::getCurrentFluxLb(ReactionPtr rxn) {
	assert(hasFluxVar(rxn));
	return SCIPvarGetLbLocal(_reactions[rxn]);
}

bool ReducedScipFluxModel::hasCurrentFlux() {
	ScipModelPtr scip = getScip();
	return scip->hasCurrentFlux();
}

double ReducedScipFluxModel::getCurrentFlux(ReactionPtr rxn) {
	assert(hasFluxVar(rxn));
	return SCIPvarGetLbLocal(_reactions[rxn]);
}

double ReducedScipFluxModel::getFlux(SolutionPtr sol, ReactionPtr rxn) {
	ScipModelPtr scip = getScip();
	if( SCIPgetStage(scip->getScip()) == SCIP_STAGE_SOLVING) {
		assert( SCIPgetLPSolstat(scip->getScip()) == SCIP_LPSOLSTAT_OPTIMAL ); // LP must be solved to optimality, else the result is rather meaningless
		return SCIPgetSolVal(scip->getScip(), NULL, getFlux(rxn));
	}
	else if(SCIPgetStage(scip->getScip()) == SCIP_STAGE_SOLVED) {
		SCIP_SOL* sol = SCIPgetBestSol(scip->getScip());
		return SCIPgetSolVal(scip->getScip(), sol, getFlux(rxn));
	}
	else {
		assert( false ); // solving has not yet started
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("Solving has not yet started!") );
		return 0; // never called
	}
}

bool ReducedScipFluxModel::hasFluxVar(ReactionPtr rxn) {
	return _reactions.find(rxn) != _reactions.end();
}

bool ReducedScipFluxModel::hasFlux() {
	return(!_reactions.empty() );
}


} /* namespace metaopt */
