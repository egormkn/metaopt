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
 * ScipModel.cpp
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#include "ScipModel.h"
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "metaopt/scip/ScipError.h"
#include "ModelAddOn.h"
#include "metaopt/Properties.h"

using namespace boost;
using namespace std;

namespace metaopt {


ScipModel::ScipModel(ModelPtr model) {
	_model = model;

	BOOST_SCIP_CALL( init_scip() );

	setPrecision(model->getFluxPrecision());
	setPotPrecision(model->getPotPrecision());
}

SCIP_RETCODE ScipModel::init_scip() {
	_scip = NULL;
	// create scip
	SCIP_CALL( SCIPcreate(&_scip) );

	// create a problem
	SCIP_CALL( SCIPcreateProb(_scip, "metaopt", 0,0,0,0,0,0,0) );

	// set obj sense to maximize
	SCIP_CALL( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );

	// include default plugins
	SCIP_CALL( SCIPincludeDefaultPlugins(_scip) );

	//SCIP_CALL( SCIPsetIntParam(_scip, "presolving/maxrounds",0) );

	return SCIP_OKAY;
}

SCIP_RETCODE ScipModel::free_scip() {
	SCIP_CALL( SCIPfree(&_scip) );
	return SCIP_OKAY;
}

ScipModel::~ScipModel() {
	// destruct addons first, because they may require a still-working ScipModel for destruction
	// destruct in reverse order, so that no dependencies are broken
	while(!_addons.empty()) {
		ModelAddOnPtr addon = _addons.back();
		addon->destroy(this);
		_addons.pop_back();
		assert(addon.unique()); // after the model is destroyed, the addon is useless
	}

	// free vars
	typedef std::pair<ReactionPtr, SCIP_VAR*> ReactionVar;
	typedef std::pair<MetabolitePtr, SCIP_VAR*> MetaboliteVar;
	foreach(ReactionVar r, _reactions) {
		SCIP_VAR* var = r.second;
		int code = SCIPreleaseVar(_scip, &var);
		assert(code == SCIP_OKAY);
	}
	_reactions.clear();

	foreach(MetaboliteVar r, _metabolites) {
		SCIP_VAR* var = r.second;
		int code = SCIPreleaseVar(_scip, &var);
		assert(code == SCIP_OKAY);
	}
	_metabolites.clear();

	// free scip
	int code = free_scip();
	assert(code == SCIP_OKAY);
}


void ScipModel::setPrecision(PrecisionPtr precision) {
	assert(SCIPgetStage(_scip) == SCIP_STAGE_INIT);
	_precision = precision;
	BOOST_SCIP_CALL( SCIPsetRealParam(_scip, "numerics/epsilon", _precision->getPrimalFeasTol()*1e-3) ); // absolute values smaller than this are considered zero by the solver
	BOOST_SCIP_CALL( SCIPsetRealParam(_scip, "numerics/sumepsilon", _precision->getPrimalFeasTol()) ); // absolute values of sums smaller than this are considered zero by the solver
	BOOST_SCIP_CALL( SCIPsetRealParam(_scip, "numerics/feastol", _precision->getPrimalFeasTol()) ); // feasibility tolerance for constraints
	BOOST_SCIP_CALL( SCIPsetRealParam(_scip, "numerics/lpfeastol", _precision->getPrimalFeasTol()) ); // primal feasibility tolerance of LP solver
	BOOST_SCIP_CALL( SCIPsetRealParam(_scip, "numerics/dualfeastol", _precision->getDualFeasTol()) ); // feasibility tolerance for reduced costs in LP solution
	BOOST_SCIP_CALL( SCIPsetRealParam(_scip, "numerics/barrierconvtol", _precision->getPrimalFeasTol()*1e-4) ); // LP convergence tolerance used in barrier algorithm (probably not used)
}

const PrecisionPtr& ScipModel::getPrecision() const {
	return _precision;
}

void ScipModel::setPotPrecision(PrecisionPtr precision) {
	assert(SCIPgetStage(_scip) == SCIP_STAGE_INIT);
	_potprecision = precision;
	// don't change scip, this is only important for implicit potential variables
}

const PrecisionPtr& ScipModel::getPotPrecision() const {
	return _potprecision;
}



SCIP_VAR* ScipModel::getFlux(ReactionPtr rxn) {
	unordered_map<ReactionPtr, SCIP_VAR*>::iterator iter = _reactions.find(rxn);
	if(iter != _reactions.end()) {
		return iter->second;
	}
	else {
		SCIP_VAR* var = NULL;
		BOOST_SCIP_CALL( SCIPcreateVar(_scip,
				&var,
				rxn->getCName(),
				rxn->getLb(),
				rxn->getUb(),
				rxn->getObj(),
				SCIP_VARTYPE_CONTINUOUS,
				true,
				false,
				0,0,0,0,0) );

		BOOST_SCIP_CALL( SCIPaddVar(_scip, var) );

		BOOST_SCIP_CALL( SCIPmarkDoNotMultaggrVar(_scip, var) ); //TODO: Is this really necessary?

		_reactions[rxn] = var;
		return var;
	}
}

SCIP_VAR* ScipModel::getPotential(MetabolitePtr met) {
	unordered_map<MetabolitePtr, SCIP_VAR*>::iterator iter = _metabolites.find(met);
	if(iter != _metabolites.end()) {
		return iter->second;
	}
	else {
		SCIP_VAR* var = NULL;
		BOOST_SCIP_CALL( SCIPcreateVar(_scip,
				&var,
				met->getCName(),
				met->getPotLb(),
				met->getPotUb(),
				met->getPotObj(),
				SCIP_VARTYPE_CONTINUOUS,
				true,
				false,
				0,0,0,0,0) );

		BOOST_SCIP_CALL( SCIPaddVar(_scip, var) );

		BOOST_SCIP_CALL( SCIPmarkDoNotMultaggrVar(_scip, var) ); //TODO: Is this really necessary?

		_metabolites[met] = var;
		return var;
	}
}

double ScipModel::getCurrentFluxUb(ReactionPtr rxn) {
	unordered_map<ReactionPtr, SCIP_VAR*>::iterator iter = _reactions.find(rxn);
	if(iter != _reactions.end()) {
		return SCIPvarGetUbLocal(iter->second);
	}
	else {
		return rxn->getUb();
	}
}

double ScipModel::getCurrentFluxLb(ReactionPtr rxn) {
	unordered_map<ReactionPtr, SCIP_VAR*>::iterator iter = _reactions.find(rxn);
	if(iter != _reactions.end()) {
		return SCIPvarGetLbLocal(iter->second);
	}
	else {
		return rxn->getLb();
	}
}

double ScipModel::getCurrentPotentialUb(MetabolitePtr met) {
	unordered_map<MetabolitePtr, SCIP_VAR*>::iterator iter = _metabolites.find(met);
	if(iter != _metabolites.end()) {
		return SCIPvarGetUbLocal(iter->second);
	}
	else {
		return met->getPotUb();
	}
}

double ScipModel::getCurrentPotentialLb(MetabolitePtr met) {
	unordered_map<MetabolitePtr, SCIP_VAR*>::iterator iter = _metabolites.find(met);
	if(iter != _metabolites.end()) {
		return SCIPvarGetLbLocal(iter->second);
	}
	else {
		return met->getPotLb();
	}
}

bool ScipModel::hasFluxVar(ReactionPtr rxn) {
	unordered_map<ReactionPtr, SCIP_VAR*>::iterator iter = _reactions.find(rxn);
	return iter != _reactions.end();
}

bool ScipModel::hasPotentialVar(MetabolitePtr met) {
	unordered_map<MetabolitePtr, SCIP_VAR*>::iterator iter = _metabolites.find(met);
	return iter != _metabolites.end();
}

bool ScipModel::hasCurrentFlux() {
	return(!_reactions.empty()
		&& (
				(SCIPgetStage(_scip) == SCIP_STAGE_SOLVING && (SCIPgetLPSolstat(_scip) == SCIP_LPSOLSTAT_OPTIMAL || SCIPhasPrimalRay(_scip)) )
			 || SCIPgetStage(_scip) == SCIP_STAGE_SOLVED
		   )
		  );
}

bool ScipModel::hasCurrentPotentials() {
	return(!_metabolites.empty()
		&& (
				(SCIPgetStage(_scip) == SCIP_STAGE_SOLVING && SCIPgetLPSolstat(_scip) == SCIP_LPSOLSTAT_OPTIMAL)
			 || SCIPgetStage(_scip) == SCIP_STAGE_SOLVED
		   )
		  );
}

bool ScipModel::hasFlux() {
	return(!_reactions.empty() );
}

bool ScipModel::hasPotentials() {
	return(!_metabolites.empty() );
}

double ScipModel::getCurrentFlux(ReactionPtr rxn) {
	assert( hasFluxVar(rxn) );
	if( SCIPgetStage(_scip) == SCIP_STAGE_SOLVING) {
		if(SCIPhasPrimalRay(_scip)) {
			// if it is unbounded, set value of ray instead
			return SCIPgetPrimalRayVal(_scip, getFlux(rxn));
		}
		else {
			assert( SCIPgetLPSolstat(_scip) == SCIP_LPSOLSTAT_OPTIMAL ); // LP must be solved to optimality, else the result is rather meaningless
			/*if(oldnode == SCIPgetCurrentNode(_scip)) {
				// return cached value
			}*/
			return SCIPgetSolVal(_scip, NULL, getFlux(rxn));
		}
	}
	else if(SCIPgetStage(_scip) == SCIP_STAGE_SOLVED) {
		if(SCIPhasPrimalRay(_scip)) {
			return SCIPgetPrimalRayVal(_scip, getFlux(rxn));
		}
		else {
			SCIP_SOL* sol = SCIPgetBestSol(_scip);
			return SCIPgetSolVal(_scip, sol, getFlux(rxn));
		}
	}
	else {
		assert( false ); // solving has not yet started
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("Solving has not yet started!") );
		return 0; // never called
	}
}


double ScipModel::getCurrentPotential(MetabolitePtr met) {
	assert( hasPotentialVar(met) );
	if( SCIPgetStage(_scip) == SCIP_STAGE_SOLVING) {
		assert( SCIPgetLPSolstat(_scip) == SCIP_LPSOLSTAT_OPTIMAL ); // LP must be solved to optimality, else the result is rather meaningless
		return SCIPgetSolVal(_scip, NULL, getPotential(met));
	}
	else if(SCIPgetStage(_scip) == SCIP_STAGE_SOLVED) {
		SCIP_SOL* sol = SCIPgetBestSol(_scip);
		return SCIPgetSolVal(_scip, sol, getPotential(met));
	}
	else {
		assert( false ); // solving has not yet started
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("Solving has not yet started!") );
		return 0; // never called
	}
}


double ScipModel::getFlux(SolutionPtr sol, ReactionPtr rxn) {
	assert( hasFluxVar(rxn) );
	return SCIPgetSolVal(_scip, sol.get(), getFlux(rxn));
}


double ScipModel::getPotential(SolutionPtr sol, MetabolitePtr met) {
	assert( hasPotentialVar(met) );
	return SCIPgetSolVal(_scip, sol.get(), getPotential(met));
}

void ScipModel::solve() {
	BOOST_SCIP_CALL( SCIPsolve(_scip) );
}

double ScipModel::getObjectiveValue() {
	assert( SCIPgetStage(_scip) == SCIP_STAGE_SOLVED );  // Problem has not yet been solved!";
	SCIP_SOL* sol = SCIPgetBestSol(_scip);
	return SCIPgetSolOrigObj(_scip, sol);
}

bool ScipModel::isOptimal() {
	assert( SCIPgetStage(_scip) == SCIP_STAGE_SOLVED );  // Problem has not yet been solved!";
	return SCIPgetStatus(_scip) == SCIP_STATUS_OPTIMAL;
}

bool ScipModel::isUnbounded() {
	assert( SCIPgetStage(_scip) == SCIP_STAGE_SOLVED );  // Problem has not yet been solved!";
	return SCIPgetStatus(_scip) == SCIP_STATUS_UNBOUNDED;
}

bool ScipModel::hasPrimalRay() {
	return SCIPhasPrimalRay(_scip);
}

bool ScipModel::isInfeasible() {
	assert( SCIPgetStage(_scip) == SCIP_STAGE_SOLVED );  // Problem has not yet been solved!";
	return SCIPgetStatus(_scip) == SCIP_STATUS_INFEASIBLE;
}

void ScipModel::addAddOn(ModelAddOnPtr addon) {
	_addons.push_back(addon);
}

bool ScipModel::computeAddOnValues(SolutionPtr sol) {
	bool success = true;
	for(vector<ModelAddOnPtr>::iterator iter = _addons.begin(); success && iter != _addons.end(); iter++) {
		success = (*iter)->computeSolutionVals(sol);
	}
	return success;
}

void ScipModel::setDirection(SCIP_NODE* node, ReactionPtr rxn, bool fwd) {
	assert( hasFluxVar(rxn) );
	// set direction on the flux vars
	SCIP_VAR* var = getFlux(rxn);
	if(fwd) {
		BOOST_SCIP_CALL( SCIPchgVarLbNode(_scip, node, var, 0) ); // exclude negative flux
	}
	else {
		BOOST_SCIP_CALL( SCIPchgVarUbNode(_scip, node, var, 0) ); // exclude positive flux
	}
	// notify addons
	for(vector<ModelAddOnPtr>::iterator iter = _addons.begin(); iter != _addons.end(); iter++) {
		(*iter)->setDirection(node, rxn, fwd);
	}
}

void ScipModel::setBlockedFlux(SCIP_NODE* node, ReactionPtr rxn, bool fwd) {
	assert( hasFluxVar(rxn) );
	// set direction on the flux vars
	SCIP_VAR* var = getFlux(rxn);
	if(fwd) {
		// exclude positive flux
		if(node == NULL) {
			SCIPchgVarUb(_scip, var, 0); // this works in every stage
		}
		else {
			BOOST_SCIP_CALL( SCIPchgVarUbNode(_scip, node, var, 0) );
		}
	}
	else {
		// exclude negative flux
		if(node == NULL) {
			SCIPchgVarLb(_scip, var, 0);
		}
		else {
			BOOST_SCIP_CALL( SCIPchgVarLbNode(_scip, node, var, 0) );
		}
	}
}

shared_ptr<unordered_set<ReactionPtr> > ScipModel::getFixedDirections() {
	shared_ptr<unordered_set<ReactionPtr> > result(new unordered_set<ReactionPtr>());
	for(vector<ModelAddOnPtr>::iterator iter = _addons.begin(); iter != _addons.end(); iter++) {
		shared_ptr<const unordered_set<ReactionPtr> > fixed = (*iter)->getFixedDirections();
		result->insert(fixed->begin(), fixed->end());
	}

	// alternatively the model may also know about fixed reactions
	// But, this is currently only implemented in a very rudimentary way, since we would need to track this information in the nodes.
	// And I do not know how to store user-data in the search-tree nodes.
	// Ofcourse it would be possible to store a hash-map that performs the lookup,
	// however, this solution would not know about nodes that have already been cut off.

	// only, if we know that a reaction is flux forcing, it is fixed.
	//TODO: change, if we create positive upper bounds etc.
	const unordered_set<ReactionPtr>& fluxforcing = getModel()->getFluxForcingReactions();
	result->insert(fluxforcing.begin(), fluxforcing.end());

	return result;
}


} /* namespace metaopt */
