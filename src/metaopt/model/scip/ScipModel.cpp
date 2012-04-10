/*
 * ScipModel.cpp
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#include "ScipModel.h"
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include <iostream>

#include "metaopt/scip/ScipError.h"
#include "metaopt/Properties.h"

using namespace boost;

namespace metaopt {

ScipModel::ScipModel(ModelPtr model) {
	_model = model;

	BOOST_SCIP_CALL( init_scip() );
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

	return SCIP_OKAY;
}

SCIP_RETCODE ScipModel::free_scip() {
	SCIP_CALL( SCIPfree(&_scip) );
	return SCIP_OKAY;
}

ScipModel::~ScipModel() {

	// free scip
	if(free_scip() != SCIP_OKAY) {
		std::cerr << "Failed to free Scip, destroying anyways." << std::endl;
	}
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
				rxn->getName().c_str(),
				rxn->getLb(),
				rxn->getUb(),
				rxn->getObj(),
				SCIP_VARTYPE_CONTINUOUS,
				true,
				false,
				0,0,0,0,0) );
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
				met->getName().c_str(),
				met->getPotLb(),
				met->getPotUb(),
				met->getPotObj(),
				SCIP_VARTYPE_CONTINUOUS,
				true,
				false,
				0,0,0,0,0) );
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
	return(!_reactions.empty() && SCIPgetLPSolstat(_scip) == SCIP_LPSOLSTAT_OPTIMAL);
}

bool ScipModel::hasCurrentPotentials() {
	return(!_metabolites.empty() && SCIPgetLPSolstat(_scip) == SCIP_LPSOLSTAT_OPTIMAL);
}

double ScipModel::getCurrentFlux(ReactionPtr rxn) {
	if( SCIPgetLPSolstat(_scip) != SCIP_LPSOLSTAT_OPTIMAL) {
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("LP not solved to optimality") );
	}
	else {
		return SCIPgetVarSol(_scip, getFlux(rxn));
	}
}


double ScipModel::getCurrentPotential(MetabolitePtr met) {
	if( SCIPgetLPSolstat(_scip) != SCIP_LPSOLSTAT_OPTIMAL) {
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("LP not solved to optimality") );
	}
	else {
		return SCIPgetVarSol(_scip, getPotential(met));
	}
}


} /* namespace metaopt */
