/*
 * LPFlux.cpp
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#include "LPFlux.h"
#include "objscip/objscip.h"
#include "metaopt/Properties.h"
#include "metaopt/scip/ScipError.h"
#include <map>
#include <vector>
#include <iostream>

using namespace std;

namespace metaopt {

typedef pair<ReactionPtr, int> VarAssign;

LPFlux::LPFlux(ModelPtr model, bool exchange) {
	_model = model;
	BOOST_SCIP_CALL( init_lp(exchange) );
}

SCIP_RETCODE LPFlux::init_lp(bool exchange) {
	_lpi = NULL;
	SCIP_CALL( SCIPlpiCreate(&_lpi, "LPFlux", SCIP_OBJSEN_MAXIMIZE) );

	// create metabolite -> row_index map
	int metabolite_index = 0;
	foreach(const MetabolitePtr m, _model->getMetabolites()) {
		if(!exchange || !m->hasBoundaryCondition()) {
			// we
			_metabolites[m] = metabolite_index;
			metabolite_index++;
		}
	}
	_num_metabolites = metabolite_index;

	// create flux variables,
	// create reaction -> variable_index map
	// create stoichiometric matrix
	int reaction_var = 0;
	foreach(ReactionPtr r, _model->getReactions()) {
		vector<int> ind;
		vector<double> coef;
		if(exchange || !r->isExchange()) {
			_reactions[r] = reaction_var++;
			foreach(Stoichiometry m, r->getStoichiometries()) {
				if(!exchange || m.first->hasBoundaryCondition()) {
					ind.push_back(_metabolites[m.first]);
					coef.push_back(m.second);
				}
			}
			double obj = r->getObj();
			double lb = r->getLb();
			double ub = r->getUb();
			int beg = 0;
			// actually we have a nice name for the column, but it wants a char* instead of a const char*. I don't think it is worth copying names ;)
			SCIP_CALL( SCIPlpiAddCols(_lpi, 1, &obj, &lb, &ub, NULL, ind.size(), &beg, ind.data(), coef.data()) );
		}
	}
	_num_reactions = reaction_var;

	// set coefficients of rhs and lhs of every row to zero (steady state assupmtion)
	// metabolites with boundary condition are already excluded
	double zeros[_num_metabolites];
	int ind[_num_metabolites];
	for(int i = 0; i < _num_metabolites; i++) {
		ind[i] = i;
		zeros[i] = 0;
	}
	SCIP_CALL( SCIPlpiChgSides(_lpi, _num_metabolites, ind, zeros, zeros));

	return SCIP_OKAY;
}

LPFlux::~LPFlux() {
	int code = free_lp();
	// in case of error make sure that the object is destroyed anyways to keep harm as little as possible.
	if(code != SCIP_OKAY) {
		cerr << "Failed to free LPFlux, this may lead to a memory leak" << endl;
	}
}

SCIP_RETCODE LPFlux::free_lp() {
	SCIP_CALL( SCIPlpiFree(&_lpi) );
	return SCIP_OKAY;
}

void LPFlux::setLb(ReactionPtr r, double lb) {
	int ind = _reactions.at(r);
	double oub;
	BOOST_SCIP_CALL( SCIPlpiGetBounds(_lpi, ind, ind, NULL, &oub) );
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, 1, &ind, &lb, &oub) );
}

void LPFlux::setUb(ReactionPtr r, double ub) {
	int ind = _reactions.at(r);
	double olb;
	BOOST_SCIP_CALL( SCIPlpiGetBounds(_lpi, ind, ind, &olb, NULL) );
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, 1, &ind, &olb, &ub) );
}

void LPFlux::setObj(ReactionPtr r, double obj) {
	int ind = _reactions.at(r);
	BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, 1, &ind, &obj));
}

void LPFlux::setZeroObj() {
	double zeros[_num_reactions];
	int ind[_num_reactions];
	for(int i = 0; i < _num_reactions; i++) {
		ind[i] = i;
		zeros[i] = 0;
	}
	BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, _num_reactions, ind, zeros) );
}

void LPFlux::setBounds(LPFluxPtr other) {
	// the data may be stored in a different order, so we cannot simply do a batch copy, but have to translate the indices.
	// oind stores the desired permutation
	double olb[other->_num_reactions];
	double oub[other->_num_reactions];
	int oind[_num_reactions];
	foreach(VarAssign v, _reactions) {
		oind[v.second] = other->_reactions.at(v.first);
	}
	BOOST_SCIP_CALL( SCIPlpiGetBounds(other->_lpi, 0, other->_num_reactions-1, olb, oub) );
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	for(int i = 0; i < _num_reactions; i++) {
		ind[i] = i;
		lb[i] = olb[oind[i]];
		ub[i] = oub[oind[i]];
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setBounds(ScipModelPtr other) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		lb[i] = other->getCurrentFluxLb(v.first);
		ub[i] = other->getCurrentFluxUb(v.first);
		i++;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionBounds(LPFluxPtr flux) {
	// the data may be stored in a different order, so we cannot simply do a batch copy, but have to translate the indices.
	// oind stores the desired permutation
	double oflux[flux->_num_reactions];
	int oind[_num_reactions];
	foreach(VarAssign v, _reactions) {
		oind[v.second] = flux->_reactions.at(v.first);
	}
	BOOST_SCIP_CALL( SCIPlpiGetSol(flux->_lpi, NULL, oflux, NULL, NULL, NULL) );
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	for(int i = 0; i < _num_reactions; i++) {
		ind[i] = i;
		lb[i] = oflux[oind[i]] < -EPSILON ? -1 : 0;
		ub[i] = oflux[oind[i]] >  EPSILON ?  1 : 0;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionBounds(ScipModelPtr flux) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		lb[i] = flux->getCurrentFlux(v.first) < -EPSILON ? -1 : 0;
		ub[i] = flux->getCurrentFlux(v.first) >  EPSILON ?  1 : 0;
		i++;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionObj(LPFluxPtr flux) {
	// the data may be stored in a different order, so we cannot simply do a batch copy, but have to translate the indices.
	// oind stores the desired permutation
	double oflux[flux->_num_reactions];
	int oind[_num_reactions];
	foreach(VarAssign v, _reactions) {
		oind[v.second] = flux->_reactions.at(v.first);
	}
	BOOST_SCIP_CALL( SCIPlpiGetSol(flux->_lpi, NULL, oflux, NULL, NULL, NULL) );
	double obj[_num_reactions];
	int ind[_num_reactions];
	for(int i = 0; i < _num_reactions; i++) {
		ind[i] = i;
		double val = oflux[oind[i]];
		if(val < -EPSILON) obj[i] = -1;
		else if(val > EPSILON) obj[i] = 1;
		else obj[i] = 0;
	}
	BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, _num_reactions, ind, obj) );
}

void LPFlux::setDirectionObj(ScipModelPtr flux) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double obj[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		double val = flux->getCurrentFlux(v.first);
		if(val < -EPSILON) obj[i] = -1;
		else if(val > EPSILON) obj[i] = 1;
		else obj[i] = 0;
		i++;
	}
	BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, _num_reactions, ind, obj) );
}

void LPFlux::setObjSense(bool maximize) {
	if(maximize) {
		BOOST_SCIP_CALL( SCIPlpiChgObjsen(_lpi, SCIP_OBJSEN_MAXIMIZE) );
	}
	else {
		BOOST_SCIP_CALL( SCIPlpiChgObjsen(_lpi, SCIP_OBJSEN_MINIMIZE) );
	}
}

void LPFlux::solvePrimal() {
	BOOST_SCIP_CALL( SCIPlpiSolvePrimal(_lpi) );
}

void LPFlux::solveDual() {
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
}

void LPFlux::solve() {
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
}

double LPFlux::getFlux(ReactionPtr rxn) {
	boost::unordered_map<ReactionPtr, int>::iterator iter = _reactions.find(rxn);
	if(iter == _reactions.end()) {
		return 0;
	}
	else {
		double primsol[_num_reactions];
		BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, primsol, NULL, NULL, NULL) );
		return primsol[iter->second];
	}
}

double LPFlux::getDual(MetabolitePtr met) {
	boost::unordered_map<MetabolitePtr, int>::iterator iter = _metabolites.find(met);
	if(iter == _metabolites.end()) {
		return 0; // the constraint does not exist, hence it is never active, hence its dual value is always 0
	}
	else {
		double dualsol[_num_metabolites];
		BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, NULL, dualsol, NULL, NULL) );
		return dualsol[iter->second];
	}
}


bool LPFlux::isOptimal() {
	return SCIPlpiIsOptimal(_lpi);
}


bool LPFlux::isFeasible() {
	return SCIPlpiIsPrimalFeasible(_lpi);
}


} /* namespace metaopt */
