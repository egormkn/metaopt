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
 * LPFlux.cpp
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#include "metaopt/Properties.h"

#include <map>
#include <vector>
#include <iostream>
#include "LPFlux.h"
#include "objscip/objscip.h"
#include "metaopt/scip/ScipError.h"

using namespace std;
using namespace boost;

namespace metaopt {

typedef pair<ReactionPtr, int> VarAssign;

LPFlux::LPFlux(ModelPtr model, bool exchange) {
	_model = model;
	BOOST_SCIP_CALL( init_lp(exchange) );
}

SCIP_RETCODE LPFlux::init_lp(bool exchange) {
	_lpi = NULL;
	SCIP_CALL( SCIPlpiCreate(&_lpi, NULL, "LPFlux", SCIP_OBJSEN_MAXIMIZE) );
	SCIPlpiSetRealpar(_lpi, SCIP_LPPAR_FEASTOL, 1e-10);
	SCIPlpiSetRealpar(_lpi, SCIP_LPPAR_DUALFEASTOL, 1e-10);
	SCIPlpiSetIntpar(_lpi, SCIP_LPPAR_PRESOLVING, 0);

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

	// create (empty) rows
	{
		// the new scip version seems to require that rows are created beforehand
		double lhs = 0;
		double rhs = 0;
		int beg[0];
		int ind[0];
		double var[0];
		for(int i = 0; i < _num_metabolites; i++) {
			SCIP_CALL( SCIPlpiAddRows(_lpi, 1, &lhs, &rhs, NULL, 0, beg, ind, var) );
		}
	}


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
				if(!exchange || !m.first->hasBoundaryCondition()) {
					ind.push_back(_metabolites.at(m.first));
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
	_primsol.resize(_num_reactions,0); // allocate sufficient memory
	_cstat_computed = false;
	_redcost_computed = false;

	int nrows;
	SCIP_CALL( SCIPlpiGetNRows(_lpi, &nrows) );

#ifndef NDEBUG
	SCIP_CALL( SCIPlpiGetNRows(_lpi, &nrows) );
	assert(nrows == _num_metabolites);
#endif

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
	assert(code == SCIP_OKAY);
}

SCIP_RETCODE LPFlux::free_lp() {
	SCIP_CALL( SCIPlpiFree(&_lpi) );
	return SCIP_OKAY;
}

void LPFlux::setLb(ReactionPtr r, double lb) {
	try {
		int ind = _reactions.at(r);
		double oub;
		BOOST_SCIP_CALL( SCIPlpiGetBounds(_lpi, ind, ind, NULL, &oub) );
		BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, 1, &ind, &lb, &oub) );
	}
	catch(std::exception &ex) {
		BOOST_SCIP_CALL( SCIP_ERROR );
	}
}

void LPFlux::setUb(ReactionPtr r, double ub) {
	try {
		int ind = _reactions.at(r);
		double olb;
		BOOST_SCIP_CALL( SCIPlpiGetBounds(_lpi, ind, ind, &olb, NULL) );
		BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, 1, &ind, &olb, &ub) );
	}
	catch(std::exception &ex) {
		BOOST_SCIP_CALL( SCIP_ERROR );
	}
}

void LPFlux::setObj(ReactionPtr r, double obj) {
	try {
		int ind = _reactions.at(r);
		double olb;
		BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, 1, &ind, &obj));
	}
	catch(std::exception &ex) {
		BOOST_SCIP_CALL( SCIP_ERROR );
	}
}

double LPFlux::getLb(ReactionPtr r) {
	try {
		int ind = _reactions.at(r);
		double lb;
		BOOST_SCIP_CALL( SCIPlpiGetBounds(_lpi, ind, ind, &lb, NULL) );
		return lb;
	}
	catch(std::exception &ex) {
		BOOST_SCIP_CALL( SCIP_ERROR );
		return NAN;
	}
}

double LPFlux::getUb(ReactionPtr r) {
	try {
		int ind = _reactions.at(r);
		double ub;
		BOOST_SCIP_CALL( SCIPlpiGetBounds(_lpi, ind, ind, NULL, &ub) );
		return ub;
	}
	catch(std::exception &ex) {
		BOOST_SCIP_CALL( SCIP_ERROR );
		return NAN;
	}
}

double LPFlux::getObj(ReactionPtr r) {
	try {
		int ind = _reactions.at(r);
		double obj;
		BOOST_SCIP_CALL( SCIPlpiGetObj(_lpi, ind, ind, &obj) );
		return obj;
	}
	catch(std::exception &ex) {
		BOOST_SCIP_CALL( SCIP_ERROR );
		return NAN;
	}
}

void LPFlux::setObjective(LPFluxPtr other) {
	/**
	 * reallign values
	 */
	double oobj[other->_num_reactions];
	BOOST_SCIP_CALL( SCIPlpiGetObj(other->_lpi, 0, other->_num_reactions-1, oobj) );
	int ind[_num_reactions];
	double obj[_num_reactions];
	foreach(VarAssign v, _reactions) {
		ind[v.second] = v.second;
		obj[v.second] = oobj[other->_reactions.at(v.first)];
	}
	// set
	BOOST_SCIP_CALL(SCIPlpiChgObj(_lpi, _num_reactions, ind, obj));
}


bool LPFlux::isMaximize() {
	SCIP_OBJSEN objsen;
	BOOST_SCIP_CALL(SCIPlpiGetObjsen(_lpi, &objsen));
	return objsen == SCIP_OBJSEN_MAXIMIZE;
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

void LPFlux::setBounds(AbstractScipFluxModelPtr other) {
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
	vector<double>& oflux = flux->_primsol;
	int oind[_num_reactions];
	foreach(VarAssign v, _reactions) {
		oind[v.second] = flux->_reactions.at(v.first);
	}
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

void LPFlux::setDirectionBounds(AbstractScipFluxModelPtr flux) {
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

void LPFlux::setDirectionBounds(SolutionPtr sol, AbstractScipFluxModelPtr flux) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		lb[i] = flux->getFlux(sol, v.first) < -EPSILON ? -1 : 0;
		ub[i] = flux->getFlux(sol, v.first) >  EPSILON ?  1 : 0;
		i++;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionBoundsInfty(LPFluxPtr flux) {
	// the data may be stored in a different order, so we cannot simply do a batch copy, but have to translate the indices.
	// oind stores the desired permutation
	vector<double>& oflux = flux->_primsol;
	int oind[_num_reactions];
	foreach(VarAssign v, _reactions) {
		oind[v.second] = flux->_reactions.at(v.first);
	}
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	for(int i = 0; i < _num_reactions; i++) {
		ind[i] = i;
		lb[i] = oflux[oind[i]] < -EPSILON ? -INFINITY : 0;
		ub[i] = oflux[oind[i]] >  EPSILON ?  INFINITY : 0;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionBoundsInfty(AbstractScipFluxModelPtr flux) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		lb[i] = flux->getCurrentFlux(v.first) < -EPSILON ? -INFINITY : 0;
		ub[i] = flux->getCurrentFlux(v.first) >  EPSILON ?  INFINITY : 0;
		i++;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionBoundsInfty(SolutionPtr sol, AbstractScipFluxModelPtr flux) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double lb[_num_reactions];
	double ub[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		assert(flux->getCurrentFluxLb(v.first)-EPSILON < flux->getFlux(sol, v.first) && flux->getCurrentFluxUb(v.first)+EPSILON > flux->getFlux(sol, v.first));
		lb[i] = flux->getFlux(sol, v.first) < -EPSILON ? -INFINITY : 0;
		ub[i] = flux->getFlux(sol, v.first) >  EPSILON ?  INFINITY : 0;
		i++;
	}
	BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_reactions, ind, lb, ub) );
}

void LPFlux::setDirectionObj(LPFluxPtr flux) {
	// the data may be stored in a different order, so we cannot simply do a batch copy, but have to translate the indices.
	// oind stores the desired permutation
	vector<double>& oflux = flux->_primsol;
	int oind[_num_reactions];
	foreach(VarAssign v, _reactions) {
		oind[v.second] = flux->_reactions.at(v.first);
	}
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

void LPFlux::setDirectionObj(AbstractScipFluxModelPtr flux) {
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

void LPFlux::setDirectionObj(SolutionPtr sol, AbstractScipFluxModelPtr flux) {
	// Here, we have no other option than to iterate through all reactions and do a seperate function call to get the bounds
	double obj[_num_reactions];
	int ind[_num_reactions];
	int i = 0;
	foreach(VarAssign v, _reactions) {
		ind[i] = v.second;
		double val = flux->getFlux(sol, v.first);
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
	_cstat_computed = false;
	_redcost_computed = false;
	BOOST_SCIP_CALL( SCIPlpiSolvePrimal(_lpi) );
#ifdef LPSOLVER_SOPLEX
	if(SCIPlpiGetInternalStatus(_lpi) == -4) {
		// basis is singular, we should resolve from scratch
		BOOST_SCIP_CALL( SCIPlpiClearState(_lpi) );
		// resolve
		BOOST_SCIP_CALL( SCIPlpiSolvePrimal(_lpi) );
	}
#endif
	if(SCIPlpiIsPrimalFeasible(_lpi)) {
		// capture solution in _primsol
		BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
	}
	if(!SCIPlpiIsStable(_lpi)) {
		cout << "warning unstable solution" << endl;
	}
	if(!SCIPlpiIsOptimal(_lpi)) {
		cout << "warning not solved to optimality" << endl;
	}
}

void LPFlux::solveDual() {
	_cstat_computed = false;
	_redcost_computed = false;
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
#ifdef LPSOLVER_SOPLEX
	if(SCIPlpiGetInternalStatus(_lpi) == -4) {
		// basis is singular, we should resolve from scratch
		BOOST_SCIP_CALL( SCIPlpiClearState(_lpi) );
		// resolve
		BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
	}
#endif
	if(SCIPlpiIsPrimalFeasible(_lpi)) {
		// capture solution in _primsol
		BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
	}
	if(!SCIPlpiIsStable(_lpi)) {
		cout << "warning unstable solution" << endl;
	}
	if(!SCIPlpiIsOptimal(_lpi)) {
		cout << "warning not solved to optimality" << endl;
	}
}

void LPFlux::solve() {
	_cstat_computed = false;
	_redcost_computed = false;
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
#ifdef LPSOLVER_SOPLEX
	if(SCIPlpiGetInternalStatus(_lpi) == -4) {
		// basis is singular, we should resolve from scratch
		BOOST_SCIP_CALL( SCIPlpiClearState(_lpi) );
		// resolve
		BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
	}
#endif
	if(SCIPlpiIsPrimalFeasible(_lpi)) {
		// capture solution in _primsol
		BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
	}
	if(!SCIPlpiIsStable(_lpi)) {
		cout << "warning unstable solution" << endl;
	}
	if(!SCIPlpiIsOptimal(_lpi)) {
		cout << "warning not solved to optimality" << endl;
	}
}

void LPFlux::resetState() {
	BOOST_SCIP_CALL( SCIPlpiClearState(_lpi) );
}

double LPFlux::getFlux(ReactionPtr rxn) {
	boost::unordered_map<ReactionPtr, int>::iterator iter = _reactions.find(rxn);
	if(iter == _reactions.end()) {
		return 0;
	}
	else {
		return _primsol[iter->second];
	}
}

double LPFlux::getAlpha(ReactionPtr rxn) {
	return getFlux(rxn);
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

int LPFlux::getColumnStatus(ReactionPtr rxn) {
	if(!_cstat_computed) {
		int ncols;
		BOOST_SCIP_CALL( SCIPlpiGetNCols(_lpi, &ncols) );
		_cstat.resize(ncols, 0);
		BOOST_SCIP_CALL( SCIPlpiGetBase(_lpi, _cstat.data(), NULL) );
		_cstat_computed = true;
	}
	boost::unordered_map<ReactionPtr, int>::iterator iter = _reactions.find(rxn);
	if(iter == _reactions.end()) {
		return SCIP_BASESTAT_ZERO;
	}
	else {
		return _cstat[iter->second];
	}
}

double LPFlux::getReducedCost(ReactionPtr rxn) {
	if(!_redcost_computed) {
		int ncols;
		BOOST_SCIP_CALL( SCIPlpiGetNCols(_lpi, &ncols) );
		_redcost.resize(ncols, 0);
		BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, NULL, NULL, NULL, _redcost.data()) );
		_redcost_computed = true;
	}
	boost::unordered_map<ReactionPtr, int>::iterator iter = _reactions.find(rxn);
	if(iter == _reactions.end()) {
		return SCIP_BASESTAT_ZERO;
	}
	else {
		return _redcost[iter->second];
	}
}


void LPFlux::set(AbstractScipFluxModelPtr smodel) {
	foreach(VarAssign v, _reactions) {
		_primsol[v.second] = smodel->getCurrentFlux(v.first);
	}
}

void LPFlux::set(SolutionPtr sol, AbstractScipFluxModelPtr smodel) {
	foreach(VarAssign v, _reactions) {
		_primsol[v.second] = smodel->getFlux(sol, v.first);
	}
}

void LPFlux::subtract(LPFluxPtr flux, double scale) {
	foreach(VarAssign v, flux->_reactions) {
		_primsol[_reactions.at(v.first)] -= flux->_primsol[v.second] * scale;
	}
}

double LPFlux::getSubScale(LPFluxPtr source) {
	double scale = 10000000;
	bool found = false;
	foreach(VarAssign v, _reactions) {
		double val = _primsol[v.second];
		double subVal = source->getFlux(v.first);
		if(val > EPSILON && subVal > EPSILON ) {
			if(scale > val/subVal) {
				scale = val/subVal;
				found = true;
			}
		}
		else if(val < -EPSILON && subVal < -EPSILON ) {
			if(scale > val/subVal) {
				scale = val/subVal;
				found = true;
			}
		}
	}
	if(found) return scale;
	else return -1;
}

bool LPFlux::isOptimal() {
	return SCIPlpiIsOptimal(_lpi);
}


bool LPFlux::isFeasible() {
	return SCIPlpiIsPrimalFeasible(_lpi);
}

double LPFlux::getObjVal() {
	double objval;
	BOOST_SCIP_CALL( SCIPlpiGetObjval(_lpi, &objval) );
	return objval;
}

void LPFlux::print() {
	int ncols = 0;
	int nrows = 0;
	int nnonz = 0;
	BOOST_SCIP_CALL( SCIPlpiGetNCols(_lpi, &ncols) );
	BOOST_SCIP_CALL( SCIPlpiGetNCols(_lpi, &nrows) );
	BOOST_SCIP_CALL( SCIPlpiGetNNonz(_lpi, &nnonz) );
	cout << "ncols: " << ncols << "  nrows: " << nrows << "  nnonz: " << nnonz << endl;

	cout << "feasible: " << isFeasible() << endl;
	cout << "optimal: " << isOptimal() << endl;

	if(isOptimal()) { // then, it is also dual feasible
		foreach(MetabolitePtr m, _model->getMetabolites()) {
			double p = getDual(m);
			if(p != 0) cout << m->getName() << ": " << p << endl;
		}
		cout << endl;
	}
	if(isFeasible()) {
		cout << "current flux val: " << getObjVal() << endl;

		foreach(ReactionPtr r, _model->getReactions()) {
			double f = getFlux(r);
			if(f > EPSILON || f < -EPSILON) cout << r->getName() << ": " << f << endl;
		}
	}
}

void LPFlux::write(const char* filename) {
	BOOST_SCIP_CALL(SCIPlpiWriteLP(_lpi, filename));
}

void LPFlux::writeState(const char* filename) {
	BOOST_SCIP_CALL(SCIPlpiWriteState(_lpi, filename));
}


void LPFlux::setExtraPotConstraints(unordered_set<PotSpaceConstraintPtr>& psc) {
	// pot constraints are variables (we are working in the dual!)
	int columns = _reactions.size() + _extraConstraints.size();
	int dstat[columns];
	//first _reactions.size() columns are variables of proper reactions and must be maintained
	// mark all variables as to be retained and then mark variables that are to be deleted
	for(int i = 0; i < columns; i++) dstat[i] = 0;
	for(unordered_map<PotSpaceConstraintPtr, int>::iterator iter = _extraConstraints.begin(); iter != _extraConstraints.end(); ) {
		if(psc.find(iter->first) == psc.end()) {
			// constraint not included anymore -> delete it
			dstat[iter->second] = 1;
			iter = _extraConstraints.erase(iter); // increases to next entry
		}
		else {
			// keep it
			iter++;
		}
	}
	BOOST_SCIP_CALL( SCIPlpiDelColset(_lpi, dstat) );
#ifndef NDEBUG
	for(unsigned int i = 0; i < _reactions.size(); i++) {
		assert(dstat[i] == (int) i); // reactions should keep indices
	}
#endif
	unsigned int end = _reactions.size() + _extraConstraints.size();
	foreach(PotSpaceConstraintPtr p, psc) {
		unordered_map<PotSpaceConstraintPtr, int>::iterator iter = _extraConstraints.find(p);
		if(iter != _extraConstraints.end()) {
			iter->second = dstat[iter->second]; // update to new index
		}
		else {
			_extraConstraints[p] = end++;
			// constraint on mu is p->_coef * mu >= 0
			// dual changes sign
			double lb = -INFINITY;
			double ub = 0;
			double obj = 0;
			vector<int> ind(p->_coef.size());
			vector<double> coef(p->_coef.size());
			foreach(Stoichiometry s, p->_coef) {
				ind.push_back(_metabolites.at(s.first));
				coef.push_back(s.second);
			}
			int beg = 0;
			BOOST_SCIP_CALL( SCIPlpiAddCols(_lpi, 1, &obj, &lb, &ub, NULL, coef.size(), &beg, ind.data(), coef.data()) );
		}
	}
	assert(end == _reactions.size() + _extraConstraints.size());
	SCIPlpiClearState(_lpi); // hmm... very ugly
	// we also have to adjust the size of the primsol vector
	_primsol.resize(end, 0);
}


vector<PotSpaceConstraintPtr> LPFlux::getActivePotConstraints() {
	typedef pair<PotSpaceConstraintPtr, int> PSCEntry;

	vector<PotSpaceConstraintPtr> out;

	foreach(PSCEntry e, _extraConstraints) {
		// potSpaceConstraints only allow reverse flux
		if(_primsol[e.second] < -EPSILON) {
			out.push_back(e.first);
		}
	}
	return out;
}


#if 0
	// these methods are for debugging only! A state of the LP can be stored and fetched later on
	void LPFlux::loadState() {
		if(SCIPlpiWasSolved(_lpi)) {
			cout << "lp was solved, storing basis" << endl;
			SCIPlpiGetBase(_lpi, cstat.data(), rstat.data());
		}
		else {
			cout << "lp was not solved, stroing no basis" << endl;
		}
	}

	void LPFlux::setOldState() {
		SCIPlpiSetBase(_lpi, cstat.data(), rstat.data());
	}

	void LPFlux::initStateInfo() {
		cstat.reserve(10000);
		rstat.reserve(10000);
	}
#endif

	SCIP_LPI* LPFlux::getLPI() {
		return _lpi;
	}

	int LPFlux::getIndex(ReactionPtr rxn) {
		return _reactions[rxn];
	}

	int LPFlux::getIndex(MetabolitePtr met) {
		return _metabolites[met];
	}

} /* namespace metaopt */
