/*
 * LPPotentials.cpp
 *
 *  Created on: 23.04.2012
 *      Author: arnem
 */

#include <vector>
#include "LPPotentials.h"
#include "metaopt/Properties.h"

using namespace std;

namespace metaopt {

// For the strict inequalities, we have to introduce an epsilon for now
// this constant must fit to other constraints that also use epsilons
//#define REACTION_DIRECTIONS_EPSILON 1

#define FEASTEST_VAR 0

LPPotentials::LPPotentials(ModelPtr model) {
	_model = model;
	BOOST_SCIP_CALL( init_lp() );

	// we have the lp, we now know the size of the base
	init(_feasTest);
}

SCIP_RETCODE LPPotentials::init_lp() {
	// we initially build the feas-test LP, because feas test should always be called before the optimization step

	_lpi = NULL;
	SCIP_CALL( SCIPlpiCreate(&_lpi, "LPFlux", SCIP_OBJSEN_MAXIMIZE) );

	// create metabolite -> row_index map
	int metabolite_var = 1;
	foreach(const MetabolitePtr m, _model->getMetabolites()) {
		_metabolites[m] = metabolite_var;
		metabolite_var++;
	}
	_num_metabolites = metabolite_var-1; // we have the feastest var at index 0

	// create flux variables,
	// create reaction -> variable_index map
	// create stoichiometric matrix
	int reaction_index = 0;
	foreach(ReactionPtr r, _model->getReactions()) {
		vector<int> ind;
		vector<double> coef;
		if(!r->isExchange()) { // potential differences of exchange reactions are meaningless
			_reactions[r] = reaction_index++;
			foreach(Stoichiometry m, r->getStoichiometries()) {
				ind.push_back(_metabolites.at(m.first));
				coef.push_back(m.second);
			}
			double lhs = -INFINITY;
			double rhs = INFINITY;
			if(r->getLb() > EPSILON) {
				//rhs = -REACTION_DIRECTIONS_EPSILON;
				rhs = 0;
				// also add the variable to test strict feasibility
				ind.push_back(FEASTEST_VAR);
				coef.push_back(1);
			}
			if(r->getUb() < -EPSILON) {
				//lhs = REACTION_DIRECTIONS_EPSILON;
				lhs = 0;
				// also add the variable to test strict feasibility
				ind.push_back(FEASTEST_VAR);
				coef.push_back(-1);
			}
			int beg = 0;
			// actually we have a nice name for the row, but it wants a char* instead of a const char*. I don't think it is worth copying names ;)
			SCIP_CALL( SCIPlpiAddRows(_lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
		}
	}
	_num_reactions = reaction_index;
	_primsol.resize(_num_metabolites+1, 0); // allocate sufficient memory

	// set bounds on metabolite potentials
	double lb[_num_metabolites+1];
	double ub[_num_metabolites+1];
	int ind[_num_metabolites+1];

	// set vals for feastest var
	lb[0] = -INFINITY;
	ub[0] = 1; // else we may get unbounded feas-test solutions
	ind[0] = FEASTEST_VAR;

	typedef pair<MetabolitePtr, int> PotVar;
	foreach(PotVar v, _metabolites) {
		ind[v.second] = v.second;
		assert(v.first->getPotLb() <= v.first->getPotUb()); // am I allowed to check for equality?
		lb[v.second] = v.first->getPotLb();
		ub[v.second] = v.first->getPotUb();
		double vobj = v.first->getPotObj();
		if( vobj < -EPSILON || vobj > EPSILON) {
			_orig_obj.push_back(vobj);
			_obj_ind.push_back(v.second);
			_zero_obj.push_back(0);
		}
	}
	SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_metabolites+1, ind, lb, ub));
	double obj = 1;
	// ind already has the correct value at index 0
	SCIP_CALL( SCIPlpiChgObj(_lpi, 1, ind, &obj));

	return SCIP_OKAY;
}

void LPPotentials::init(Basis& b) {
#ifndef NDEBUG
	int ncols, nrows;
	BOOST_SCIP_CALL( SCIPlpiGetNCols(_lpi, &ncols) );
	BOOST_SCIP_CALL( SCIPlpiGetNRows(_lpi, &nrows) );
	assert(ncols == _num_metabolites+1);
	assert(nrows == _num_reactions);
#endif

	b.rstat.resize(_num_reactions,0);
	b.cstat.resize(_num_metabolites+1,0); // don't forget the extra variable for the strict feasibility test variable
	b.initialized = false;
}

SCIP_RETCODE LPPotentials::free_lp() {
	SCIP_CALL( SCIPlpiFree(&_lpi) );
	return SCIP_OKAY;
}

LPPotentials::~LPPotentials() {
	int code = free_lp();
	// in case of error make sure that the object is destroyed anyways to keep harm as little as possible or stop the program if in debug mode.
	assert(code == SCIP_OKAY);
}

void LPPotentials::setDirections(LPFluxPtr other) {
	double lhs[_num_reactions];
	double rhs[_num_reactions];
	int ind[_num_reactions];

	typedef pair<ReactionPtr, int> RxnCon;

	foreach(RxnCon c, _reactions) {
		double val = other->getFlux(c.first);
		lhs[c.second] = -INFINITY;
		rhs[c.second] = INFINITY;
		if(val > EPSILON) {
			//rhs[c.second] = -REACTION_DIRECTIONS_EPSILON;
			rhs[c.second] = 0;
			BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, c.second, 0, 1) );
		}
		else if(val < -EPSILON) {
			//lhs[c.second] = REACTION_DIRECTIONS_EPSILON;
			lhs[c.second] = 0;
			BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, c.second, 0, -1) );
		}
		else {
			BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, c.second, 0, 0) );
		}
		ind[c.second] = c.second;
	}
	BOOST_SCIP_CALL( SCIPlpiChgSides(_lpi, _num_reactions, ind, lhs, rhs));
}

void LPPotentials::setDirections(SolutionPtr sol, ScipModelPtr smodel) {
	double lhs[_num_reactions];
	double rhs[_num_reactions];
	int ind[_num_reactions];

	typedef pair<ReactionPtr, int> RxnCon;

	foreach(RxnCon c, _reactions) {
		double val = smodel->getFlux(sol, c.first);
		lhs[c.second] = -INFINITY;
		rhs[c.second] = INFINITY;
		if(val > EPSILON) {
			//rhs[c.second] = -REACTION_DIRECTIONS_EPSILON;
			rhs[c.second] = 0;
			BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, c.second, 0, 1) );
		}
		else if(val < -EPSILON) {
			//lhs[c.second] = REACTION_DIRECTIONS_EPSILON;
			lhs[c.second] = 0;
			BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, c.second, 0, -1) );
		}
		else {
			BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, c.second, 0, 0) );
		}
		ind[c.second] = c.second;
	}
	BOOST_SCIP_CALL( SCIPlpiChgSides(_lpi, _num_reactions, ind, lhs, rhs));
}

void LPPotentials::setDirection(ReactionPtr rxn, bool fwd) {
	double lhs = -INFINITY;
	double rhs = INFINITY;
	if(fwd) {
		//rhs = -REACTION_DIRECTIONS_EPSILON;
		rhs = 0;
	}
	else {
		//lhs = REACTION_DIRECTIONS_EPSILON;
		lhs = 0;
	}
	int ind = _reactions.at(rxn);
	BOOST_SCIP_CALL( SCIPlpiChgSides(_lpi, 1, &ind, &lhs, &rhs));
}

bool LPPotentials::optimize() {
	// set original objective
	if(!_obj_ind.empty()) { // (if we have no objective potentials, we must not do this, else we get NULL-pointer issues)
		BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, _obj_ind.size(), _obj_ind.data(), _orig_obj.data()));
	}
	// we only changed objective, so use primal simplex
	BOOST_SCIP_CALL( SCIPlpiSolvePrimal(_lpi) );

	if(! SCIPlpiIsOptimal(_lpi) ) {
		return false; // we somehow failed to solve the LP. Thus, we cannot determine if it is strictly feasible
	}

	// capture solution in _primsol
	BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
	return true;
}

bool LPPotentials::testStrictFeasible(bool& result) {
	// actually solve and test feasibility

	// set feastest objective
	if(!_obj_ind.empty()) { // (if we have no objective potentials, we must not do this, else we get NULL-pointer issues)
		BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, _obj_ind.size(), _obj_ind.data(), _zero_obj.data()));
	}
	int ind = FEASTEST_VAR;
	double obj = 1;
	BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, 1, &ind, &obj));

	// ignore basis information stuff, because we will deactivate some of the constraints (by setting bounds to inf) from time to time
	// set base
	//if(_feasTest.initialized) {
		//BOOST_SCIP_CALL( SCIPlpiSetBase(_lpi, _feasTest.cstat.data(), _feasTest.rstat.data()) );
	//}

	// solve
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );

#ifdef LPSOLVER_SOPLEX
	// the status code -4 has only this meaning for soplex
	if(SCIPlpiGetInternalStatus(_lpi) == -4) {
		// basis is singular, we should resolve from scratch
		BOOST_SCIP_CALL( SCIPlpiClearState(_lpi) );
		// resolve
		BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
	}
#endif

	// by the structure of the problem, the problem should always be primal feasible, however numerical issues may prevent the solver from seeing this
	// in this case, we will not want to throw a runtime error
	// hence we do not check by assert
	if(! SCIPlpiIsPrimalFeasible(_lpi) && ! SCIPlpiIsOptimal(_lpi) ) {
		//SCIPlpiWriteLP(_lpi, "debug.lp");
		return false; // we somehow failed to solve the LP. Thus, we cannot determine if it is strictly feasible
	}

	// ignore basis information stuff, only causes problems
	// store base
	//BOOST_SCIP_CALL( SCIPlpiGetBase(_lpi, _feasTest.cstat.data(), _feasTest.rstat.data()) );
	//_feasTest.initialized = true;

	BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );

	double val; // objective value
	BOOST_SCIP_CALL( SCIPlpiGetObjval(_lpi, &val) );

	result = val >= EPSILON; // test against epsilon (because of rounding issues)
	return true;
}

double LPPotentials::getPotential(MetabolitePtr met) {
	int i = _metabolites.at(met);
	return _primsol.at(i);
}

bool LPPotentials::isCurrentSolutionFeasible() {
	return SCIPlpiIsPrimalFeasible(_lpi);
}

} /* namespace metaopt */
