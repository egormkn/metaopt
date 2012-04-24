/*
 * LPPotentials.cpp
 *
 *  Created on: 23.04.2012
 *      Author: arnem
 */

#include <vector>
#include "LPPotentials.h"

using namespace std;

namespace metaopt {

// For the strict inequalities, we have to introduce an epsilon for now
// this constant must fit to other constraints that also use epsilons
#define REACTION_DIRECTIONS_EPSILON 1

LPPotentials::LPPotentials(ModelPtr model) {
	_model = model;
	BOOST_SCIP_CALL( init_lp() );
}

SCIP_RETCODE LPPotentials::init_lp() {
	_lpi = NULL;
	SCIP_CALL( SCIPlpiCreate(&_lpi, "LPFlux", SCIP_OBJSEN_MAXIMIZE) );

	// create metabolite -> row_index map
	int metabolite_index = 0;
	foreach(const MetabolitePtr m, _model->getMetabolites()) {
		_metabolites[m] = metabolite_index;
		metabolite_index++;
	}
	_num_metabolites = metabolite_index;

	// create flux variables,
	// create reaction -> variable_index map
	// create stoichiometric matrix
	int reaction_var = 0;
	foreach(ReactionPtr r, _model->getReactions()) {
		vector<int> ind;
		vector<double> coef;
		if(!r->isExchange()) { // potential differences of exchange reactions are meaningless
			_reactions[r] = reaction_var++;
			foreach(Stoichiometry m, r->getStoichiometries()) {
				ind.push_back(_metabolites.at(m.first));
				coef.push_back(m.second);
			}
			double lhs = -INFINITY;
			double rhs = INFINITY;
			if(r->getLb() > EPSILON) {
				rhs = -REACTION_DIRECTIONS_EPSILON;
			}
			if(r->getUb() < -EPSILON) {
				lhs = REACTION_DIRECTIONS_EPSILON;
			}
			int beg = 0;
			// actually we have a nice name for the row, but it wants a char* instead of a const char*. I don't think it is worth copying names ;)
			SCIP_CALL( SCIPlpiAddRows(_lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
		}
	}
	_num_reactions = reaction_var;
	_primsol.resize(_num_metabolites, 0); // allocate sufficient memory

	// set bounds on metabolite potentials
	double lb[_num_metabolites];
	double ub[_num_metabolites];
	double obj[_num_metabolites];
	int ind[_num_metabolites];

	typedef pair<MetabolitePtr, int> PotVar;
	foreach(PotVar v, _metabolites) {
		ind[v.second] = v.second;
		lb[v.second] = v.first->getPotLb();
		ub[v.second] = v.first->getPotUb();
		obj[v.second] = v.first->getPotObj();
	}
	SCIP_CALL( SCIPlpiChgBounds(_lpi, _num_metabolites, ind, lb, ub));
	SCIP_CALL( SCIPlpiChgObj(_lpi, _num_metabolites, ind, obj));

	return SCIP_OKAY;
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
			rhs[c.second] = -REACTION_DIRECTIONS_EPSILON;
		}
		if(val < -EPSILON) {
			lhs[c.second] = REACTION_DIRECTIONS_EPSILON;
		}
		ind[c.second] = c.second;
	}
	BOOST_SCIP_CALL( SCIPlpiChgSides(_lpi, _num_reactions, ind, lhs, rhs));
}

void LPPotentials::setDirection(ReactionPtr rxn, bool fwd) {
	double lhs = -INFINITY;
	double rhs = INFINITY;
	if(fwd) {
		rhs = -REACTION_DIRECTIONS_EPSILON;
	}
	else {
		lhs = REACTION_DIRECTIONS_EPSILON;
	}
	int ind = 0;
	BOOST_SCIP_CALL( SCIPlpiChgSides(_lpi, 1, &ind, &lhs, &rhs));
}

void LPPotentials::solvePrimal() {
	BOOST_SCIP_CALL( SCIPlpiSolvePrimal(_lpi) );
	// capture solution in _primsol
	BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
}

void LPPotentials::solveDual() {
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );
	// capture solution in _primsol
	BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
}

bool LPPotentials::isFeasible() {
	return SCIPlpiIsPrimalFeasible(_lpi);
}

double LPPotentials::getPotential(MetabolitePtr met) {
	int i = _metabolites.at(met);
	return _primsol.at(i);
}


} /* namespace metaopt */
