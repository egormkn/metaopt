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
 * DualPotentials.cpp
 *
 *  Created on: 02.05.2012
 *      Author: arnem
 */

#include "DualPotentials.h"

using namespace std;
using namespace boost;

namespace metaopt {

// first we need to clarify some indexing conventions

// Variables
// for each reaction, we have an \alpha variable
#define ALPHA_START _num_beta_vars
// for each metabolite, we have a \beta^+ and a \beta^- variable, indexing alternates (\beta^+_1, \beta^-_1, \beta^+_2, \beta^-_2, ...)
#define BETA_START 0
// we have a single gamma variable
#define GAMMA (_num_reactions+_num_beta_vars)
// after gamma comes the extra variables corresponding to additional potential space constraints
#define EXTRA (GAMMA + 1)
// Constraints
// for each metabolite m, we have the "\mu" constraint S \alpha + \beta^+_m - \beta^-_m = 0
#define MU_START 0
// we have one "X" constraint u \beta^+ - \ell \beta^- + \gamma \leq 0
#define X_CONSTRAINT _num_metabolites
// we have one "z" constraint \gamma + \Eins \alpha_N - \Eins \alpha_P = 1
#define Z_CONSTRAINT (_num_metabolites+1)

DualPotentials::DualPotentials(ModelPtr model) {
	_model = model;
	BOOST_SCIP_CALL( init_lp() );
}

SCIP_RETCODE DualPotentials::init_lp() {
	// we initially build the feas-test LP, because feas test should always be called before the optimization step

	_lpi = NULL;
	SCIP_CALL( SCIPlpiCreate(&_lpi, NULL, "LPFlux", SCIP_OBJSEN_MINIMIZE) );
	SCIPlpiSetIntpar(_lpi, SCIP_LPPAR_PRESOLVING, 0);

	// create metabolite -> index map
	// initialize beta variables
	int metabolite_index = 0;
	_num_beta_vars = 0;
	_num_metabolites = _model->getMetabolites().size();

	// create (empty) rows
	{
		// the new scip version seems to require that rows are created beforehand
		double lhs = 0;
		double rhs = 0;
		int beg[0];
		int ind[0];
		double var[0];
		for(int i = MU_START; i < MU_START+_num_metabolites; i++) {
			SCIP_CALL( SCIPlpiAddRows(_lpi, 1, &lhs, &rhs, NULL, 0, beg, ind, var) );
		}
		// X
		SCIP_CALL( SCIPlpiAddRows(_lpi, 1, &lhs, &rhs, NULL, 0, beg, ind, var) );
		// Z
		SCIP_CALL( SCIPlpiAddRows(_lpi, 1, &lhs, &rhs, NULL, 0, beg, ind, var) );
	}



	foreach(const MetabolitePtr m, _model->getMetabolites()) {
		_metabolites[m] = metabolite_index;
		vector<int> ind;
		// for \beta^+ var
		vector<double> coefpos;
		// for \beta^- var
		vector<double> coefneg;

		// X constraint
		ind.push_back(X_CONSTRAINT);
		// if PotUb or PotLb is INFINITY, we will not add the corresponding beta variable (but, for simplicity, we will create it)
		coefpos.push_back(m->getPotUb());
		coefneg.push_back(-m->getPotLb());

		// MU constraint
		ind.push_back(MU_START + metabolite_index);
		coefpos.push_back(1);
		coefneg.push_back(-1);

		double obj = 0;
		double lb = 0;
		double ub = INFINITY;
		int beg = 0;
		if(isinf(m->getPotUb()) == 0) {
			SCIP_CALL( SCIPlpiAddCols(_lpi, 1, &obj, &lb, &ub, NULL, ind.size(), &beg, ind.data(), coefpos.data()) );
			_num_beta_vars++;
		}
		if(isinf(m->getPotLb()) == 0) {
			SCIP_CALL( SCIPlpiAddCols(_lpi, 1, &obj, &lb, &ub, NULL, ind.size(), &beg, ind.data(), coefneg.data()) );
			_num_beta_vars++;
		}
		metabolite_index++;
	}
	assert(_num_metabolites == metabolite_index);

	// create reaction -> index map
	// initialize alpha variables
	int reaction_index = 0;
	foreach(ReactionPtr r, _model->getReactions()) {
		vector<int> ind;
		vector<double> coef;
		if(!r->isExchange()) { // potential differences of exchange reactions are meaningless
			_reactions[r] = reaction_index++;
			foreach(Stoichiometry m, r->getStoichiometries()) {
				ind.push_back(MU_START + _metabolites.at(m.first));
				coef.push_back(m.second);
			}
			// set preliminary values for bounds.
			// we require a call to setDirections before running the first optimization
			double lb = 0;
			double ub = 0;
			int beg = 0;
			double obj = 0;
			// actually we have a nice name for the row, but it wants a char* instead of a const char*. I don't think it is worth copying names ;)
			SCIP_CALL( SCIPlpiAddCols(_lpi, 1, &obj, &lb, &ub, NULL, ind.size(), &beg, ind.data(), coef.data()) );
		}
	}
	_num_reactions = reaction_index;

	// create gamma variable
	{
		vector<int> ind;
		vector<double> coef;

		ind.push_back(X_CONSTRAINT);
		coef.push_back(1);

		ind.push_back(Z_CONSTRAINT);
		coef.push_back(1);

		double lb = 0;
		double ub = INFINITY;
		double obj = 0;
		int beg = 0;

		SCIP_CALL( SCIPlpiAddCols(_lpi, 1, &obj, &lb, &ub, NULL, ind.size(), &beg, ind.data(), coef.data()) );
	}

	_primsol.resize(_num_beta_vars+_num_reactions +1, 0); // allocate sufficient memory (\beta^+, \beta^-, \alpha, \gamma)

	// set sides of constraints
	// MU
	{
		double lhs[_num_metabolites];
		double rhs[_num_metabolites];
		int ind[_num_metabolites];
		for(int i = 0; i < _num_metabolites; i++) {
			lhs[i] = 0;
			rhs[i] = 0;
			ind[i] = MU_START + i;
		}

		SCIP_CALL( SCIPlpiChgSides(_lpi, _num_metabolites, ind, lhs, rhs) );
	}

	// X
	{
		double lhs = -INFINITY;
		double rhs = 0;
		int ind = X_CONSTRAINT;

		SCIP_CALL( SCIPlpiChgSides(_lpi, 1, &ind, &lhs, &rhs) );
	}

	// Z
	{
		double lhs = 1;
		double rhs = 1;
		int ind = Z_CONSTRAINT;

		SCIP_CALL( SCIPlpiChgSides(_lpi, 1, &ind, &lhs, &rhs) );
	}

	return SCIP_OKAY;
}

DualPotentials::~DualPotentials() {
	int code = free_lp();
	// in case of error make sure that the object is destroyed anyways to keep harm as little as possible or stop the program if in debug mode.
	assert(code == SCIP_OKAY);
}

SCIP_RETCODE DualPotentials::free_lp() {
	SCIP_CALL( SCIPlpiFree(&_lpi) );
	return SCIP_OKAY;
}

void DualPotentials::setDirections(LPFluxPtr flux, shared_ptr<unordered_set<ReactionPtr> > fixed_rxns) {
	/* we have to adjust
	 * * the bounds of the \alpha variables
	 * * the \alpha coefficients of the Z constraint
	 * * the objective function
	 */

	typedef pair<ReactionPtr, int> RxnIndex;

	foreach(RxnIndex ri, _reactions) {
		double val = flux->getFlux(ri.first);
		double lb, ub, obj, coef;
		if(val > EPSILON) {
			lb = 0;
			ub = INFINITY;
			coef = 1;
			if(fixed_rxns->find(ri.first) != fixed_rxns->end()) {
				// reaction direction is already fixed, hence unimportant
				obj = 0;
			}
			else {
				obj = 1;
			}
		}
		else if(val < -EPSILON) {
			lb = -INFINITY;
			ub = 0;
			coef = -1;
			if(fixed_rxns->find(ri.first) != fixed_rxns->end()) {
				// reaction direction is already fixed, hence unimportant
				obj = 0;
			}
			else {
				obj = -1;
			}
		}
		else {
			// no flux -> potential is not constrained -> dual variable is zero
			lb = 0;
			ub = 0;
			obj = 0;
			coef = 0;
		}

		int ind = ALPHA_START + ri.second;

		BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, 1, &ind, &lb, &ub) );
		BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, 1, &ind, &obj) );
		BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, Z_CONSTRAINT, ind, coef));
	}
}

void DualPotentials::setDirections(LPFluxPtr flux, boost::unordered_map<ReactionPtr, ReactionPtr>& toFluxRxn, shared_ptr<unordered_set<ReactionPtr> > fixed_rxns) {
	/* we have to adjust
	 * * the bounds of the \alpha variables
	 * * the \alpha coefficients of the Z constraint
	 * * the objective function
	 */

	typedef pair<ReactionPtr, int> RxnIndex;

	foreach(RxnIndex ri, _reactions) {
		ReactionPtr fluxRxn;
		try {
			fluxRxn = toFluxRxn.at(ri.first);
		}
		catch(std::exception& ex) {
			BOOST_THROW_EXCEPTION(UnknownReactionError() << reaction_name(ri.first->getName()));
		}
		double val = flux->getFlux(fluxRxn);
		double lb, ub, obj, coef;
		if(val > EPSILON) {
			lb = 0;
			ub = INFINITY;
			coef = 1;
			if(fixed_rxns->find(ri.first) != fixed_rxns->end()) {
				// reaction direction is already fixed, hence unimportant
				obj = 0;
			}
			else {
				obj = 1;
			}
		}
		else if(val < -EPSILON) {
			lb = -INFINITY;
			ub = 0;
			coef = -1;
			if(fixed_rxns->find(ri.first) != fixed_rxns->end()) {
				// reaction direction is already fixed, hence unimportant
				obj = 0;
			}
			else {
				obj = -1;
			}
		}
		else {
			// no flux -> potential is not constrained -> dual variable is zero
			lb = 0;
			ub = 0;
			obj = 0;
			coef = 0;
		}

		int ind = ALPHA_START + ri.second;

		BOOST_SCIP_CALL( SCIPlpiChgBounds(_lpi, 1, &ind, &lb, &ub) );
		BOOST_SCIP_CALL( SCIPlpiChgObj(_lpi, 1, &ind, &obj) );
		BOOST_SCIP_CALL( SCIPlpiChgCoef(_lpi, Z_CONSTRAINT, ind, coef));
	}
}

void DualPotentials::optimize() {
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
}

bool DualPotentials::isFeasible() {
	return SCIPlpiIsPrimalFeasible(_lpi);
}

shared_ptr<unordered_set<ReactionPtr> > DualPotentials::getIS() {
	// TODO: does not necessarily compute a minimal infeasible set, but for now it should be ok.

	shared_ptr<unordered_set<ReactionPtr> > result(new unordered_set<ReactionPtr>());

	typedef pair<ReactionPtr, int> RxnIndex;
	foreach(RxnIndex ri, _reactions) {
		double val = _primsol.at(ALPHA_START + ri.second);
		if(val > EPSILON || val < -EPSILON) {
			result->insert(ri.first);
		}
	}
	return result;
}

double DualPotentials::getAlpha(ReactionPtr rxn) {
	unordered_map<ReactionPtr, int>::iterator iter = _reactions.find(rxn);
	if(iter == _reactions.end()) {
		return 0;
	}
	else {
		return _primsol.at(ALPHA_START + iter->second);
	}
}

void DualPotentials::setExtraPotConstraints(unordered_set<PotSpaceConstraintPtr>& psc) {
	// pot constraints are variables (we are working in the dual!)
	int columns = EXTRA + _extraConstraints.size();
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
	for(int i = 0; i < _reactions.size(); i++) {
		assert(dstat[i] == i); // reactions should keep indices
	}
#endif
	int end = EXTRA + _extraConstraints.size();
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
	assert(end == EXTRA + _extraConstraints.size());
	SCIPlpiClearState(_lpi); // hmm... very ugly
	_primsol.resize(end, 0);
}


vector<PotSpaceConstraintPtr> DualPotentials::getActivePotConstraints() {
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

} /* namespace metaopt */
