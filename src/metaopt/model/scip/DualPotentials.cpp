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
	SCIP_CALL( SCIPlpiCreate(&_lpi, "LPFlux", SCIP_OBJSEN_MINIMIZE) );

	// create metabolite -> index map
	// initialize beta variables
	int metabolite_index = 0;
	_num_beta_vars = 0;
	_num_metabolites = _model->getMetabolites().size();
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

void DualPotentials::optimize() {
	BOOST_SCIP_CALL( SCIPlpiSolveDual(_lpi) );

	BOOST_SCIP_CALL( SCIPlpiGetSol(_lpi, NULL, _primsol.data(), NULL, NULL, NULL) );
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


} /* namespace metaopt */
