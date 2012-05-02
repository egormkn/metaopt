/*
 * AntiCycleConstraint.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: arne
 */

#include "ThermoConstraintHandler.h"
#include "metaopt/scip/ScipError.h"

#define CONSTRAINT_NAME "ThermoConstraint"

using namespace boost;

namespace metaopt {

ThermoConstraintHandler::ThermoConstraintHandler(ScipModelPtr model) :
	_model(model)
{
	// nothing else to do
}

ThermoConstraintHandler::~ThermoConstraintHandler() {
	// nothing to do
}

SCIP_RESULT ThermoConstraintHandler::enforceObjectiveCycles(ConstraintData &data, SolutionPtr sol) {
	ScipModelPtr model = getScip();
	data.cycle_find->setDirectionBoundsInfty(sol, model);

	data.cycle_find->setDirectionObj(sol, model);
	// only include preference on variables that are not yet fixed to one sign
	shared_ptr<unordered_set<ReactionPtr> > fixedDirs = model->getFixedDirections(model->getModel()->getInternalReactions());
	foreach(ReactionPtr rxn, fixedDirs) {
		data.cycle_find->setObj(rxn, 0);
	}

	foreach(ReactionPtr rxn, model->getModel()->getObjectiveReactions()) {
		double val = model->getFlux(sol, rxn);
		// force a tiny flow through the reaction in the current direction
		if(val > EPSILON) {
			data.cycle_find->setLb(rxn, 1);
			data.cycle_find->solveDual();
			if(data.cycle_find->getFlux(rxn) > EPSILON) {
				return branchCycle(data, sol);
			}
			data.cycle_find->setLb(rxn, 0); //undo the change
		}
		else if( val < -EPSILON) {
			data.cycle_find->setUb(rxn, -1);
			data.cycle_find->solveDual();
			if(data.cycle_find->getFlux(rxn) < -EPSILON) {
				return branchCycle(data, sol);
			}
			data.cycle_find->setUb(rxn, 0); //undo the change
		}
	}
	return SCIP_FEASIBLE;
}

SCIP_RESULT ThermoConstraintHandler::branchCycle(ConstraintData& data, SolutionPtr sol) {
	int count = 0;
	ScipModelPtr model = getScip();
	foreach(ReactionPtr rxn, model->getModel()->getReactions()) {
		double val = data.cycle_find->getFlux(rxn);
		double lb;
		double ub;
		lb = model->getCurrentFluxLb(rxn);
		ub = model->getCurrentFluxUb(rxn);
		if(val > EPSILON) {
			if(lb < EPSILON) { // ub must be positive, since positive flow is not allowed else, restriction to zero must also be allowed
				//cout << "branching ub "<<iter.getId() << endl;
				SCIP_NODE* node;
				double prio = val/model->getFlux(sol, rxn); // idea: small reductions are better than large ones
				//double prio = 1.0;
#ifdef ESTIMATE_PESSIMISTIC
				BOOST_SCIP_CALL( SCIPcreateChild(scip, &node, prio, SCIPtransformObj(scip,SCIPgetSolOrigObj(scip, NULL)-lpobjval/prio)) ); // estimate must SCIPgetSolTransObj(scip, NULL)-lpobjval/prio)be for transformed node, sp transform estimated value for orig prob
#else
				BOOST_SCIP_CALL( SCIPcreateChild(model->getScip(), &node, prio, SCIPtransformObj(model->getScip(),SCIPgetSolOrigObj(model->getScip(), sol.get()))) ); // estimate must SCIPgetSolTransObj(scip, NULL)-lpobjval/prio)be for transformed node, sp transform estimated value for orig prob
#endif
				model->setDirection(node, rxn, false); // restrict reaction to backward direction
				count ++;
			}
		}
		else if(val < -EPSILON) {
			if(ub > -EPSILON) { // lb must be negative, since negative flow is not allowed else, restriction to zero must also be allowed
				//cout << "branching lb "<<iter.getId() << endl;
				SCIP_NODE* node;
				double prio = val/model->getFlux(sol, rxn); // idea: small reductions are better than large ones
				//double prio = 1.0;
#ifdef ESTIMATE_PESSIMISTIC
				BOOST_SCIP_CALL( SCIPcreateChild(scip, &node, prio, SCIPtransformObj(scip,SCIPgetSolOrigObj(scip, NULL)-lpobjval/prio)) );
#else
				BOOST_SCIP_CALL( SCIPcreateChild(model->getScip(), &node, prio, SCIPtransformObj(model->getScip(),SCIPgetSolOrigObj(model->getScip(), sol.get()))) ); // estimate must SCIPgetSolTransObj(scip, NULL)-lpobjval/prio)be for transformed node, sp transform estimated value for orig prob
#endif
				model->setDirection(node, rxn, true); // restrict reaction to forward direction
				count ++;
			}
		}
	}

	if(count == 0) {
		//cout << "cut off" << endl;
		return SCIP_CUTOFF;
	}
	else {
		//cout << "branched " << count  << endl;
		return SCIP_BRANCHED;
	}
}

SCIP_RESULT ThermoConstraintHandler::enforceNonSimple(ConstraintData& data, SolutionPtr sol) {
	ScipModelPtr model = getScip();
	// so lets first preprocess and get rid of the easy cycles
	// to do so, we will simply remove all the easy cycles (this is actually what CycleDeletionHeur does )
	// so, we will modify the current solution, hence, we have to copy it
	data.flux_simpl->set(sol, model);
	// use cycle_test to find internal cycles
	data.cycle_test->setDirectionBounds(data.flux_simpl);
	// we don't want to find flux through flux forcing reactions
	// so set bounds there already to zero
	foreach(ReactionPtr rxn, model->getModel()->getFluxForcingReactions()) {
		data.cycle_test->setLb(rxn,0);
		data.cycle_test->setUb(rxn,0);
	}
	// iteratively search for a cycle and subtract that cycle
	bool hasFlux;
	do {
		// only allow flux through reactions that still carry flux
		data.cycle_test->setDirectionBounds(data.flux_simpl);
		// we don't want to find flux through flux forcing reactions
		// so set bounds there already to zero
		// (we also don't want to find flux through objective reactions, but that is already taken care of by step 1)
		// TODO: Rechanging the bounds of flux forcing reactions may be quite inefficient, if we have many flux forcing reactions
		foreach(ReactionPtr rxn, model->getModel()->getFluxForcingReactions()) {
			data.cycle_test->setLb(rxn,0);
			data.cycle_test->setUb(rxn,0);
		}
		data.cycle_test->solve();
		if(!data.cycle_test->isFeasible()) {
			// usually, this should never happen, but we may run into numerical issues
			// these numerical issues should not kill the program
			// exceptions in this method are caught and forwarded as a SCIP_ERROR,
			// so we'll simply throw an error
			BOOST_SCIP_CALL( -6 ); // Error in LP solver
		}
		hasFlux = data.cycle_test->getObjVal() > EPSILON;
		if(hasFlux) {
			// now subtract computed flux from solution flux
			double scale = data.flux_simpl->getSubScale(data.cycle_test);
			assert(scale > 0);
			data.flux_simpl->subtract(data.cycle_test, scale);
		}
	}
	while(hasFlux);

	// data.flux_simpl now contains a still valid solution without easy cycles.
	// we now start looking for infeasible sets

	shared_ptr<unordered_set<ReactionPtr> > fixedDirs = model->getFixedDirections(model->getModel()->getInternalReactions());
	data.is_find->setDirections(data.flux_simpl, fixedDirs);
}

/**
 * main method for enforcing this constraint
 */
SCIP_RETCODE ThermoConstraintHandler::enforce(SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_RESULT* result) {

	// 1st step: try finding objective cycles
	for(int i = 0; i < nconss; i++) {
		ConstraintData* consdata = (ConstraintData*) SCIPconsGetData(conss[i]);
		try {
			SolutionPtr solptr = wrap(sol, getScip());
			*result = enforceObjectiveCycles(*consdata, solptr);
			if(*result != SCIP_FEASIBLE) return SCIP_OKAY;
		}
		catch( boost::exception& ex) {
			return SCIP_ERROR;
		}
	}

	//2nd step: try finding infeasible sets that are either not cycles or flux-forcing cycles
	// the other cycles we can get rid of easily.
	for(int i = 0; i < nconss; i++) {
		ConstraintData* consdata = (ConstraintData*) SCIPconsGetData(conss[i]);
		try {
			SolutionPtr solptr = wrap(sol, getScip());
			*result = enforceNonSimple(*consdata, solptr);
			if(*result != SCIP_FEASIBLE) return SCIP_OKAY;
		}
		catch( boost::exception& ex) {
			return SCIP_ERROR;
		}

	}
}

SCIP_RETCODE ThermoConstraintHandler::check(SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_RESULT* result) {

}

/// Constraint enforcing method of constraint handler for LP solutions.
SCIP_RETCODE ThermoConstraintHandler::scip_enfolp(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS**        conss,              /**< array of constraints to process */
		int                nconss,             /**< number of constraints to process */
		int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
		SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
		SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
) {
	assert(scip == getScip()->getScip());
	if(solinfeasible) {
		*result = SCIP_DIDNOTRUN;
		return SCIP_OKAY;
	}
	else {
		return( enforce(conss, nconss, NULL, result) );
	}
}

/// Constraint enforcing method of constraint handler for pseudo solutions.
SCIP_RETCODE ThermoConstraintHandler::scip_enfops(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS**        conss,              /**< array of constraints to process */
		int                nconss,             /**< number of constraints to process */
		int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
		SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
		SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
		SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
) {
	assert(scip == getScip()->getScip());
	if(solinfeasible || objinfeasible) {
		*result = SCIP_DIDNOTRUN;
		return SCIP_OKAY;
	}
	else {
		return( enforce(conss, nconss, NULL, result) );
	}
}

/// Feasibility check method of constraint handler for primal solutions.
SCIP_RETCODE ThermoConstraintHandler::scip_check(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS**        conss,              /**< array of constraints to process */
		int                nconss,             /**< number of constraints to process */
		SCIP_SOL*          sol,                /**< the solution to check feasibility for */
		SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
		SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
		SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
		SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
) {
	assert(scip == getScip()->getScip());
	return check(conss, nconss, sol, result);
}


/// Constraint display method of constraint handler.
SCIP_RETCODE ThermoConstraintHandler::scip_print(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS*         cons,               /**< the constraint that should be displayed */
		FILE*              file                /**< the text file to store the information into */
) {
	fprintf(file,CONSTRAINT_NAME);
	return SCIP_OKAY;
}

} /* namespace metaopt */
