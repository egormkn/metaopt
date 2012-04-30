/*
 * AntiCycleConstraint.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: arne
 */

#include "AntiCycleConstraint.h"
#include "metaopt/scip/ScipError.h"

#define CONSTRAINT_NAME "AntiCycleConstraint"

namespace metaopt {

AntiCycleConstraint::AntiCycleConstraint(ScipModelPtr model, bool objective) :
	_model(model), _objective(objective)
{
	// nothing else to do
}

AntiCycleConstraint::~AntiCycleConstraint() {
	// nothing to do
}

SCIP_RESULT AntiCycleConstraint::enforce(ConstraintData &data, SolutionPtr sol) {
	ScipModelPtr model = getScip();
	data.flux->setDirectionBounds(sol, model);
	data.flux->setDirectionObj(sol, model);
	foreach(ReactionPtr rxn, getImportantReactions()) {
		double val = model->getFlux(sol, rxn);
		// force a tiny flow through the reaction in the current direction
		if(val > EPSILON) {
			data.flux->setLb(rxn, 2*EPSILON);
			data.flux->solveDual();
			if(data.flux->getFlux(rxn) > EPSILON) {
				analyzeCycle();
			}
			data.flux->setLb(rxn, 0); //undo the change
		}
		else if( val < -EPSILON) {
			data.flux->setUb(rxn, -2*EPSILON);
			data.flux->solveDual();
			if(data.flux->getFlux(rxn) < -EPSILON) {
				analyzeCycle();
			}
			data.flux->setUb(rxn, 0); //undo the change
		}
	}

}

/**
 * main method for enforcing this constraint
 */
SCIP_RETCODE AntiCycleConstraint::enforce(SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_RESULT* result) {
	for(int i = 0; i < nconss; i++) {
		ConstraintData* consdata = (ConstraintData*) SCIPconsGetData(conss[i]);
		try {
			SolutionPtr solptr = wrap(sol, getScip());
			*result = enforce(*consdata, solptr);
			return SCIP_OKAY;
		}
		catch( boost::exception& ex) {
			return SCIP_ERROR;
		}
	}
}

SCIP_RETCODE AntiCycleConstraint::check(SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_RESULT* result) {

}

/// Constraint enforcing method of constraint handler for LP solutions.
SCIP_RETCODE AntiCycleConstraint::scip_enfolp(
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
SCIP_RETCODE AntiCycleConstraint::scip_enfops(
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
SCIP_RETCODE AntiCycleConstraint::scip_check(
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
SCIP_RETCODE AntiCycleConstraint::scip_print(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS*         cons,               /**< the constraint that should be displayed */
		FILE*              file                /**< the text file to store the information into */
) {
	assert(scip == getScip()->getScip());
	if(objective) {
		fprintf(file,CONSTRAINT_NAME " (only cycles with objective reactions)");
	}
	else {
		fprintf(file,CONSTRAINT_NAME " (only cycles with flux forcing reactions)");
	}
	return SCIP_OKAY;
}

} /* namespace metaopt */
