/*
 * AntiCycleConstraint.h
 *
 *  Created on: Apr 27, 2012
 *      Author: arne
 */

#ifndef THERMOCONSTRAINTHANDLER_H_
#define THERMOCONSTRAINTHANDLER_H_

#include "objscip/objconshdlr.h"
#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/model/scip/DualPotentials.h"
#include "metaopt/model/scip/LPPotentials.h"

namespace metaopt {

class ThermoConstraintHandler: public scip::ObjConshdlr {
public:

	/**
	 * Creates a new AntiCycleConstraint handler.
	 *
	 * Generally, only important cycle, i.e. cycle that either contain objective reactions or flux forcing reactions, are checked.
	 *
	 * However, cycles that contain objective reactions should be priorized.
	 *
	 */
	ThermoConstraintHandler(ScipModelPtr model);
	virtual ~ThermoConstraintHandler();

	/**
	 * branch on the cycle of the current solution
	 */
	SCIP_RESULT branchCycle(SolutionPtr sol);

	/**
	 * enforce that the constraint contains no objective cycles (if it has, branch)
	 */
	SCIP_RESULT enforceObjectiveCycles(SolutionPtr sol);

	/**
	 * enforce infeasible sets that are either not cycles or flux-forcing cycles
	 */
	SCIP_RESULT enforceNonSimple(SolutionPtr sol);

	/**
	 * main method for enforcing constraints of this constrainthandler
	 */
	SCIP_RETCODE enforce(SCIP_CONS** cons, int nconss, SCIP_SOL* sol, SCIP_RESULT* result);

	/**
	 * check a single constraint
	 */
	SCIP_RESULT check(SolutionPtr sol);

	/**
	 * main method for checking constraints of the constrainthandler
	 */
	SCIP_RETCODE check(SCIP_CONS** cons, int nconss, SCIP_SOL* sol, SCIP_RESULT* result);

	/// Constraint enforcing method of constraint handler for LP solutions.
	virtual SCIP_RETCODE scip_enfolp(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS**        conss,              /**< array of constraints to process */
			int                nconss,             /**< number of constraints to process */
			int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
			SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
			SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
	);

	/// Constraint enforcing method of constraint handler for pseudo solutions.
	virtual SCIP_RETCODE scip_enfops(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS**        conss,              /**< array of constraints to process */
			int                nconss,             /**< number of constraints to process */
			int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
			SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
			SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
			SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
	);

	/// Feasibility check method of constraint handler for primal solutions.
	virtual SCIP_RETCODE scip_check(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS**        conss,              /**< array of constraints to process */
			int                nconss,             /**< number of constraints to process */
			SCIP_SOL*          sol,                /**< the solution to check feasibility for */
			SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
			SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
			SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
			SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
	);


	/// Variable rounding lock method of constraint handler.
	virtual SCIP_RETCODE scip_lock(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
			 *   constraint handler does not need constraints */
			int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
			int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
	);



	/// Constraint display method of constraint handler.
	virtual SCIP_RETCODE scip_print(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS*         cons,               /**< the constraint that should be displayed */
			FILE*              file                /**< the text file to store the information into */
	);

	/// Transforms constraint data into data belonging to the transformed problem.
	virtual SCIP_RETCODE scip_trans(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS*         sourcecons,         /**< source constraint to transform */
			SCIP_CONS**        targetcons          /**< pointer to store created target constraint */
	);



private:
	inline ScipModelPtr getScip();

	const boost::weak_ptr<ScipModel> _model;

	///////////////////////////////////////////////////
	// Helper variables
	///////////////////////////////////////////////////

	// for step 1
	LPFluxPtr _cycle_find; // for finding small violating cycles

	// for step 2
	LPFluxPtr _flux_simpl; // for eliminating unimportant cycles
	DualPotentialsPtr _is_find; // for finding general infeasible sets

	// for finding any violating cycle
	LPFluxPtr _cycle_test; // for testing feasibility

	// for checking feasibility
	LPPotentialsPtr _pot_test;
};

inline ScipModelPtr ThermoConstraintHandler::getScip() {
	assert(!_model.expired());
	return _model.lock();
}

inline void createThermoConstraint(ScipModelPtr model) {
	ThermoConstraintHandler* handler = new ThermoConstraintHandler(model);
	BOOST_SCIP_CALL( SCIPincludeObjConshdlr( model->getScip(), handler, TRUE ) );
}

} /* namespace metaopt */
#endif /* THERMOCONSTRAINTHANDLER_H_ */
