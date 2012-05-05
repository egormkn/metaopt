/*
 * AntiCycleConstraint.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: arne
 */

#include <iostream>

#include "ThermoConstraintHandler.h"
#include "metaopt/scip/ScipError.h"

#define CONSTRAINT_NAME "ThermoConstraint"
#define CONSTRAINT_DESCRIPTION "This Constraint enforces thermodynamic feasible from the perspective of the flux variables. It is most efficient, if we do not optimize on the potentials"

// A note on priorities:
// Large priorities are executed first

// separation sucks for this kind of constraint (at least the ones I test), this is why the sepa value is low (and no sepa routine is implemented)
#define SEPAPRIORITY -100000
// we can already enforce fractional solutions (integer variables are only supplied by an addon), so here we need a high priority
// but enforcing may be expensive
// integrality enforcement is at priority 0
#define ENFOPRIORITY 2000
// checking in general involves solving an LP, hence it is quite expensive and should not be done to often
// indicator constraint has a priority of -1000000
// so, we should come definitely later
#define CHECKPRIORITY -2000000
// we don't have any good separation method
#define SEPA_FREQ -1
// we don't propagate yet
#define PROP_FREQ -1
// always work on all constraints (since we usually only have one)
#define EAGER_FREQ 1
// no presolving implemented yet
#define MAX_PRESOLVER_ROUNDS 0
// actually unimportant, since separation is not implemented,
// but if it were, we only want to run separation if all other separators didn't find anything
#define DELAY_SEPA TRUE
// also currently unimportant, since separation is not implemented.
// However propagation may be more useful in the future
#define DELAY_PROP FALSE
// also not implemented
#define DELAY_PRESOL FALSE
// also currently unimportant, since no propagation is implemented
#define PROPAGATION_TIMING SCIP_PROPTIMING_ALWAYS

using namespace boost;

namespace metaopt {

ThermoConstraintHandler::ThermoConstraintHandler(ScipModelPtr model) :
	scip::ObjConshdlr(model->getScip(),
			CONSTRAINT_NAME,
			CONSTRAINT_DESCRIPTION,
			SEPAPRIORITY,
			ENFOPRIORITY,
			CHECKPRIORITY,
			SEPA_FREQ,
			PROP_FREQ,
			EAGER_FREQ,
			MAX_PRESOLVER_ROUNDS,
			DELAY_SEPA, DELAY_PROP, DELAY_PRESOL, FALSE, PROPAGATION_TIMING),
	_smodel(model),
	_model(model->getModel())
{
	// init helper variables
	_cycle_find = LPFluxPtr( new LPFlux(model->getModel(), false));
	_cycle_find->setObjSense(false); // minimize
	_cycle_test = LPFluxPtr( new LPFlux(model->getModel(), false));
	_cycle_test->setObjSense(true); // maximize
	_flux_simpl = LPFluxPtr( new LPFlux(model->getModel(), true));
	_is_find = DualPotentialsPtr( new DualPotentials(model->getModel()));
	_pot_test = LPPotentialsPtr( new LPPotentials(model->getModel()));
}

ThermoConstraintHandler::~ThermoConstraintHandler() {
	// nothing to do
}

SCIP_RESULT ThermoConstraintHandler::enforceObjectiveCycles(SolutionPtr& sol) {
	ScipModelPtr model = getScip();
	_cycle_find->setDirectionBoundsInfty(sol, model);

	// _cycle_find is initialized to minimize
	_cycle_find->setDirectionObj(sol, model);
	// only include preference on variables that are not yet fixed to one sign
	shared_ptr<unordered_set<ReactionPtr> > fixedDirs = model->getFixedDirections();
	foreach(ReactionPtr rxn, *fixedDirs) {
		if(!rxn->isExchange())
			_cycle_find->setObj(rxn, 0);
	}

	foreach(ReactionPtr rxn, model->getModel()->getObjectiveReactions()) {
		double val = model->getFlux(sol, rxn);
		// force a tiny flow through the reaction in the current direction
		if(val > EPSILON) {
			_cycle_find->setLb(rxn, 1);
			_cycle_find->solveDual();
			if(_cycle_find->getFlux(rxn) > EPSILON) {
				return branchCycle(sol);
			}
			_cycle_find->setLb(rxn, 0); //undo the change
		}
		else if( val < -EPSILON) {
			_cycle_find->setUb(rxn, -1);
			_cycle_find->solveDual();
			if(_cycle_find->getFlux(rxn) < -EPSILON) {
				return branchCycle(sol);
			}
			_cycle_find->setUb(rxn, 0); //undo the change
		}
	}
	return SCIP_FEASIBLE;
}

SCIP_RESULT ThermoConstraintHandler::branchCycle(SolutionPtr& sol) {
	int count = 0;
	ScipModelPtr model = getScip();
	foreach(ReactionPtr rxn, model->getModel()->getReactions()) {
		double val = _cycle_find->getFlux(rxn);
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

SCIP_RESULT ThermoConstraintHandler::enforceNonSimple(SolutionPtr& sol) {
	ScipModelPtr model = getScip();
	// so lets first preprocess and get rid of the easy cycles
	// to do so, we will simply remove all the easy cycles (this is actually what CycleDeletionHeur does )
	// so, we will modify the current solution, hence, we have to copy it
	_flux_simpl->set(sol, model);
	// use cycle_test to find internal cycles, first set the objective, because it will stay the same
	// _cycle_test is initialized to maximize
	// in the loop, we will adopt the bounds
	_cycle_test->setDirectionObj(_flux_simpl);
	// iteratively search for a cycle and subtract that cycle
	bool hasFlux;
	do {
		// only allow flux through reactions that still carry flux
		_cycle_test->setDirectionBounds(_flux_simpl);
		// we don't want to find flux through flux forcing reactions
		// so set bounds there already to zero
		// (we also don't want to find flux through objective reactions, but that is already taken care of by step 1)
		// TODO: Rechanging the bounds of flux forcing reactions may be quite inefficient, if we have many flux forcing reactions
		foreach(ReactionPtr rxn, model->getModel()->getFluxForcingReactions()) {
			_cycle_test->setLb(rxn,0);
			_cycle_test->setUb(rxn,0);
		}
		_cycle_test->solve();
		if(!_cycle_test->isFeasible()) {
			// usually, this should never happen, but we may run into numerical issues
			// these numerical issues should not kill the program
			// exceptions in this method are caught and forwarded as a SCIP_ERROR,
			// so we'll simply throw an error
			BOOST_SCIP_CALL( -6 ); // Error in LP solver
		}
		hasFlux = _cycle_test->getObjVal() > EPSILON;
		if(hasFlux) {
			// now subtract computed flux from solution flux
			double scale = _flux_simpl->getSubScale(_cycle_test);
			assert(scale > 0);
			_flux_simpl->subtract(_cycle_test, scale);
		}
	}
	while(hasFlux);

	// data.flux_simpl now contains a still valid solution without easy cycles.

	// we now start looking for infeasible sets that we will branch on

	return branchIS(sol);
}

SCIP_RESULT ThermoConstraintHandler::branchIS(SolutionPtr& sol) {
	ScipModelPtr model = getScip();

	shared_ptr<unordered_set<ReactionPtr> > fixedDirs = model->getFixedDirections();
	_is_find->setDirections(_flux_simpl, fixedDirs);

	_is_find->optimize();

	if(!_is_find->isFeasible()) {
		return SCIP_FEASIBLE;
	}
	else {
		// we found an infeasible set we have to get rid by branching
		shared_ptr<unordered_set<ReactionPtr> > is = _is_find->getIS();

		int count = 0;
		foreach(ReactionPtr rxn, *is) {
			// reaction is in the basic solution, so its value is nonzero and of the same sign as in flux_simpl
			double val = _flux_simpl->getFlux(rxn);
			double lb = model->getCurrentFluxLb(rxn);
			double ub = model->getCurrentFluxUb(rxn);
			if(val > 0) {
				if(lb < EPSILON) {
					SCIP_NODE* node;
					double prio = val/model->getFlux(sol, rxn); // idea: small reductions are better than large ones
					//double prio = 1.0;
					BOOST_SCIP_CALL( SCIPcreateChild(model->getScip(), &node, prio, SCIPtransformObj(model->getScip(),SCIPgetSolOrigObj(model->getScip(), sol.get()))) ); // estimate must SCIPgetSolTransObj(scip, NULL)-lpobjval/prio)be for transformed node, sp transform estimated value for orig prob
					model->setDirection(node, rxn, false); // restrict reaction to backward direction
					count++;
				}
			}
			else {
				if(ub > -EPSILON) {
					SCIP_NODE* node;
					double prio = val/model->getFlux(sol, rxn); // idea: small reductions are better than large ones
					//double prio = 1.0;
					BOOST_SCIP_CALL( SCIPcreateChild(model->getScip(), &node, prio, SCIPtransformObj(model->getScip(),SCIPgetSolOrigObj(model->getScip(), sol.get()))) ); // estimate must SCIPgetSolTransObj(scip, NULL)-lpobjval/prio)be for transformed node, sp transform estimated value for orig prob
					model->setDirection(node, rxn, true); // restrict reaction to forward direction
					count++;
				}
			}
		}
		if(count >= 1) {
			return SCIP_BRANCHED;
		}
		else {
			return SCIP_CUTOFF;
		}
	}
}

SCIP_RESULT ThermoConstraintHandler::enforceLastResort(SolutionPtr& sol) {

	std::cout << "Warning: ThermoConstraintHandler is now checking for infeasible simple cycles (CycleDeletionHeur can deal with them!)." << std::endl;

	_flux_simpl->set(sol, getScip());
	// don't subtract easy cycles, find infeasible set directly

	SCIP_RESULT res = branchIS(sol);
	if(res == SCIP_FEASIBLE) {
		std::cout << "thermo constraint feasible" << std::endl;
	}
	return res;
}

/**
 * main method for enforcing this constraint
 */
SCIP_RETCODE ThermoConstraintHandler::enforce(SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_RESULT* result) {
	// default is feasible, unless we find something that contradicts this
	*result = SCIP_FEASIBLE;
	try {
		// 1st step: try finding objective cycles
		SolutionPtr solptr = wrap_weak(sol);
		*result = enforceObjectiveCycles(solptr);
		assert(solptr.unique());
		if(*result != SCIP_FEASIBLE) return SCIP_OKAY;

		//2nd step: try finding infeasible sets that are either not cycles or flux-forcing cycles
		// the other cycles we can get rid of easily.
		*result = enforceNonSimple(solptr);
		assert(solptr.unique());
		if(*result != SCIP_FEASIBLE) return SCIP_OKAY;

		//3rd step: Last resort, we didn't find anything difficult, so branch on something that the CycleDeletionHeur could also deal with.
		// Usually, we should not need to do this.
		// Hence, this will also print a warning
		*result = enforceLastResort(solptr);
		assert(solptr.unique());
	}
	catch( boost::exception& ex) {
		return SCIP_ERROR;
	}

	return SCIP_OKAY;
}

SCIP_RESULT ThermoConstraintHandler::check(SolutionPtr& sol) {
	_pot_test->setDirections(sol, getScip());
	bool result = false;
	if(_pot_test->testStrictFeasible(result)) {
		if(result) return SCIP_FEASIBLE;
		else return SCIP_INFEASIBLE;
	}
	else {
		BOOST_SCIP_CALL( -6 ); // LP error, was unable to check
		return SCIP_FEASIBLE; // is never executed
	}
}

SCIP_RETCODE ThermoConstraintHandler::check(SCIP_CONS** conss, int nconss, SCIP_SOL* sol, SCIP_RESULT* result) {
	// default is feasible, unless we find something that contradicts this
	*result = SCIP_FEASIBLE;
	// in contrast to the enforce method, we do not first check for objective cycles etc.
	// We can directly check feasibility using LPPotentials.
	// This is also the reason, why all those different kinds of enforcement methods are combined in one constraint handler.
	try {
		SolutionPtr solptr = wrap_weak(sol);
		*result = check(solptr);
		assert(solptr.unique());
	}
	catch( boost::exception& ex) {
		return SCIP_ERROR;
	}
	return SCIP_OKAY;
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

//----------------------------------------
SCIP_RETCODE ThermoConstraintHandler::scip_lock(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
                                              *   constraint handler does not need constraints */
      int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
      int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
      )
{
	/*
	 * The problem of this method is that it is also called during the destruction process of scip.
	 * In this case the _smodel weak_ptr is already invalidated!
	 */

	if(!_smodel.expired()) {

		ScipModelPtr smodel = _smodel.lock();

		// capture results for the time when _smodel expires
		_lockingInfos.clear();

		//I'd love to incorporate current fixing data, but I don't know if this method is supposed to return a local or a global result.
		// since this method is also called at destruction, I guess I am supposed to return globally valid information.

		//shared_ptr<unordered_set<ReactionPtr> > fixed_dirs = smodel->getFixedDirections();

		foreach(ReactionPtr rxn, smodel->getModel()->getInternalReactions()) {
			// see above
			//if(fixed_dirs->find(rxn) == fixed_dirs->end()) {
				// reaction direction is not yet fixed

				LockingInfo li;

				// if we cannot round down to zero, rounding down is safe
				li.down_safe = rxn->getUb() < EPSILON || rxn->getLb() > EPSILON;
				// if we cannot round up to zero, rounding up is safe
				li.up_safe = rxn->getLb() > -EPSILON || rxn->getUb() < -EPSILON;

				li.var = smodel->getFlux(rxn);

				if(!li.down_safe) {
					if(!li.up_safe) {
						BOOST_SCIP_CALL( SCIPaddVarLocks(scip, li.var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
					}
					else {
						BOOST_SCIP_CALL( SCIPaddVarLocks(scip, li.var, nlockspos, nlocksneg) );
					}
				}
				else if(!li.up_safe) {
					BOOST_SCIP_CALL( SCIPaddVarLocks(scip, li.var, nlocksneg, nlockspos) );
				}

				// store locking data
				_lockingInfos.push_back(li);
			//}
		}
	}
	else {
		foreach(LockingInfo li, _lockingInfos) {
			if(!li.down_safe) {
				if(!li.up_safe) {
					BOOST_SCIP_CALL( SCIPaddVarLocks(scip, li.var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
				}
				else {
					BOOST_SCIP_CALL( SCIPaddVarLocks(scip, li.var, nlockspos, nlocksneg) );
				}
			}
			else if(!li.up_safe) {
				BOOST_SCIP_CALL( SCIPaddVarLocks(scip, li.var, nlocksneg, nlockspos) );
			}
		}
	}

	return SCIP_OKAY;
}

SCIP_RETCODE ThermoConstraintHandler::scip_trans(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS*         sourcecons,         /**< source constraint to transform */
		SCIP_CONS**        targetcons          /**< pointer to store created target constraint */
		)
{
	// we don't have any constraints that need transforming
	return SCIP_OKAY;
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
