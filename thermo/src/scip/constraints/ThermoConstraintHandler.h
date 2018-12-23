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
 * AntiCycleConstraint.h
 *
 *  Created on: Apr 27, 2012
 *      Author: arne
 */

#ifndef THERMOCONSTRAINTHANDLER_H_
#define THERMOCONSTRAINTHANDLER_H_

#include "objscip/objconshdlr.h"
#include "model/impl/FullModel.h"
#include "model/scip/ReducedScipFluxModel.h"
#include "model/scip/ScipModel.h"
#include "model/scip/DualPotentials.h"
#include "model/scip/LPPotentials.h"
#include "Uncopyable.h"
#include "scip/constraints/PotBoundPropagation2.h"
#include "model/Coupling.h"
#include "model/scip/ISSupply.h"
#include "model/scip/PotSpaceConstraint.h"
#include "scip/constraints/ThermoInfeasibleSetPool.h"
#include "Properties.h"

// set to 1 to use aggregated reactions instead of the original reactions (not correctly implemented yet)
#define THERMOCONS_USE_AGGR_RXN 0

namespace metaopt {

#define THERMO_CONSTRAINT_NAME "ThermoConstraint"

class ThermoConstraintHandler: public scip::ObjConshdlr, Uncopyable {
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
	 * Add precomputed coupling information to the constraint handler.
	 * The coupling information may be modified by the constraint handler using further information form the problem.
	 * If this is not desired, please make sure to copy the coupling information before passing it to this method.
	 */
	void setCouplingHint(CouplingPtr coupling);

	/**
	 * branch on the cycle of the current solution of _cycle_find
	 */
	SCIP_RESULT branchCycle(SolutionPtr& sol);

	/**
	 * enforce that the constraint contains no objective cycles (if it has, branch)
	 */
	SCIP_RESULT enforceObjectiveCycles(SolutionPtr& sol);

	/**
	 * find an infeasible set in _flux_simpl and branch on it
	 */
	SCIP_RESULT branchIS(SolutionPtr& sol);

	/**
	 * enforce infeasible sets that are either not cycles or flux-forcing cycles
	 */
	SCIP_RESULT enforceNonSimple(SolutionPtr& sol);

	/**
	 * Look generally for infeasible sets, including cycles that are neither flux-forcing nor objective.
	 * If this constraint is infeasible, this method will not return feasible!
	 */
	SCIP_RESULT enforceLastResort(SolutionPtr& sol);

	/**
	 * main method for enforcing constraints of this constrainthandler
	 */
	SCIP_RETCODE enforce(SCIP_CONS** cons, int nconss, SCIP_SOL* sol, SCIP_RESULT* result);

	/**
	 * check a single constraint
	 */
	SCIP_RESULT check(SolutionPtr& sol);

	/**
	 * main method for checking constraints of the constrainthandler
	 */
	SCIP_RETCODE check(SCIP_CONS** cons, int nconss, SCIP_SOL* sol, SCIP_RESULT* result);

	/**
	 * uses current potential bounds to infer blocked reactions
	 */
	SCIP_RESULT propagate();

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
            SCIP_Bool completely,
            SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
	);

	virtual SCIP_RETCODE scip_prop(
			SCIP *  			scip,
			SCIP_CONSHDLR *  	conshdlr,
			SCIP_CONS **  		conss,
			int  				nconss,
			int  				nusefulconss,
			int nmarkedconss,           // New SCIP api
			SCIP_PROPTIMING proptiming, // New SCIP api
			SCIP_RESULT *  		result
	);

	/// Variable rounding lock method of constraint handler.
	virtual SCIP_RETCODE scip_lock(
			SCIP*              scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
			SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
			 *   constraint handler does not need constraints */
			SCIP_LOCKTYPE locktype,                // New argument in SCIP API
			int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
			int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
	); // Old SCIP API, see below

	/// Initializes helper variables necessary for solve. Uses presolving improvements.
	virtual SCIP_RETCODE scip_exitpre(
			SCIP* 			scip,
			SCIP_CONSHDLR*  conshdlr,
			SCIP_CONS**		conss,
			int				ncons);

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

	virtual SCIP_RETCODE scip_delete(
			SCIP* 				scip,               /**< SCIP data structure */
			SCIP_CONSHDLR*		conshdlr,           /**< the constraint handler itself */
			SCIP_CONS*			cons,				/**< constraint that will be deleted */
			SCIP_CONSDATA**		data				/**< constraint data to free */
	);


private:
	inline ScipModelPtr getScip();

	const boost::weak_ptr<ScipModel> _smodel;
	const ModelPtr _model;

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

	// this pool is used to store already found infeasible sets, so that we don't have to go looking again.
	ThermoInfeasibleSetPool _infeas_pool;


	// propagates potential bounds that can be used to detect disabled reactions
#if 0
	PotBoundPropagation2 _pbp;
#endif

#if THERMOCONS_USE_AGGR_RXN
	// reduced model after presolve
	FullModelPtr _reduced;

	// reduced scip model after presolve
	ReducedScipFluxModelPtr _reducedScip;

	boost::unordered_map<ReactionPtr, ReactionPtr> _toReducedRxn; // map translating original reactions to reduced reactions
	boost::unordered_map<ReactionPtr, ReactionPtr> _toOriginalRxn; // map translating reduced reactions to the generating original reactions
	boost::unordered_map<MetabolitePtr, MetabolitePtr> _toReducedMet; // map translating original metabolites to reduced metabolites
#endif

	CouplingPtr _coupling; // stores coupling information

	//////////////////////////////////////////////////
	// ugly hack for locking numbers
	//////////////////////////////////////////////////

	struct LockingInfo {
		SCIP_VAR* var; // weak pointer to variable with locking information (must not free, lives as long as _smodel)
		bool down_safe; // true if it is save to round variable down
		bool up_safe; // true if it is save to round variable up
	};

	std::vector<LockingInfo> _lockingInfos;

	/**
	 * Given a list of candidates to branch on, performs some simplifications,
	 * and then performs the branching.
	 * The candidate reactions are oriented in the direction that should be blocked.
	 */
	SCIP_RESULT branch(boost::unordered_set<DirectedReaction>& candidates, ISSupplyPtr iss, SolutionPtr sol);

	/**
	 * Method called by branch to setup the node with the branching decisions.
	 * Sets the direction of the potential difference, which is exactly the opposite of the flux direction.
	 */
	void setDirection(ISSupplyPtr& iss, SCIP_NODE* node, CoverReaction& c);

	/**
	 * similar as setDirection, but does not block flux.
	 * Only creates a PotSpaceConstraint that blocks flow through the combined reaction in the CoverReaction.
	 * This is done by constraining the potential space
	 */
	void reducePotSpace(ISSupplyPtr& iss, SCIP_NODE* node, CoverReaction& c);

	/**
	 * Adds and creates a new local PotSpace constraint generated by the given cover.
	 * Pot space constraints constrain the potential space.
	 * They may not be stored as explicit constraints, since the potential space may be represented implicitly.
	 */
	void addPotSpaceConstraint(PotSpaceConstraintPtr psc, SCIP_NODE* node);

};

inline ScipModelPtr ThermoConstraintHandler::getScip() {
	assert(!_smodel.expired());
	return _smodel.lock();
}

/**
 * creates default Thermo Constraint
 */
inline void createThermoConstraint(ScipModelPtr model) {
	ThermoConstraintHandler* handler = new ThermoConstraintHandler(model);
	SCIP* scip = model->getScip();
	BOOST_SCIP_CALL( SCIPincludeObjConshdlr( scip, handler, TRUE ) );
	// we now have to add the default constraint
	SCIP_CONSHDLR* hdlr = SCIPfindConshdlr( scip, THERMO_CONSTRAINT_NAME);
	SCIP_CONS* cons;
	BOOST_SCIP_CALL( SCIPcreateCons(scip, &cons, "default thermo constraint", hdlr, NULL, true, true, true, true, true, false, false, false, false, false) );
	BOOST_SCIP_CALL( SCIPaddCons(scip, cons) );
	BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/**
 * Creates Thermo constraint with hint on flux coupled reactions
 */
inline void createThermoConstraint(ScipModelPtr model, CouplingPtr c) {
	ThermoConstraintHandler* handler = new ThermoConstraintHandler(model);
	handler->setCouplingHint(c);
	SCIP* scip = model->getScip();
	BOOST_SCIP_CALL( SCIPincludeObjConshdlr( scip, handler, TRUE ) );
	// we now have to add the default constraint
	SCIP_CONSHDLR* hdlr = SCIPfindConshdlr( scip, THERMO_CONSTRAINT_NAME);
	SCIP_CONS* cons;
	BOOST_SCIP_CALL( SCIPcreateCons(scip, &cons, "default thermo constraint", hdlr, NULL, true, true, true, true, true, false, false, false, false, false) );
	BOOST_SCIP_CALL( SCIPaddCons(scip, cons) );
	BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}


} /* namespace metaopt */
#endif /* THERMOCONSTRAINTHANDLER_H_ */
