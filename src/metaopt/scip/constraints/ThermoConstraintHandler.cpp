/*
 * AntiCycleConstraint.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: arne
 */
#include "metaopt/Properties.h"

#include <iostream>
#include <vector>
#include "scip/scipdefplugins.h"
#include "scip/scip.h"

#include "metaopt/scip/ScipError.h"
#include "ThermoConstraintHandler.h"
#include "metaopt/model/scip/ReducedScipFluxModel.h"
#include "metaopt/model/impl/FullModel.h"
#include "boost/shared_ptr.hpp"
#include "metaopt/model/Metabolite.h"
#include "metaopt/model/scip/PotSpaceConstraint.h"

#define CONSTRAINT_NAME THERMO_CONSTRAINT_NAME
#define CONSTRAINT_DESCRIPTION "This Constraint enforces thermodynamic feasible from the perspective of the flux variables. It is most efficient, if we do not optimize on the potentials"

//#define FINDBUG
#define LOGBRANCHING

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
// the propagation method is quite a beast... (solving several LPs)
//#define PROP_FREQ 2
#define PROP_FREQ -1
// always work on all constraints (since we usually only have one)
#define EAGER_FREQ 1
// no presolving implemented yet
#define MAX_PRESOLVER_ROUNDS 0
// actually unimportant, since separation is not implemented,
// but if it were, we only want to run separation if all other separators didn't find anything
#define DELAY_SEPA TRUE
// propagation method is very expensive, so only execute it, if no other propagators found anything
#define DELAY_PROP TRUE
// also not implemented
#define DELAY_PRESOL FALSE
// I don't know how the LP-loop can help the propagator, so only execute it before the LP-loop
#define PROPAGATION_TIMING SCIP_PROPTIMING_BEFORELP

using namespace boost;
using namespace std;

struct SCIP_ConsData {
	metaopt::PotSpaceConstraintPtr val;
	SCIP_ConsData(metaopt::PotSpaceConstraintPtr v) : val(v) {}
};

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
			DELAY_SEPA, DELAY_PROP, DELAY_PRESOL, TRUE, PROPAGATION_TIMING), // set it needs cons (although it actually doesn't), since we want to add constraints on the virtual potential space to this handler
	_smodel(model),
	_model(model->getModel()) //,
#if 0
	_pbp(model->getModel())
#endif
{
	// initialization of helper variables is not done after presolving
	// currently we do not account for improvements of the presolver that way.

	// init helper variables
	_cycle_find = LPFluxPtr( new LPFlux(_model, false));
	_cycle_find->setObjSense(false); // minimize
	_cycle_test = LPFluxPtr( new LPFlux(_model, false));
	_cycle_test->setObjSense(true); // maximize
	_flux_simpl = LPFluxPtr( new LPFlux(_model, true));
	_is_find = DualPotentialsPtr( new DualPotentials(_model));
	_pot_test = LPPotentialsPtr( new LPPotentials(_model)); // this cannot be initialized after presolving, because check may already be run earlier
}

ThermoConstraintHandler::~ThermoConstraintHandler() {
	// nothing to do
}

void ThermoConstraintHandler::setCouplingHint(CouplingPtr coupling) {
	_coupling = coupling;
}

SCIP_RESULT ThermoConstraintHandler::enforceObjectiveCycles(SolutionPtr& sol) {
	ScipModelPtr model = getScip();
#if THERMOCONS_USE_AGGR_RXN
	_cycle_find->setDirectionBoundsInfty(sol, _reducedScip);
#else
	_cycle_find->setDirectionBoundsInfty(sol, model);
#endif
	// _cycle_find is initialized to minimize
#if THERMOCONS_USE_AGGR_RXN
	_cycle_find->setDirectionObj(sol, _reducedScip);
#else
	_cycle_find->setDirectionObj(sol, model);
#endif
	// only include preference on variables that are not yet fixed to one sign
	shared_ptr<unordered_set<ReactionPtr> > fixedDirs = model->getFixedDirections();
	foreach(ReactionPtr rxn, *fixedDirs) {
		if(!rxn->isExchange())
#if THERMOCONS_USE_AGGR_RXN
			_cycle_find->setObj(_toReducedRxn[rxn], 0);
#else
			_cycle_find->setObj(rxn, 0);
#endif
	}

#if THERMOCONS_USE_AGGR_RXN
	foreach(ReactionPtr rxn, _reduced->getObjectiveReactions()) {
#else
	foreach(ReactionPtr rxn, _model->getObjectiveReactions()) {
#endif
		if(!rxn->isExchange()) {
#if THERMOCONS_USE_AGGR_RXN
			double val = _reducedScip->getFlux(sol, rxn);
#else
			double val = model->getFlux(sol, rxn);
#endif
			// force a tiny flow through the reaction in the current direction
			if(val > EPSILON) {
				_cycle_find->setLb(rxn, 1);
				_cycle_find->solveDual();
				if(_cycle_find->isFeasible()) {
					return branchCycle(sol);
				}
				_cycle_find->setLb(rxn, 0); //undo the change
			}
			else if( val < -EPSILON) {
				_cycle_find->setUb(rxn, -1);
				_cycle_find->solveDual();
				if(_cycle_find->isFeasible()) {
					return branchCycle(sol);
				}
				_cycle_find->setUb(rxn, 0); //undo the change
			}
		}
	}
	return SCIP_FEASIBLE;
}

SCIP_RESULT ThermoConstraintHandler::branchCycle(SolutionPtr& sol) {
	ScipModelPtr model = getScip();

	//TODO: its a waste computing this twice
	shared_ptr<unordered_set<ReactionPtr> > fixedDirs = model->getFixedDirections();

	unordered_set<DirectedReaction> branchingCandidates;

#if THERMOCONS_USE_AGGR_RXN
	foreach(ReactionPtr rxn, _reduced->getReactions()) {
		if(fixedDirs->find(_toOriginalRxn[rxn]) == fixedDirs->end()) { // not fixed
			double val = _cycle_find->getFlux(rxn);
			double lb;
			double ub;
			lb = _reducedScip->getCurrentFluxLb(rxn);
			ub = _reducedScip->getCurrentFluxUb(rxn);
#else
	foreach(ReactionPtr rxn, _model->getReactions()) {
		if(fixedDirs->find(rxn) == fixedDirs->end()) { // not fixed
			double val = _cycle_find->getFlux(rxn);
			double lb;
			double ub;
			lb = model->getCurrentFluxLb(rxn);
			ub = model->getCurrentFluxUb(rxn);
#endif
			if(val > EPSILON) {
				if(lb < EPSILON) { // ub must be positive, since positive flow is not allowed else, restriction to zero must also be allowed
					//cout << "branching ub "<<iter.getId() << endl;
					branchingCandidates.insert(DirectedReaction(rxn, true));
				}
			}
			else if(val < -EPSILON) {
				if(ub > -EPSILON) { // lb must be negative, since negative flow is not allowed else, restriction to zero must also be allowed
					//cout << "branching lb "<<iter.getId() << endl;
					branchingCandidates.insert(DirectedReaction(rxn, false));
				}
			}
		}
	}

	return branch(branchingCandidates, _cycle_find, sol);
}

SCIP_RESULT ThermoConstraintHandler::enforceNonSimple(SolutionPtr& sol) {
	// so lets first preprocess and get rid of the easy cycles
	// to do so, we will simply remove all the easy cycles (this is actually what CycleDeletionHeur does )
	// so, we will modify the current solution, hence, we have to copy it
#if THERMOCONS_USE_AGGR_RXN
	_flux_simpl->set(sol, _reducedScip);
#else
	ScipModelPtr scip = getScip();
	_flux_simpl->set(sol, scip);
#endif
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
#if THERMOCONS_USE_AGGR_RXN
		foreach(ReactionPtr rxn, _reduced->getFluxForcingReactions()) {
#else
			foreach(ReactionPtr rxn, _model->getFluxForcingReactions()) {
#endif
			//foreach(ReactionPtr rxn, _reduced->getProblematicReactions()) {
			if(!rxn->isExchange()) {
				_cycle_test->setLb(rxn,0);
				_cycle_test->setUb(rxn,0);
			}
		}
#if THERMOCONS_USE_AGGR_RXN
		foreach(ReactionPtr rxn, _reduced->getProblematicReactions()) {
#else
			foreach(ReactionPtr rxn, _model->getProblematicReactions()) {
#endif
			if(!rxn->isExchange()) {
				_cycle_test->setLb(rxn,0);
				_cycle_test->setUb(rxn,0);
			}
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
	// _flux_simpl is in the reduced space, but _is_find is not.
	// _is_find is not in the reduced space, because else we might loose metabolites and thus may loose infeasible sets.
#if THERMOCONS_USE_AGGR_RXN
	_is_find->setDirections(_flux_simpl, _toReducedRxn, fixedDirs);
#else
	_is_find->setDirections(_flux_simpl, fixedDirs);
#endif

	cout << "general is" << endl;

	_is_find->optimize();

	if(!_is_find->isFeasible()) {
		return SCIP_FEASIBLE;
	}
	else {
		// we found an infeasible set we have to get rid by branching
		shared_ptr<unordered_set<ReactionPtr> > is = _is_find->getIS();

#if 0
		LPFlux debugFlux(_model, true);
		debugFlux.setObjSense(true);
		bool good = false;
#endif

		unordered_set<DirectedReaction> branchingCandidates;

		foreach(ReactionPtr rxn, *is) {
			assert(!rxn->isExchange());
			if(fixedDirs->find(rxn) == fixedDirs->end()) { // not fixed
				// reaction is in the basic solution, so its value is nonzero and of the same sign as in flux_simpl
#if THERMOCONS_USE_AGGR_RXN
				ReactionPtr aggrrxn = _toReducedRxn[rxn];
				double val = _flux_simpl->getFlux(aggrrxn);
				double lb = _reducedScip->getCurrentFluxLb(aggrrxn);
				double ub = _reducedScip->getCurrentFluxUb(aggrrxn);
#else
				double val = _flux_simpl->getFlux(rxn);
				double lb = model->getCurrentFluxLb(rxn);
				double ub = model->getCurrentFluxUb(rxn);
#endif
				//std::cout << rxn->toString() << "; " << lb << " (" << rxn->getLb() << ") <= " << val << " <= " << ub << " (" << rxn->getUb() << ")" << std::endl;
				if(val > 0) {
					if(lb < EPSILON) {
						branchingCandidates.insert(DirectedReaction(rxn, true));
#if 0
						debugFlux.setBounds(model);
						debugFlux.setUb(rxn, 0);
						debugFlux.solve();
						double debugVal = debugFlux.getObjVal();
						std::cout << debugVal << std::endl;
						good = good || debugVal > EPSILON;
#endif

					}
				}
				else {
					if(ub > -EPSILON) {
						branchingCandidates.insert(DirectedReaction(rxn, false));
#if 0
						debugFlux.setBounds(model);
						debugFlux.setLb(rxn, 0);
						debugFlux.solve();
						double debugVal = debugFlux.getObjVal();
						std::cout << debugVal << std::endl;
						good = good || debugVal > EPSILON;
#endif
					}
				}
			}
		}

#if 0
		if(!good) {
			std::cout << "not good!, current sol: " << SCIPgetSolOrigObj(model->getScip(), sol.get()) << std::endl;
		}
#endif

		return branch(branchingCandidates, _is_find, sol);
	}
}

void ThermoConstraintHandler::reduceFluxSpace(ISSupplyPtr& iss, SCIP_NODE*& node, CoverReaction& c) {
	ScipModelPtr model = getScip();
	PotSpaceConstraintPtr psc(new PotSpaceConstraint());
	psc->_coef = unordered_map<MetabolitePtr, double>(c.reaction._rxn->getStoichiometries());
	double a = iss->getAlpha(c.reaction._rxn);
	foreach(DirectedReaction d, *(c.covered)) {
		double da = iss->getAlpha(d._rxn);
		foreach(Stoichiometry s, d._rxn->getStoichiometries()) {
			double sval = 0;
			unordered_map<MetabolitePtr, double>::iterator iter = psc->_coef.find(s.first);
			if(iter != psc->_coef.end()) sval = iter->second;
			sval += (da / a) * s.second;
			psc->_coef[s.first] = sval;
		}
	}
	// the pot space constraint will always be formulated as _coef * mu >= 0
	if(!c.reaction._fwd) {
		// if we block backward flux, the potential difference should be negative, hence we have to switch the signs
		foreach(Stoichiometry& s, psc->_coef) {
			s.second *= -1;
		}
	}
	// reduce the indirect potential space
	addPotSpaceConstraint(psc, node);

	// if we have potential variables, we can also reduce the potential space directly
	if(model->hasPotentials()) {
		// if we block forward flux, this means that the potential difference must be nonnegative
		double lhs = 0;
		double rhs = INFINITY;

		// transform data into usable format
		vector<SCIP_VAR*> vars;
		vector<double> coef;
		foreach(Stoichiometry s, psc->_coef) {
			assert(model->hasPotentialVar(s.first)); // if we have some potential variables, we assume, that we have all potential variables
			vars.push_back(model->getPotential(s.first));
			coef.push_back(s.second);
		}

		SCIP_CONS* cons;
		BOOST_SCIP_CALL( SCIPcreateConsLinear(model->getScip(), &cons, "branchCover", vars.size(), vars.data(), coef.data(), lhs, rhs, false, true, true, true, true, true, false, false, true, false) );
		BOOST_SCIP_CALL( SCIPaddConsNode(model->getScip(), node, cons, NULL) );
		BOOST_SCIP_CALL( SCIPreleaseCons(model->getScip(), &cons) );
	}
}

void ThermoConstraintHandler::setDirection(ISSupplyPtr& iss, SCIP_NODE*& node, CoverReaction& c) {
	ScipModelPtr model = getScip();
	if(c.covered->empty()) {
		model->setDirection(node, c.reaction._rxn, !c.reaction._fwd);
		// actually, we should test, if the model has some method to store fixed directions
		// if it does, we do not need to add the extra PotSpaceConstraint
		// since we don't check, we'll always add the constraint
		PotSpaceConstraintPtr psc(new PotSpaceConstraint());
		if(c.reaction._fwd) {
			foreach(Stoichiometry s, c.reaction._rxn->getStoichiometries()) {
				psc->_coef[s.first] = s.second;
			}
		}
		else {
			foreach(Stoichiometry s, c.reaction._rxn->getStoichiometries()) {
				psc->_coef[s.first] = -s.second;
			}
		}
		// reduce the indirect potential space
		addPotSpaceConstraint(psc, node);
	}
	else {
		// we definitely can block the flux through the covering reaction
		model->setBlockedFlux(node, c.reaction._rxn, c.reaction._fwd);

		// on the branching step, we will not set the direction of a reaction,
		// but only disable flux through the covering reaction.
		// However, we can also constrain the flux space (and should do so),
		// by enforcing that the potential difference of the sum of the covering and the covered reactions to have a specific sign (depending on fwd).
		reduceFluxSpace(iss, node, c);
	}
}

SCIP_RESULT ThermoConstraintHandler::branch(unordered_set<DirectedReaction>& branchingCandidates, ISSupplyPtr iss, SolutionPtr sol) {
	ScipModelPtr model = getScip();

	// compute a cover
	shared_ptr<vector<CoverReaction> > cover = _coupling->computeCover(branchingCandidates);

#ifdef LOGBRANCHING
	cout << "Input Branching Set: ";
	foreach(DirectedReaction d, branchingCandidates) {
		cout << d._rxn->getName() << (d._fwd?"_fwd ":"_bwd ");
	}
	cout << endl;
	cout << "Computed Cover: ";
	foreach(CoverReaction c, *cover) {
		cout << c.reaction._rxn->getName() << " ( ";
		foreach(DirectedReaction d, *c.covered) {
			cout << d._rxn->getName() << " ";
		}
		cout << ") ";
	}
	cout << endl;
	SCIP_NODE* current = SCIPgetCurrentNode(model->getScip());
	cout << "children of: " << SCIPnodeGetNumber(current) << " are ";
#endif

	if(cover->size() == 0) {

#ifdef LOGBRANCHING
		cout << " cutoff" << endl;
#endif
		return SCIP_CUTOFF;
	}
	else if(cover->size() == 1) {
		CoverReaction c = *(cover->begin());
		SCIP_NODE* node = SCIPgetCurrentNode(model->getScip());
		setDirection(iss, node, c); // restrict reaction to backward direction

#ifdef LOGBRANCHING
		cout << " reduced dom" << endl;
#endif
		return SCIP_REDUCEDDOM;
	}
	else {
		// we have to branch

		// to ensure that we create a partition, we don't only block directions, but branch on the following:
		// Dmu_r >= 0 or
		// Dmu_s >= 0 and Dmu_r <= 0 or
		// Dmu_r >= 0 and Dmu_s <= 0 and Dmu_r <= 0 and so forth
		// so we have to remember the old reactions we already created branching decisions for
		// but we have to store them in reverse direction
		vector<CoverReaction> additional;

		foreach(CoverReaction c, *cover) {
			ReactionPtr& rxn = c.reaction._rxn;
			//double val = iss->getAlpha(rxn);
			SCIP_NODE* node;
#if THERMOCONS_USE_AGGR_RXN
			double prio = val/_reducedScip->getFlux(sol, aggrrxn); // idea: small reductions are better than large ones
#else
			//double prio = val/model->getFlux(sol, rxn); // idea: small reductions are better than large ones
			double prio = model->getFlux(sol, rxn); // idea: small reductions are better than large ones
#endif
			//double prio = 1.0;
			double estimate = SCIPgetLocalTransEstimate(model->getScip()); // it may not give any improvement
			//double estimate = SCIPtransformObj(model->getScip(),SCIPgetSolOrigObj(model->getScip(), sol.get())); // estimate must SCIPgetSolTransObj(scip, NULL)-lpobjval/prio)be for transformed node, sp transform estimated value for orig prob
			BOOST_SCIP_CALL( SCIPcreateChild(model->getScip(), &node, prio, estimate) );
#ifdef LOGBRANCHING
			cout << SCIPnodeGetNumber(node) << " ";
#endif
			setDirection(iss, node, c); // restrict reaction to backward direction
			// also enforce potential directions on additional reactions
#ifndef FINDBUG
			foreach(CoverReaction ca, additional) {
				if(ca.covered->empty()) {
					// also fix flux direction
					// if the reaction is already fixed to one direction, don't fix the flux direction again
					if(ca.reaction._fwd) {
						if(model->getCurrentFluxUb(ca.reaction._rxn) > EPSILON) {
							setDirection(iss, node, ca);
						}
						else {
							// only fix pot difference sign
							reduceFluxSpace(iss, node, ca);
						}
					}
					else {
						if(model->getCurrentFluxLb(ca.reaction._rxn) < -EPSILON) {
							setDirection(iss, node, ca);
						}
						else {
							// only fix pot difference sign
							reduceFluxSpace(iss, node, ca);
						}
					}
				}
				else {
					// only fix pot difference sign
					reduceFluxSpace(iss, node, ca);
				}
			}

			// add new additional reaction
			CoverReaction a = c;
			a.reaction._fwd = !a.reaction._fwd;
			foreach(DirectedReaction& d, *(a.covered)) {
				d._fwd = !d._fwd;
			}
			additional.push_back(a);
#endif
		}

#ifdef LOGBRANCHING
		cout << endl;
#endif
		return SCIP_BRANCHED;
	}
}

SCIP_RESULT ThermoConstraintHandler::enforceLastResort(SolutionPtr& sol) {

	std::cout << "Warning: ThermoConstraintHandler is now checking for infeasible simple cycles (CycleDeletionHeur can deal with them!)." << std::endl;

#if THERMOCONS_USE_AGGR_RXN
	_flux_simpl->set(sol, _reducedScip);
#else
	ScipModelPtr scip = getScip();
	_flux_simpl->set(sol, scip);
#endif
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
	if(!getScip()->hasCurrentFlux()) {
		* result = SCIP_DIDNOTRUN;
		return SCIP_OKAY;
	}
	// default is feasible, unless we find something that contradicts this
	*result = SCIP_FEASIBLE;
	try {

		//configure additional cons to helper variables
		unordered_set<PotSpaceConstraintPtr> extra;
		for(int i = 0; i < nconss; i++) {
			SCIP_ConsData* data = SCIPconsGetData(conss[i]);
			if(data != NULL && data->val != NULL) { // there is exactly one constraint with ConsData == NULL, which is the initial constraint
				extra.insert(data->val);
			}
		}
		_cycle_find->setExtraPotConstraints(extra); // used to find objecitve cycles
		_is_find->setExtraPotConstraints(extra); // used to find general infeasible sets
		// the other LPFluxS are not used as ISSuply


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
	catch( std::exception& ex) {
		std::cerr << boost::diagnostic_information(ex) << std::endl;
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

SCIP_RESULT ThermoConstraintHandler::propagate() {
#if 0
	ScipModelPtr scip = getScip();
	_pbp.update(scip);
	shared_ptr<std::vector<std::pair<ReactionPtr,bool> > > blocked = _pbp.getBlockedReactions();

	bool propagated = false;
	for(unsigned int i = 0; i < blocked->size(); i++) {
		ReactionPtr r = blocked->at(i).first;
		bool fwd = blocked->at(i).second;
		if(fwd && scip->getCurrentFluxUb(r) > EPSILON) {
			if(scip->getCurrentFluxLb(r) > EPSILON) {
				return SCIP_CUTOFF; // blocking flux forcing reactions is not allowed!
			}
			scip->setBlockedFlux(NULL, r, true); // change the current node
			propagated = true;
		}
		if(!fwd && scip->getCurrentFluxLb(r) < -EPSILON) {
			if(scip->getCurrentFluxUb(r) < -EPSILON) {
				return SCIP_CUTOFF; // blocking flux forcing reactions is not allowed!
			}
			scip->setBlockedFlux(NULL, r, false); // change the current node
			propagated = false;
		}
	}
	if(propagated) {
		return SCIP_REDUCEDDOM;
	}
	else {
		return SCIP_DIDNOTFIND;
	}
#else
	return SCIP_DIDNOTRUN;
#endif
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


SCIP_RETCODE ThermoConstraintHandler::scip_prop(
		SCIP *  			scip,
		SCIP_CONSHDLR *  	conshdlr,
		SCIP_CONS **  		conss,
		int  				nconss,
		int  				nusefulconss,
		SCIP_RESULT *  		result
) {
	assert(scip == getScip()->getScip());
	try {
		*result = propagate();
		return SCIP_OKAY;
	}
	catch( boost::exception& ex) {
		return SCIP_ERROR;
	}
	return SCIP_ERROR;
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

				// if we cannot round down away zero, rounding down is safe
				li.down_safe = rxn->getUb() < -EPSILON || rxn->getLb() > -EPSILON;
				// if we cannot round up away zero, rounding up is safe
				li.up_safe = rxn->getLb() > EPSILON || rxn->getUb() < EPSILON;

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

SCIP_RETCODE ThermoConstraintHandler::scip_exitpre(
		SCIP* 			scipscip,
		SCIP_CONSHDLR*  conshdlr,
		SCIP_CONS**		conss,
		int				ncons)
{
	// here, we create all the helper variables,
	// however we also want to account for model improvements the presolver derived
	// these may, for example, reduce the sizes of the circuits

	ScipModelPtr scip = getScip();
	assert(scipscip == scip->getScip());

#if THERMOCONS_USE_AGGR_RXN
	_reducedScip = ReducedScipFluxModelPtr( new ReducedScipFluxModel(scip) );
	_reduced = FullModelPtr(new FullModel());
	// first create list of metabolites in the reduced model
	_toReducedMet.clear(); // clear in case of repeated use
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		_toReducedMet[met] = _reduced->createMetabolite("t_"+met->getName());
	}

	boost::unordered_map<SCIP_VAR*, ReactionPtr> aggrReactions;
	_toReducedRxn.clear(); // clear in case of repeated use
	_toOriginalRxn.clear(); // clear in case of repeated use
	// include only those reactions in the reduced model that are not aggregated
	foreach(ReactionPtr rxn, _model->getReactions()) {
		if(scip->hasFluxVar(rxn)) {
			SCIP_VAR* var = scip->getFlux(rxn);
			if(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL) {
				// look at transformed var
				var = SCIPvarGetTransVar(var);
			}
			double scalar = 1;
			while(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED) {
				// var is aggregated, so the reaction is fully coupled to another reaction.

				// the presolver will make sure that the bounds are correctly maintained,
				// hence, we will not think about bounds.
				// however we cannot fetch the correct stoichiometric coefficients, so for this we have to do the backtrace.
				scalar *= SCIPvarGetAggrScalar(var);
				var = SCIPvarGetAggrVar(var);
			}
			// var is now the root-var that is used in the presolved problem.
			ReactionPtr aggrrxn;
			if(aggrReactions.find(var) == aggrReactions.end()) {
				aggrrxn = _reduced->createReaction("t_"+rxn->getName());
				aggrrxn->setLb(SCIPvarGetLbGlobal(var));
				aggrrxn->setUb(SCIPvarGetUbGlobal(var));
				aggrrxn->setObj(SCIPvarGetObj(var));
				_reducedScip->setFlux(aggrrxn, var);
				aggrReactions[var] = aggrrxn;
				_toOriginalRxn[aggrrxn] = rxn;
			}
			else {
				aggrrxn = aggrReactions[var];
			}
#ifndef NDEBUG
			if(aggrrxn->getName().compare("t_IMPDH") == 0) {
				std::cout << aggrrxn->getName() << " <- " << rxn->getName() << std::endl;
			}
#endif

			if(!aggrrxn->isExchange()) aggrrxn->setExchange(rxn->isExchange()); // if one of the aggregated reactions is an exchange reaction, the aggregated reaction is an exchange reaction
			if(!aggrrxn->isProblematic()) aggrrxn->setProblematic(rxn->isProblematic());// if one of the aggregated reactions is problematic, the aggregated reaction is problematic
			_toReducedRxn[rxn] = aggrrxn;
			foreach(Stoichiometry s, rxn->getStoichiometries()) {
				MetabolitePtr met = _toReducedMet[s.first];
				double stoich = aggrrxn->getStoichiometry(met);
				stoich += scalar*s.second;
				aggrrxn->setStoichiometry(met, stoich);
			}
		}
	}

#ifndef NDEBUG
	//foreach(ReactionPtr aggrrxn, _reduced->getReactions()) {
	{
		ReactionPtr aggrrxn = _reduced->getReaction("t_IMPDH");
		std::cout << aggrrxn->getName() << "("<< aggrrxn->getLb() << "," << aggrrxn->getUb() << "," << aggrrxn->getObj() << (aggrrxn->isExchange()?",ex":"")<< "): ";
		foreach(ReactionPtr rxn, _model->getReactions()) {
			if(_toReducedRxn[rxn] == aggrrxn) {
				std::cout << rxn->getName() << "("<<rxn->getLb() << "," << rxn->getUb() << "," << rxn->getObj() << (rxn->isExchange()?",ex":"")<< ") ";
			}
		}
		std::cout << std::endl;
		foreach(Stoichiometry s, aggrrxn->getStoichiometries()) {
			std::cout << s.second << "*" << s.first->getName() << " ";
		}
		std::cout << std::endl;
	}
#endif

	// init helper variables
	_cycle_find = LPFluxPtr( new LPFlux(_reduced, false));
	_cycle_find->setObjSense(false); // minimize
	_cycle_test = LPFluxPtr( new LPFlux(_reduced, false));
	_cycle_test->setObjSense(true); // maximize
	_flux_simpl = LPFluxPtr( new LPFlux(_reduced, true));
	_is_find = DualPotentialsPtr( new DualPotentials(_model)); //I cannot use the reduced model here, because I would lose infeasible sets
	_pot_test = LPPotentialsPtr( new LPPotentials(_model)); // it doesn't make sense to use the reduced model, since this is only used for testing
#else

	// even if we do not use the results of the presolver to directly simplify the model, we can use that to infer coupling relations
	if(_coupling.use_count() == 0) {
		_coupling = CouplingPtr(new Coupling()); // if not supplied with external coupling information, start from scratch
	}

#ifndef FINDBUG
	// this map maps the transformed variables to their original reaction
	unordered_map<SCIP_VAR*, ReactionPtr> toOriginal;
	foreach(ReactionPtr rxn, _model->getReactions()) {
		if(scip->hasFluxVar(rxn)) {
			SCIP_VAR* var = scip->getFlux(rxn);
			if(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL) {
				// look at transformed var
				var = SCIPvarGetTransVar(var);
				toOriginal[var] = rxn;
			}
		}
	}
	foreach(ReactionPtr rxn, _model->getReactions()) {
		if(scip->hasFluxVar(rxn)) {
			SCIP_VAR* var = scip->getFlux(rxn);
			if(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL) {
				// look at transformed var
				var = SCIPvarGetTransVar(var);
			}
			double scalar = 1;
			double constant = 0;
			while(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED) {
				// var is aggregated, so the reaction is coupled to another reaction.
				// it may not be fully coupled, because existing flows must not be neglected
				// but even in those cases, we can infer directional couplings

				// if y is the aggregation variable, we have
				// var = a*y + c
				// if x is the aggregation variable of y, we get
				// var = a_1*(a_2*x + c_2) + c_1 = a_1*a_2*x + a_1*c_2 + c_1
				constant += scalar * SCIPvarGetAggrConstant(var);
				scalar *= SCIPvarGetAggrScalar(var);
				var = SCIPvarGetAggrVar(var);

				// if the transformed variable can be associated to a reaction, we can infer coupling information
				unordered_map<SCIP_VAR*, ReactionPtr>::iterator iter = toOriginal.find(var);
				if(iter != toOriginal.end()) {
					ReactionPtr other = iter->second;
					if(constant > -EPSILON) { // also include the case where constant == 0
						if(scalar > 0) { // TODO: numerical troubles?
							// positive flux on other implies positive flux on rxn
							_coupling->addCoupled(DirectedReaction(other, true), DirectedReaction(rxn, true));
							// negative flux un rxn implies negative flux on other
							_coupling->addCoupled(DirectedReaction(rxn, false), DirectedReaction(other, false));
						}
						else {
							// negative flux on other implies positive flux on rxn
							_coupling->addCoupled(DirectedReaction(other, false), DirectedReaction(rxn, true));
							// negative flux on rxn imples positive flux on other
							_coupling->addCoupled(DirectedReaction(rxn, false), DirectedReaction(other, true));

						}
					}
					if(constant < EPSILON) { // also include the case where constant == 0
						if(scalar > 0) { // TODO: numerical troubles?
							// negative flux on other implies negative flux on rxn
							_coupling->addCoupled(DirectedReaction(other, false), DirectedReaction(rxn, false));
							// positive flux un rxn implies positive flux on other
							_coupling->addCoupled(DirectedReaction(rxn, true), DirectedReaction(other, true));
						}
						else {
							// positive flux on other implies negative flux on rxn
							_coupling->addCoupled(DirectedReaction(other, true), DirectedReaction(rxn, false));
							// positive flux on rxn imples negative flux on other
							_coupling->addCoupled(DirectedReaction(rxn, true), DirectedReaction(other, false));

						}
					}
				}
			}
		}
	}
#endif
	_coupling->computeClosure();

#endif

	return SCIP_OKAY;
}

SCIP_RETCODE ThermoConstraintHandler::scip_trans(
		SCIP*              scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
		SCIP_CONS*         sourcecons,         /**< source constraint to transform */
		SCIP_CONS**        targetcons          /**< pointer to store created target constraint */
		)
{
	// constraint data does not need to be transformed, we can use the original
    // create new constraint with a value copy of the constraint data
	// we have to do a value copy, else we get problems at deletion
	SCIP_ConsData* data = SCIPconsGetData(sourcecons);
	SCIP_ConsData* copy = NULL;
	if(data != NULL) {
		copy = new SCIP_ConsData(PotSpaceConstraintPtr(data->val));
	}
	SCIP_CALL(  SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, copy,
            SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
            SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
            SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
            SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons))  );
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

void ThermoConstraintHandler::addPotSpaceConstraint(PotSpaceConstraintPtr psc, SCIP_NODE* node) {
#ifndef FINDBUG
	SCIP_CONS* cons;
	SCIP* scip = getScip()->getScip();
	SCIP_CONSHDLR* hdlr = SCIPfindConshdlr(scip, CONSTRAINT_NAME);

	//SCIP_ConsData* data = new SCIP_ConsData(psc); // will be freed by ThermoConstraintHandler::scip_delete
	SCIP_ConsData* data = new SCIP_ConsData(PotSpaceConstraintPtr(psc));

	SCIPcreateCons(scip, &cons, "potSpaceCon", hdlr, data, true, true, true, false, true, true, false, false, false, true);
	SCIPaddConsNode(scip, node, cons, NULL);
	SCIPreleaseCons(scip, &cons);
#endif
}

SCIP_RETCODE ThermoConstraintHandler::scip_delete(
		SCIP* 				scip,               /**< SCIP data structure */
		SCIP_CONSHDLR*		conshdlr,           /**< the constraint handler itself */
		SCIP_CONS*			cons,				/**< constraint that will be deleted */
		SCIP_CONSDATA**		data				/**< constraint data to free */
) {
	if(*data != NULL) {
		// delete the pointer to the shared pointer
		delete *data;
		*data = NULL;
	}
	return SCIP_OKAY;
}


} /* namespace metaopt */
