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
 * PotBoundPropagation.cpp
 *
 *  Created on: 18.05.2012
 *      Author: arnem
 */

#include <boost/unordered_set.hpp>
#include <iostream>
#include "PotBoundPropagation2.h"
#include "metaopt/model/Reaction.h"
#include "metaopt/Properties.h"

#include "objscip/objscip.h"
#include "metaopt/scip/ScipError.h"

#define PRINT_UPDATE

using namespace boost;
using namespace std;

#define EPSILON 0.0001

namespace metaopt {

PotBoundPropagation2::PotBoundPropagation2(ModelPtr model) :
	_model(model)
{
	init_Queue();
	//int i = 0;
#if 0
	while(!_queue.empty()) {
		ArcEvent a = _queue.top();
		_queue.pop();
		update(a.arc);
		//i++;
		//if(i > 100) return;
	}
#else
	foreach(MetabolitePtr m, _model->getMetabolites()) {
		_maxBounds[m]->_bound = m->getPotUb();
		_minBounds[m]->_bound = -m->getPotLb();
	}
#endif

	// for later update steps, initialize the _vars map
	int i = 0;
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		_vars[_maxBounds[met]] = i++;
		_vars[_minBounds[met]] = i++;
	}

}

PotBoundPropagation2::~PotBoundPropagation2() {
	// TODO Auto-generated destructor stub
}


void PotBoundPropagation2::init_Queue() {
	// create bounds (nodes of the Leontief-Hypergraph)
	foreach(MetabolitePtr m, _model->getMetabolites()) {
		MetBoundPtr max(new MetBound(m, false));
		_maxBounds[m] = max;
		MetBoundPtr min(new MetBound(m, true));
		_minBounds[m] = min;
	}

	// create incidencies (only meaningful for internal creations)
	foreach(ReactionPtr r, _model->getInternalReactions()) {
		//assert(r->getProducts().size() > 0 && r->getReactants().size() > 0);
		if(r->canFwd()) { // reaction can operate in forward direction.
			foreach(Stoichiometry v, r->getStoichiometries()) {
				ArcPtr a(new Arc());
				_arcs.push_back(a);
				a->_creator = r;
				a->_fwdcreator = true;
				double vcoef;
				if(v.second > 0) {
					a->_target = _maxBounds.at(v.first);
					vcoef = 1/v.second;
				}
				else {
					a->_target = _minBounds.at(v.first);
					vcoef =-1/v.second; // make sure its positive
				}
				assert(vcoef > 0);
				foreach(Stoichiometry x, r->getProducts()) {
					if(v.first != x.first) {
						assert(x.second > 0);
						pair<MetBoundPtr, double> p(_minBounds[x.first],x.second*vcoef);
						a->_input.push_back(p);
						_incidence[_minBounds.at(x.first)].push_back(a);
					}
				}
				foreach(Stoichiometry x, r->getReactants()) {
					if(v.first != x.first) {
						assert(x.second > 0);
						pair<MetBoundPtr, double> p(_maxBounds[x.first],x.second*vcoef);
						a->_input.push_back(p);
						_incidence[_maxBounds.at(x.first)].push_back(a);
					}
				}
			}
		}
		if(r->canBwd()) { // reaction can operate in backward direction
			foreach(Stoichiometry v, r->getStoichiometries()) {
				ArcPtr a(new Arc());
				_arcs.push_back(a);
				a->_creator = r;
				a->_fwdcreator = false;
				double vcoef;
				if(v.second > 0) {
					a->_target = _minBounds.at(v.first); // since reversed reaction
					vcoef = 1/v.second;
				}
				else {
					a->_target = _maxBounds.at(v.first); // since reversed reaction
					vcoef =-1/v.second; // make sure its positive
				}
				assert(vcoef > 0);
				foreach(Stoichiometry x, r->getReactants()) { // the products of the reversed reaction are the reactants of the original reaction
					if(v.first != x.first) {
						assert(x.second > 0);
						pair<MetBoundPtr, double> p(_minBounds[x.first],x.second*vcoef);
						a->_input.push_back(p);
						_incidence[_minBounds.at(x.first)].push_back(a);
					}
				}
				foreach(Stoichiometry x, r->getProducts()) { // the reactants of the reversed reaction are the products of the original reaction
					if(v.first != x.first) {
						assert(x.second > 0);
						pair<MetBoundPtr, double> p(_maxBounds[x.first],x.second*vcoef);
						a->_input.push_back(p);
						_incidence[_maxBounds.at(x.first)].push_back(a);
					}
				}
			}
		}
	}

	foreach(ReactionPtr r, _model->getReactions()) {
		if(r->isExchange()) {
			if(r->canFwd()) { // reaction in forward direction
				foreach(Stoichiometry s, r->getProducts()) {
					MetBoundPtr met = _maxBounds.at(s.first);
					met->_bound = s.first->getPotUb();
					met->_isBoundary = true;
					//cout << "boundary met (ub)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
				foreach(Stoichiometry s, r->getReactants()) {
					MetBoundPtr met = _minBounds.at(s.first);
					met->_bound = -s.first->getPotLb(); // _bound stores transformed bound!
					met->_isBoundary = true;
					//cout << "boundary met (lb)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
			}
			if(r->canBwd()) { // reaction in backward direction
				foreach(Stoichiometry s, r->getReactants()) {
					MetBoundPtr met = _maxBounds.at(s.first);
					met->_bound = s.first->getPotUb();
					met->_isBoundary = true;
					//cout << "boundary met (ub)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
				foreach(Stoichiometry s, r->getProducts()) {
					MetBoundPtr met = _minBounds.at(s.first);
					met->_bound = -s.first->getPotLb(); // _bound stores transformed bound!
					met->_isBoundary = true;
					//cout << "boundary met (lb)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
			}
		}
	}
}

PotBoundPropagation2::Arc::Arc()
	: _active(true) {
	// nothing else
}

int PotBoundPropagation2::Arc::unknown_inputs() const {
	int res = 0;
	for(unsigned int i = 0; i < _input.size(); i++) {
		if(isinf(_input[i].first->_bound)) {
			res++;
		}
	}
	return res;
}

double PotBoundPropagation2::Arc::update_value() const {
	double res = 0;
	for(unsigned int i = 0; i < _input.size(); i++) {
		int val;
		if(isinf(_input[i].first->_bound)) {
			val = _input[i].first->getMax();
		}
		else {
			val = _input[i].first->_bound;
		}
		assert(_input[i].second > 0);
		res += _input[i].second*val;
	}
	return min(_target->getMax(), res);
}

PotBoundPropagation2::ArcPtr PotBoundPropagation2::getReversed(const ArcPtr a) const {
	ArcPtr out(new Arc());
	out->_creator = a->_creator;
	out->_fwdcreator = !a->_fwdcreator;
	if(a->_target->_isMinBound) {
		out->_target = _maxBounds.at(a->_target->_met);
	}
	else {
		out->_target = _minBounds.at(a->_target->_met);
	}
	for(unsigned int i = 0; i < a->_input.size(); i++) {
		MetBoundPtr m = a->_input[i].first;
		double coef = a->_input[i].second;
		if(m->_isMinBound) {
			pair<MetBoundPtr, double> p(_maxBounds.at(m->_met), coef);
			out->_input.push_back(p);
		}
		else {
			pair<MetBoundPtr, double> p(_minBounds.at(m->_met), coef);
			out->_input.push_back(p);
		}
	}
	return out;
}

void PotBoundPropagation2::update(ArcPtr a) {
	if(a->_target->_bound + EPSILON < a->update_value()) {
		// only in this case we have to do updates
		cout << a->_creator->getName() << " updating " << a->_target->_met->getName() << (a->_target->_isMinBound? " (min)" : " (max)") << " from " << a->_target->_bound << " to " << a->update_value() << " of at most " << a->_target->getMax() << " using " << a->unknown_inputs() << " unknown metabolite bounds" << endl;
		a->_target->_bound = a->update_value();
		a->_target->_last_update = a;
		// propagate
		foreach(ArcPtr b, _incidence[a->_target]) {
			ArcEvent be(b);
			if(be.gain_val < b->_target->getMax() - b->_target->_bound) {
				// only add to queue if the update may have an effect
				_queue.push(be);
			}
		}
	}
}

void PotBoundPropagation2::print() {
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		MetBoundPtr lb = _minBounds[met];
		MetBoundPtr ub = _maxBounds[met];


		if(isinf(lb->_bound) || isinf(ub->_bound)) {
			cout << met->getName() << ": "<< met->getPotLb() << " " << lb->orig_value() << " " << ub->orig_value() << " " << met->getPotUb() << "\t\t\t";
			boost::shared_ptr<vector<ReactionPtr> > producers = met->getProducers();
			foreach(ReactionPtr r, *producers) {
				cout << r->getName() << " ";
			}
			cout << "-> ";
			boost::shared_ptr<vector<ReactionPtr> > consumers = met->getConsumers();
			foreach(ReactionPtr r, *consumers) {
				cout << r->getName() << " ";
			}
			cout << endl;
		}

		if(met->getPotLb() +EPSILON < lb->orig_value()) {
			cout << met->getName() << ": original="<< met->getPotLb() << " induced=" << lb->orig_value() << endl << endl;
			/*
			ArcPtr src = lb->_last_update.lock();
			cout << "sourced by: " << src->_creator->getName() << endl;
			for(unsigned int i = 0; i < src->_input.size(); i++) {
				cout <<  "\t" << (src->_input[i].first->_isMinBound?"+":"-")<< src->_input[i].second << "*" << src->_input[i].first->orig_value() << "\t(" << src->_input[i].first->_met->getName() << ")" << endl;
			}
			*/
		}
		if(met->getPotUb() -EPSILON > ub->orig_value()) {
			cout << met->getName() << ": original="<< met->getPotUb() << " induced=" << ub->orig_value() << endl << endl;
			/*
			ArcPtr src = ub->_last_update.lock();
			cout << "sourced by: " << src->_creator->getName() << endl;
			for(unsigned int i = 0; i < src->_input.size(); i++) {
				cout <<  "\t" << (src->_input[i].first->_isMinBound?"-":"+")<< src->_input[i].second << "*" << src->_input[i].first->orig_value() << "\t(" << src->_input[i].first->_met->getName() << ")" << endl;
			}
			*/
		}
	}
}

shared_ptr<vector<pair<ReactionPtr,bool> > > PotBoundPropagation2::getBlockedReactions() {
	shared_ptr<vector<pair<ReactionPtr,bool> > > result(new vector<pair<ReactionPtr, bool> >());
	foreach(ReactionPtr r, _model->getReactions()) {
		if(!r->isExchange()) {
			double val_fwd = 0;
			double val_bwd = 0;
			foreach(Stoichiometry s, r->getProducts()) {
				val_fwd += s.second * _minBounds.at(s.first)->orig_value();
				val_bwd -= s.second * _maxBounds.at(s.first)->orig_value();
			}
			foreach(Stoichiometry s, r->getReactants()) {
				val_fwd -= s.second * _maxBounds.at(s.first)->orig_value(); // coefficients are all positive
				val_bwd += s.second * _minBounds.at(s.first)->orig_value(); // coefficients are all positive
			}
			//cout << r->getName() << ": " << val_fwd << " / " << val_bwd << endl;
			if(val_fwd > -EPSILON) {
				cout << "blocked reaction fwd: " << r->getName() << " by " << val_fwd << endl;
				pair<ReactionPtr, bool> p(r,true);
				result->push_back(p);
			}
			if(val_bwd > -EPSILON) {
				cout << "blocked reaction bwd: " << r->getName() << " by " << val_bwd << endl;
				pair<ReactionPtr, bool> p(r,false);
				result->push_back(p);
			}
		}
	}
	return result;
}

void PotBoundPropagation2::buildHardArcConstraint(ArcPtr a, SCIP_LPI* lpi, bool fwd) {
	double lhs, rhs;
	if(fwd) {
		lhs = -INFINITY;
		rhs = 0;
	}
	else {
		lhs = 0;
		rhs = INFINITY;
	}
	int beg = 0;
	vector<int> ind;
	vector<double> coef;
	ind.push_back(_vars.at(a->_target));
	coef.push_back(1);

//	cout << a->_target->_met->getName()<<(a->_target->_isMinBound?"(min)":"(max)") << (fwd?" < ":" > ");
	for(unsigned int j = 0; j < a->_input.size(); j++) {
		MetBoundPtr met = a->_input[j].first;
		ind.push_back(_vars.at(met));
		coef.push_back(-a->_input[j].second);
//		cout << a->_input[j].second << "*" << a->_input[j].first->_met->getName() << (a->_input[j].first->_isMinBound?"(min) ":"(max) ");
	}
//	cout << endl;
	BOOST_SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
}

void PotBoundPropagation2::updateStepHard(ScipModelPtr scip) {
#if 0
	SCIP_LPI* lpi;
	BOOST_SCIP_CALL( SCIPlpiCreate(&lpi, "potboundprop_updateStepHard", SCIP_OBJSEN_MAXIMIZE) );

	// create vars - we do this by creating empty columns, since this guarantees that the var exists
	foreach(VarEntry e, _vars) {
		// just set to arbitrary values
		double lb = -INFINITY;
		double ub = INFINITY;
		double obj = 1;
		BOOST_SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL,0,0,0,0) );
	}

	// create constraints of fixed reaction directions
	shared_ptr<unordered_set<ReactionPtr> > dirs = scip->getFixedDirections();
	cout << "fixed dirs: ";
	foreach(ReactionPtr rxn, *dirs) {
		cout << rxn->getName() << " ";
	}
	cout << endl;
	// first translate reactions to arcs
	vector<ArcPtr> fixedArcs;
	foreach(ArcPtr a, _arcs) {
		ReactionPtr r = a->_creator;
		if(dirs->find(r) != dirs->end()) {
			if((scip->getCurrentFluxLb(r) < -EPSILON) == a->_fwdcreator) {
				// reaction is fixed to reverse
				// for a = (v,A) build a constraint of the form
				// 0 \leq \mu(v) - sum_{x \in A} S_{xa}/S_{va} \mu(x)
				// the following line is disabled, because if we want to add this constraints, we will have to maximize for each metabolite separately (or prove that it is not necessary)
#if 0
				buildHardArcConstraint(a, lpi, false);
#endif
				// also do this on the reversed arc, because that involves other bound-variables
				ArcPtr rev = getReversed(a);
				//buildHardArcConstraint(rev, lpi, true);
				fixedArcs.push_back(rev);
			}
			else {
				// reaction is fixed to forward
				// for a = (v,A) build a constraint of the form
				// \mu(v) - sum_{x \in A} S_{xa}/S_{va} \mu(x) \leq 0
				//buildHardArcConstraint(a, lpi, true);
				fixedArcs.push_back(a);

				//the following lines are disabled, because if we want to add this constraints, we will have to maximize for each metabolite separately (or prove that it is not necessary)
#if 0
				// also do this on the reversed arc, because that involves other bound-variables
				ArcPtr rev = getReversed(a);
				buildHardArcConstraint(rev, lpi, false);
#endif
			}
		}
	}


	// variables with bound of -inf make trouble
	// hence we propagate those manually
	bool propagatedInf = true;
	while(propagatedInf) {
		propagatedInf = false;
		foreach(ArcPtr a, fixedArcs) {
			if(!isinf(a->_target->_bound) || a->_target->_bound > 0) {
				for(unsigned int i = 0; i < a->_input.size(); i++) {
					double bound = a->_input[i].first->_bound;
					if(isinf(bound) && bound < 0) {
#ifdef PRINT_UPDATE
						cout << "updating (hard)" << a->_target->_met->getName() << (a->_target->_isMinBound?" (min)":" (max)") << " to " << "-inf" << " was " << a->_target->_bound << endl;
#endif
						a->_target->_bound = -INFINITY;
						propagatedInf = true;
					}
				}
			}
		}
	}

	// set var-bounds of vars that need no updates anymore
	// these vars are not used in any constraints, the above step made sure of that
	// also set the bounds of the other variables (so thats in correct order)
	foreach(VarEntry e, _vars) {
		if(isinf(e.first->_bound) && e.first->_bound < 0) {
			double lb = 0;
			double ub = 0;
			double obj = 0;
			int ind = e.second;
			BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
			BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
		}
		else {
			double lb = -INFINITY;
			double ub = e.first->_bound;
			double obj = 1;
			int ind = e.second;
			BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
			BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
		}
	}

	// now create the constraints
	foreach(ArcPtr a, fixedArcs) {
		if(!isinf(a->_target->_bound) || a->_target->_bound > 0) {
			buildHardArcConstraint(a, lpi, true);
		}
	}

	//SCIPlpiWriteLP(lpi, "debug.lp");

	// solve
	BOOST_SCIP_CALL( SCIPlpiSolveDual(lpi) );
	assert( SCIPlpiWasSolved(lpi) );
	if(SCIPlpiIsPrimalFeasible(lpi)) {
		if(SCIPlpiIsPrimalUnbounded(lpi)) {
			cout << "WARNING: unbounded!!! probably will get wrong results!" << endl;
		}
		// we have a feasible solution
		double objval;
		double primsol[_vars.size()];
		BOOST_SCIP_CALL( SCIPlpiGetSol(lpi,&objval, primsol, NULL, NULL, NULL) );
		// perform update
		foreach(VarEntry e, _vars) {
			if(primsol[e.second] < e.first->_bound - EPSILON) {
#ifdef PRINT_UPDATE
				cout << "updating (hard) " << e.first->_met->getName() << (e.first->_isMinBound?" (min)":" (max)") << " to " << primsol[e.second] << " was " << e.first->_bound << endl;
#endif
				e.first->_bound = primsol[e.second];
			}
		}
	}
	else {
		// if we want to account for this we may have to add complicated constraints
		// dunno how to do it, so just ignore
		cout << "WARNING: pot bound propagation found inconsisitency with direction fixations, which will be ignored" << endl;
	}
#else
	SCIP_LPI* lpi;
	BOOST_SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "potboundprop_updateStepHard", SCIP_OBJSEN_MAXIMIZE) );

	unordered_map<MetabolitePtr, int> met_idx;
	int i = 0;
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		double lb = -_minBounds.at(met)->_bound;
		double ub = _maxBounds.at(met)->_bound;
		// actually we will have to deal with cases, where lb > ub (but this doesn't happen on the current instances, so don't think about it yet)
		double obj = 0;
		met_idx[met] = i;
		BOOST_SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL,0,0,0,0) );
		i++;
	}

	// create constraints of fixed reaction directions
	shared_ptr<unordered_set<ReactionPtr> > dirs = scip->getFixedDirections();
	cout << "fixed dirs: ";
	foreach(ReactionPtr rxn, *dirs) {
		if(!rxn->isExchange()) {
			cout << rxn->getName() << " ";
			double lhs = -INFINITY;
			double rhs = INFINITY;
			if(scip->getCurrentFluxLb(rxn) < -EPSILON) {
				// reaction is fixed to reverse
				lhs = 0;
			}
			else {
				// reaction is fixed to forward
				rhs = 0;
			}
			int beg = 0;
			vector<int> ind;
			vector<double> coef;
			foreach(Stoichiometry s, rxn->getStoichiometries()) {
				ind.push_back(met_idx.at(s.first));
				coef.push_back(s.second);
			}
			BOOST_SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
		}
	}
	cout << endl;

	// now run variability analysis on metabolite potentials
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		int ind = met_idx.at(met);
		double obj = 1;
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj));
		BOOST_SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
		assert( SCIPlpiWasSolved(lpi) );
		assert( SCIPlpiIsPrimalFeasible(lpi) );
		double objval;
		SCIPlpiGetObjval(lpi, &objval);
		if(_maxBounds.at(met)->_bound > objval) {
			cout << "updating (hard) " << met->getName() << " (max) to " << objval << " was " << _maxBounds.at(met)->_bound << endl;
		}
		_maxBounds.at(met)->_bound = objval;
		obj = 0;
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj));
	}
	// and do the same for min-bounds
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		int ind = met_idx.at(met);
		double obj = -1;
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj));
		BOOST_SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
		assert( SCIPlpiWasSolved(lpi) );
		assert( SCIPlpiIsPrimalFeasible(lpi) );
		double objval; // we don't have to multiply objval by -1, because the objective function already does this for us
		SCIPlpiGetObjval(lpi, &objval);
		if(_minBounds.at(met)->_bound > objval) {
			cout << "updating (hard) " << met->getName() << " (min) to " << objval << " was " << _minBounds.at(met)->_bound << endl;
		}
		_minBounds.at(met)->_bound = objval;
		obj = 0;
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj));
	}

#endif
}

bool PotBoundPropagation2::updateStepFlow() {
	// compute X set
	unordered_set<MetBoundPtr> X;
#if 1
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		// add all bounds that are not already -inf
		MetBoundPtr bound = _maxBounds[met];
		if(!bound->_isBoundary && (!isinf(bound->_bound) || bound->_bound > 0)) X.insert(bound);
		bound = _minBounds[met];
		if(!bound->_isBoundary && (!isinf(bound->_bound) || bound->_bound > 0)) X.insert(bound);
	}
	foreach(ArcPtr a, _arcs) {
		if(a->_active) {
			double res = 0;
			for(unsigned int i = 0; i < a->_input.size(); i++) {
				res += a->_input[i].second*a->_input[i].first->_bound;
			}
			if(res > a->_target->_bound) {
				X.erase(a->_target);
			}
		}
	}
#endif
	cout << "X: ";
	foreach(MetBoundPtr m, X) {
		cout << m->_met->getName()<< (m->_isMinBound?"(min)":"(max)");
	}
	cout << endl;
	// X is now computed
	// we now have to create the LP for solving the update step

	SCIP_LPI* lpi;
	BOOST_SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "potboundprop_updateStepFlow", SCIP_OBJSEN_MINIMIZE) );

	// create vars - we do this by creating empty columns, since this guarantees that the var exists
	foreach(VarEntry e, _vars) {
		double lb = -INFINITY;
		double ub = INFINITY;
		double obj = 1;
		BOOST_SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL,0,0,0,0) );
	}
	// set var-bounds of vars that need no updates anymore
	foreach(VarEntry e, _vars) {
		if(isinf(e.first->_bound) && e.first->_bound < 0) {
			double lb = 0;
			double ub = 0;
			double obj = 0;
			int ind = e.second;
			BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
			BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
		}
		else if(e.first->_isBoundary) {
			// we can't improve bounds of boundary metabolites.
			double lb = e.first->_bound;
			double ub = e.first->_bound;
			double obj = 0;
			int ind = e.second;
			BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
			BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
		}
	}


	// create constraints for arcs
	foreach(ArcPtr a, _arcs) {
		if(!a->_target->_isBoundary && a->_active && (!isinf(a->_target->_bound) || a->_target->_bound > 0)) { // else, updating is no use
			// for a = (v,A) build a constraint of the form
			// sum_{x \in A \setminus X} S_{xa}/S_{va} _bound(x) <= \mu(v) - sum_{x \in A \cap X} S_{xa}/S_{va} \mu(x)
			double rhs = INFINITY;
			double lhs = 0; // start with zero and increase with metabolites not in X in constraint building process
			int beg = 0;
			vector<int> ind;
			vector<double> coef;
			ind.push_back(_vars.at(a->_target));
			coef.push_back(1);
			for(unsigned int j = 0; j < a->_input.size(); j++) {
				MetBoundPtr met = a->_input[j].first;
				if(X.find(met) != X.end()) {
					ind.push_back(_vars.at(met));
					coef.push_back(-a->_input[j].second);
				}
				else {
					lhs += a->_input[j].second * met->_bound;
				}
			}
			BOOST_SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
		}
	}
	BOOST_SCIP_CALL( SCIPlpiSolveDual(lpi) );
	assert( SCIPlpiWasSolved(lpi) );
	assert( SCIPlpiIsPrimalFeasible(lpi) );
	if(SCIPlpiIsPrimalUnbounded(lpi) ) {
		// can this happen?
		// yes, it can if metabolites are not reachable.
		assert( SCIPlpiExistsPrimalRay(lpi) );
		double ray[_vars.size()];
		BOOST_SCIP_CALL( SCIPlpiGetPrimalRay(lpi, ray) );
		// perform update
		bool updated = false;
		foreach(VarEntry e, _vars) {
			if(ray[e.second] < - EPSILON) { // every ray will only have nonpositive entries
#ifdef PRINT_UPDATE
				cout << "updating " << e.first->_met->getName() << (e.first->_isMinBound?" (min)":" (max)") << " to " << "-inf" << " was " << e.first->_bound << endl;
#endif
				e.first->_bound = -INFINITY;
				foreach(ArcPtr inc, _incidence.at(e.first)) {
					inc->_active = false;
				}
				updated = true;
			}
		}
		return updated;
	}
	else {
		// we have a feasible solution
		double objval;
		double primsol[_vars.size()];
		BOOST_SCIP_CALL( SCIPlpiGetSol(lpi,&objval, primsol, NULL, NULL, NULL) );
		// perform update
		bool updated = false;
		foreach(VarEntry e, _vars) {
			if(primsol[e.second] < e.first->_bound - EPSILON) {
#ifdef PRINT_UPDATE
				cout << "updating " << e.first->_met->getName() << (e.first->_isMinBound?" (min)":" (max)") << " to " << primsol[e.second] << " was " << e.first->_bound << endl;
#endif
				e.first->_bound = primsol[e.second];
				if(e.first->_isMinBound) {
					if(_maxBounds.at(e.first->_met)->_bound + e.first->_bound < 0) {
						// lb > ub !
#ifdef PRINT_UPDATE
						cout << "update produced inconsistency  in all-flow" << endl;
#endif
						foreach(ArcPtr inc, _incidence.at(e.first)) {
							inc->_active = false;
						}
					}
				}
				else {
					if(_minBounds.at(e.first->_met)->_bound + e.first->_bound < 0) {
						// lb > ub !
#ifdef PRINT_UPDATE
						cout << "update produced inconsistency  in all-flow" << endl;
#endif
						foreach(ArcPtr inc, _incidence.at(e.first)) {
							inc->_active = false;
						}
					}
				}
				updated = true;
			}
		}
		return updated;
	}
}

void PotBoundPropagation2::update(ScipModelPtr scip) {
	// for now: always compute from scratch
	foreach(MetabolitePtr m, _model->getMetabolites()) {
		_maxBounds[m]->_bound = scip->getCurrentPotentialUb(m);
		_minBounds[m]->_bound = -scip->getCurrentPotentialLb(m);

		//+inf makes trouble, so replace by really big bounds
		if(_maxBounds[m]->_bound > 100000) _maxBounds[m]->_bound = 100000;
		if(_minBounds[m]->_bound > 100000) _minBounds[m]->_bound = 100000;
	}
	updateStepHard(scip);
	bool updated = true;
	while(updated) {
		updated = updateStepFlow(); // updateStepFlow may not give best update and may require repetitive calls
		if(updated) updateStepHard(scip); // updateStepHard always gives best update
	}
}

void PotBoundPropagation2::update() {
	bool updated = true;
	while(updated) {
		updated = updateStepFlow(); // updateStepFlow may not give best update and may require repetitive calls
	}
}

} /* namespace metaopt */
