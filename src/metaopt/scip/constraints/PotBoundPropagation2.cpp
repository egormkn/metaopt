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


using namespace boost;
using namespace std;

namespace metaopt {

PotBoundPropagation2::PotBoundPropagation2(ModelPtr model) :
	_model(model)
{
	init_Queue();
	//int i = 0;
	while(!_queue.empty()) {
		ArcEvent a = _queue.top();
		_queue.pop();
		update(a.arc);
		//i++;
		//if(i > 100) return;
	}

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

	// create incidencies
	foreach(ReactionPtr r, _model->getReactions()) {
		if(r->getUb() > EPSILON) { // reaction can operate in forward direction.
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
		if(r->getLb() < -EPSILON) { // reaction can operate in backward direction
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
			if(r->getUb() > EPSILON) { // reaction in forward direction
				foreach(Stoichiometry s, r->getProducts()) {
					MetBoundPtr met = _maxBounds.at(s.first);
					met->_bound = s.first->getPotUb();
					cout << "boundary met (ub)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
				foreach(Stoichiometry s, r->getReactants()) {
					MetBoundPtr met = _minBounds.at(s.first);
					met->_bound = -s.first->getPotLb(); // _bound stores transformed bound!
					cout << "boundary met (lb)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
			}
			if(r->getLb() < -EPSILON) { // reaction in backward direction
				foreach(Stoichiometry s, r->getReactants()) {
					MetBoundPtr met = _maxBounds.at(s.first);
					met->_bound = s.first->getPotUb();
					cout << "boundary met (ub)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
				foreach(Stoichiometry s, r->getProducts()) {
					MetBoundPtr met = _minBounds.at(s.first);
					met->_bound = -s.first->getPotLb(); // _bound stores transformed bound!
					cout << "boundary met (lb)" << s.first->getName() << endl;
					foreach(ArcPtr a, _incidence[met]) {
						_queue.push(ArcEvent(a));
					}
				}
			}
		}
	}
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
	for(int i = 0; i < a->_input.size(); i++) {
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
		/*
		cout << met->getName() << ": "<< met->getPotLb() << " " << lb->orig_value() << " " << ub->orig_value() << " " << met->getPotUb() << "\t\t\t";

		if(isinf(lb->_bound) || isinf(ub->_bound)) {
			boost::shared_ptr<vector<ReactionPtr> > producers = met->getProducers();
			foreach(ReactionPtr r, *producers) {
				cout << r->getName() << " ";
			}
			cout << "-> ";
			boost::shared_ptr<vector<ReactionPtr> > consumers = met->getConsumers();
			foreach(ReactionPtr r, *consumers) {
				cout << r->getName() << " ";
			}
		}
		cout << endl;*/

		if(met->getPotLb() +EPSILON < lb->orig_value()) {
			cout << met->getName() << ": original="<< met->getPotLb() << " induced=" << lb->orig_value() << endl;
			ArcPtr src = lb->_last_update.lock();
			cout << "sourced by: " << src->_creator->getName() << endl;
			for(unsigned int i = 0; i < src->_input.size(); i++) {
				cout <<  "\t" << (src->_input[i].first->_isMinBound?"+":"-")<< src->_input[i].second << "*" << src->_input[i].first->orig_value() << "\t(" << src->_input[i].first->_met->getName() << ")" << endl;
			}
		}
		if(met->getPotUb() -EPSILON > ub->orig_value()) {
			cout << met->getName() << ": original="<< met->getPotUb() << " induced=" << ub->orig_value() << endl;
			ArcPtr src = ub->_last_update.lock();
			cout << "sourced by: " << src->_creator->getName() << endl;
			for(unsigned int i = 0; i < src->_input.size(); i++) {
				cout <<  "\t" << (src->_input[i].first->_isMinBound?"-":"+")<< src->_input[i].second << "*" << src->_input[i].first->orig_value() << "\t(" << src->_input[i].first->_met->getName() << ")" << endl;
			}
		}
	}
}

shared_ptr<vector<ReactionPtr> > PotBoundPropagation2::getBlockedReactions() {
	shared_ptr<vector<ReactionPtr> > result(new vector<ReactionPtr>());
	foreach(ReactionPtr r, _model->getReactions()) {
		if(!r->isExchange()) {
			double val = 0;
			foreach(Stoichiometry s, r->getProducts()) {
				val += s.second * _minBounds.at(s.first)->orig_value();
			}
			foreach(Stoichiometry s, r->getReactants()) {
				val -= s.second * _maxBounds.at(s.first)->orig_value(); // coefficients are all positive
			}
			cout << r->getName() << ": " << val << endl;
			if(val > -EPSILON) {
				cout << "blocked reaction: " << r->getName() << " by " << val << endl;
				result->push_back(r);
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
	ind[0] = _vars.at(a->_target);
	coef[0] = 1;
	for(unsigned int j = 0; j < a->_input.size(); j++) {
		MetBoundPtr met = a->_input[j].first;
		ind[j+1] = _vars.at(met);
		coef[j+1] = -a->_input[j].second;
	}
	BOOST_SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
}

void PotBoundPropagation2::updateStepHard(ScipModelPtr scip) {
	SCIP_LPI* lpi;
	BOOST_SCIP_CALL( SCIPlpiCreate(&lpi, "potboundprop_updateStepHard", SCIP_OBJSEN_MAXIMIZE) );

	// create constraints of fixed reaction directions
	shared_ptr<unordered_set<ReactionPtr> > dirs = scip->getFixedDirections();
	foreach(ArcPtr a, _arcs) {
		ReactionPtr r = a->_creator;
		if(dirs->find(r) != dirs->end()) {
			if(scip->getCurrentFluxLb(r) < -EPSILON == a->_fwdcreator) {
				// reaction is fixed to reverse
				// for a = (v,A) build a constraint of the form
				// 0 \leq \mu(v) - sum_{x \in A} S_{xa}/S_{va} \mu(x)
				buildHardArcConstraint(a, lpi, false);
				// also do this on the reversed arc, because that involves other bound-variables
				ArcPtr rev = getReversed(a);
				buildHardArcConstraint(rev, lpi, true);
			}
			else {
				// reaction is fixed to forward
				// for a = (v,A) build a constraint of the form
				// \mu(v) - sum_{x \in A} S_{xa}/S_{va} \mu(x) \leq 0
				buildHardArcConstraint(a, lpi, true);
				// also do this on the reversed arc, because that involves other bound-variables
				ArcPtr rev = getReversed(a);
				buildHardArcConstraint(rev, lpi, false);
			}
		}
	}

	// now set bounds on var entries
	foreach(VarEntry v, _vars) {
		double lb = -INFINITY;
		double ub = v.first->_bound;
		double obj = 1;
		int ind = v.second;
		BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
	}

	// solve
	BOOST_SCIP_CALL( SCIPlpiSolveDual(lpi) );
	assert( SCIPlpiWasSolved(lpi) );
	if(SCIPlpiIsPrimalFeasible(lpi)) {
		// we have a feasible solution
		double objval;
		double primsol[_vars.size()];
		BOOST_SCIP_CALL( SCIPlpiGetSol(lpi,&objval, primsol, NULL, NULL, NULL) );
		// perform update
		foreach(VarEntry e, _vars) {
			if(primsol[e.second] < e.first->_bound - EPSILON) {
				e.first->_bound = primsol[e.second];
			}
		}
	}
	else {

	}
}

bool PotBoundPropagation2::updateStepFlow(ScipModelPtr scip) {
	// compute X set
	unordered_set<MetBoundPtr> X;
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		X.insert(_maxBounds[met]);
		X.insert(_minBounds[met]);
	}
	foreach(ArcPtr a, _arcs) {
		double res = 0;
		for(unsigned int i = 0; i < a->_input.size(); i++) {
			res += a->_input[i].second*a->_input[i].first->_bound;
		}
		if(res > a->_target->_bound) {
			X.erase(a->_target);
		}
	}
	// X is now computed
	// we now have to create the LP for solving the update step

	SCIP_LPI* lpi;
	BOOST_SCIP_CALL( SCIPlpiCreate(&lpi, "potboundprop_updateStepFlow", SCIP_OBJSEN_MINIMIZE) );

	// create constraints for arcs
	foreach(ArcPtr a, _arcs) {
		// for a = (v,A) build a constraint of the form
		// \mu(v) - sum_{x \in A \cap X} S_{xa}/S_{va} \mu(x) <= sum_{x \in A \setminus X} S_{xa}/S_{va} _bound(x)
		double lhs = -INFINITY;
		double rhs = 0; // start with zero and increase with metabolites not in X in constraint building process
		int beg = 0;
		vector<int> ind;
		vector<double> coef;
		ind[0] = _vars.at(a->_target);
		coef[0] = 1;
		int k = 1;
		for(unsigned int j = 0; j < a->_input.size(); j++) {
			MetBoundPtr met = a->_input[j].first;
			if(X.find(met) != X.end()) {
				ind[k] = _vars.at(met);
				coef[k] = -a->_input[j].second;
				k++;
			}
			else {
				rhs += a->_input[j].second * met->_bound;
			}
		}
		BOOST_SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, ind.size(), &beg, ind.data(), coef.data()) );
	}
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		double lb = -INFINITY;
		double ub = INFINITY;
		double obj = 1;
		int ind = _vars.at(_maxBounds.at(met));
		BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
		ind = _vars.at(_minBounds.at(met));
		BOOST_SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );
		BOOST_SCIP_CALL( SCIPlpiChgObj(lpi, 1, &ind, &obj) );
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
		bool updated = true;
		foreach(VarEntry e, _vars) {
			if(ray[e.second] < - EPSILON) { // every ray will only have nonpositive entries
				e.first->_bound = -INFINITY;
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
				e.first->_bound = primsol[e.second];
				updated = true;
			}
		}
		return updated;
	}
}

#if 1
void PotBoundPropagation2::update(ScipModelPtr scip) {
	// for now: always compute from scratch
	foreach(MetabolitePtr m, _model->getMetabolites()) {
		_maxBounds[m]->_bound = scip->getCurrentPotentialUb(m);
		_minBounds[m]->_bound = -scip->getCurrentPotentialLb(m);
	}
	bool updated = updateStepHard(scip);
	while(updated) {
		updated = updateStepFlow(scip); // updateStepFlow may not give best update and may require repetitive calls
		if(updated) updateStepHard(scip); // updateStepHard always gives best update
	}
}
#endif

} /* namespace metaopt */
