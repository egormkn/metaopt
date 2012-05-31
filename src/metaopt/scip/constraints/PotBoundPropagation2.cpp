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
		foreach(Stoichiometry v, r->getStoichiometries()) {
			ArcPtr a(new Arc());
			a->_creator = r;
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

	foreach(ReactionPtr r, _model->getReactions()) {
		if(r->isExchange()) {
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

} /* namespace metaopt */
