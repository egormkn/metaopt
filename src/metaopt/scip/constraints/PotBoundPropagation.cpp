/*
 * PotBoundPropagation.cpp
 *
 *  Created on: 18.05.2012
 *      Author: arnem
 */

#include <boost/unordered_set.hpp>
#include <iostream>
#include "PotBoundPropagation.h"
#include "metaopt/model/Reaction.h"
#include "metaopt/Properties.h"

using namespace boost;
using namespace std;

namespace metaopt {

PotBoundPropagation::PotBoundPropagation(ModelPtr model) :
	_model(model)
{
	init_Queue();
	while(!_queue.empty()) {
		MetBoundPtr m = pop();
		update(m);
	}
}

PotBoundPropagation::~PotBoundPropagation() {
	// TODO Auto-generated destructor stub
}

PotBoundPropagation::MetBound::MetBound(MetabolitePtr met, bool isMinBound) :
		_met(met), _isMinBound(isMinBound)
{
	if(isMinBound) {
		_bound = INFINITY;
	}
	else {
		_bound = -INFINITY;
	}
}

double PotBoundPropagation::MetBound::restrictivity() const {
	if(_isMinBound) {
		return _bound - _met->getPotLb();
	}
	else {
		return _met->getPotUb() - _bound;
	}
}

void PotBoundPropagation::MetBound::updateUnused() {
	_unused = shared_ptr<unordered_set<MetabolitePtr> >(new unordered_set<MetabolitePtr>());
	typedef std::pair<ReactionPtr, shared_ptr< unordered_set<MetabolitePtr> > > UnusedFromEntry;

	foreach(UnusedFromEntry e, _unusedFrom) {
		_unused->insert(e.second->begin(), e.second->end());
	}
}

void PotBoundPropagation::init_Queue() {
	foreach(MetabolitePtr m, _model->getMetabolites()) {
		MetBoundPtr max(new MetBound(m, false));
		max->_unused = shared_ptr<unordered_set<MetabolitePtr> >(new unordered_set<MetabolitePtr>( _model->getMetabolites() ));
		_maxBounds[m] = max;
		MetBoundPtr min(new MetBound(m, true));
		min->_unused = shared_ptr<unordered_set<MetabolitePtr> >(new unordered_set<MetabolitePtr>( _model->getMetabolites() ));
		_minBounds[m] = min;
	}

	foreach(ReactionPtr r, _model->getReactions()) {
		if(r->isExchange()) {
			foreach(Stoichiometry s, r->getProducts()) {
				MetBoundPtr met = _maxBounds.at(s.first);
				met->_bound = s.first->getPotUb();
				enqueue(met);
			}
			foreach(Stoichiometry s, r->getReactants()) {
				MetBoundPtr met = _minBounds.at(s.first);
				met->_bound = s.first->getPotLb();
				enqueue(met);
			}
		}
	}
}

void PotBoundPropagation::updateProducts(ReactionPtr& r, MetabolitePtr& ignore, shared_ptr<unordered_set<MetabolitePtr> >& unused, double& val) {
	foreach(Stoichiometry s, r->getProducts()) {
		MetabolitePtr m = s.first;
		if(m != ignore) {
			MetBoundPtr mb = _minBounds[m];
			for(unordered_set<MetabolitePtr>::iterator iter = unused->begin(); iter != unused->end();) {
				if(mb->_unused->find(*iter) == mb->_unused->end()) {
					iter = unused->erase(iter);
				}
				else {
					iter++;
				}
			}
			if(isinf(mb->_bound)) {
				// use lb as pessimistic guess
				val -= r->getStoichiometry(m) * m->getPotLb();
			}
			else {
				val -= r->getStoichiometry(m) * mb->_bound;
			}
		}
	}
}

void PotBoundPropagation::updateReactants(ReactionPtr& r, MetabolitePtr& ignore, shared_ptr<unordered_set<MetabolitePtr> >& unused, double& val) {
	foreach(Stoichiometry s, r->getReactants()) {
		MetabolitePtr m = s.first;
		if(m != ignore) {
			MetBoundPtr mb = _maxBounds[m];
			for(unordered_set<MetabolitePtr>::iterator iter = unused->begin(); iter != unused->end();) {
				if(mb->_unused->find(*iter) == mb->_unused->end()) {
					iter = unused->erase(iter);
				}
				else {
					iter++;
				}
			}
			if(isinf(mb->_bound)) {
				// use ub as pessimistic guess
				val -= r->getStoichiometry(m) * m->getPotUb();
			}
			else {
				val -= r->getStoichiometry(m) * mb->_bound;
			}
		}
	}
}

template<class E> bool isSubset(unordered_set<E>& subset, unordered_set<E>& set) {
	foreach(E e, subset) {
		if(set.find(e) == set.end()) return false;
	}
	return true;
}

void PotBoundPropagation::update(MetBoundPtr b) {
	// TODO: I am not 100% sure, if what I am doing is even theoretically correct!!

	if(b->_isMinBound) {
		shared_ptr<vector<ReactionPtr> > producers = b->_met->getProducers();
		foreach(ReactionPtr r, *producers) {
			// compute intersection of all unused metabolites
			// and compute new value
			// for each reactant:
			// only if the new value is less than the old value and the reactant is unused, do an update
			shared_ptr<unordered_set<MetabolitePtr> > unused(b->_unused);
			double val = - r->getStoichiometry(b->_met) * b->_bound;
			// first compute val and unused only for products (because thats the same for each reactant)
			updateProducts(r, b->_met, unused, val);

			// now do the whole thing for each reactant that we update
			foreach(Stoichiometry s, r->getReactants()) {
				MetabolitePtr x = s.first;
				// we want to update x
				if(unused->find(x) != unused->end()) {
					shared_ptr<unordered_set<MetabolitePtr> > unusedComplete(unused);
					double valComplete = val;
					updateReactants(r, x, unusedComplete, valComplete);
					valComplete /= r->getStoichiometry(x);
					MetBoundPtr xb = _minBounds.at(x);
					if(unusedComplete->find(x) != unused->end()) {
						// x is still unused

						// remove x and all other reactants to obtain final list
						foreach(Stoichiometry s, r->getReactants()) {
							unusedComplete->erase(s.first);
						}

						// check if updates are necessary
						if( valComplete < xb->_bound && xb->_bound > x->getPotLb()+EPSILON) {
							// and bound is also less, so we have to update
							cout << "[" << x->getName() << "] update min bound: " << valComplete << " was: " << xb->_bound << " of total: " << x->getPotLb() << endl;
							xb->_bound = max(valComplete, x->getPotLb());
							enqueue(xb);
						}
						if( !isSubset(*unusedComplete, *(xb->_unused)) ) {
							cout << "[" << x->getName() << "] increased unused (outflow)" << endl;
							enqueue(xb);
						}

						// and we always should update the set of unused values
						xb->_unusedFrom[r] = unusedComplete;
						xb->updateUnused();
					}
				}
			}
		}
	}
	else {
		// same thing for maxBound
		shared_ptr<vector<ReactionPtr> > consumers = b->_met->getConsumers();
		foreach(ReactionPtr r, *consumers) {
			// compute intersection of all unused metabolites
			// and compute new value
			// for each reactant:
			// only if the new value is less than the old value and the reactant is unused, do an update
			shared_ptr<unordered_set<MetabolitePtr> > unused(b->_unused);
			double val = - r->getStoichiometry(b->_met) * b->_bound;
			// first compute val and unused only for products (because thats the same for each reactant)
			updateReactants(r, b->_met, unused, val);

			// now do the whole thing for each reactant that we update
			foreach(Stoichiometry s, r->getProducts()) {
				MetabolitePtr x = s.first;
				// we want to update x
				if(unused->find(x) != unused->end()) {
					shared_ptr<unordered_set<MetabolitePtr> > unusedComplete(unused);
					double valComplete = val;
					updateProducts(r, x, unusedComplete, valComplete);
					valComplete /= r->getStoichiometry(x);
					MetBoundPtr xb = _maxBounds.at(x);
					if(unusedComplete->find(x) != unused->end()) {
						// x is still unused

						// remove x and all other reactants to obtain final list
						foreach(Stoichiometry s, r->getProducts()) {
							unusedComplete->erase(s.first);
						}

						// check if updates are necessary
						if( valComplete > xb->_bound && xb->_bound < x->getPotUb()-EPSILON) {
							// and bound is also less, so we have to update
							cout << "[" << x->getName() << "] update max bound: " << valComplete << " was: " << xb->_bound << " of total: " << x->getPotUb()<< endl;
							xb->_bound = min(valComplete, x->getPotUb());
							enqueue(xb);
						}
						if( !isSubset(*unusedComplete, *(xb->_unused)) ) {
							cout << "[" << x->getName() << "] increased unused (inflow)" << endl;
							enqueue(xb);
						}

						// and we always should update the set of unused values
						xb->_unusedFrom[r] = unusedComplete;
						xb->updateUnused();
					}
				}
			}
		}
	}
}

void PotBoundPropagation::print() {
	foreach(MetabolitePtr met, _model->getMetabolites()) {
		MetBoundPtr lb = _minBounds[met];
		MetBoundPtr ub = _maxBounds[met];
		cout << met->getPotLb() << " " << lb->_bound << " " << ub->_bound << " " << met->getPotUb() << endl;
	}
}

} /* namespace metaopt */
