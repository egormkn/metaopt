/*
 * Reaction.cpp
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#include <boost/throw_exception.hpp>
#include <string>
#include <iostream>

#include "Reaction.h"
#include "Model.h"
#include "../Properties.h"

using namespace boost;
using namespace std;

namespace metaopt {

Reaction::Reaction(weak_ptr<Model> model) {
	_model = model;
}

Reaction::~Reaction() {
	// weak pointer auto-destructs
}

inline void Reaction::notifyChange() {
	if(_model.expired()) {
		BOOST_THROW_EXCEPTION( ModelOwnershipError() << reaction_name(getName()) );
	}
	_model.lock()->updateReaction(shared_from_this());
}

const unordered_map<MetabolitePtr, double> Reaction::getStoichiometries() const {
	return _stoichiometries;
}

const unordered_map<MetabolitePtr, double> Reaction::getReactants() const {
	return _reactants;
}

const unordered_map<MetabolitePtr, double> Reaction::getProducts() const {
	return _products;
}

double Reaction::getStoichiometry(MetabolitePtr m) const {
	assert(m->getOwner() == getOwner());
	unordered_map<MetabolitePtr, double>::const_iterator i = _stoichiometries.find(m);
	if(i == _stoichiometries.end()) return 0;
	else return i->second;
}

double Reaction::getReactant(MetabolitePtr m) const {
	assert(m->getOwner() == getOwner());
	unordered_map<MetabolitePtr, double>::const_iterator i = _reactants.find(m);
	if(i == _reactants.end()) return 0;
	else return i->second;
}

double Reaction::getProduct(MetabolitePtr m) const {
	assert(m->getOwner() == getOwner());
	unordered_map<MetabolitePtr, double>::const_iterator i = _products.find(m);
	if(i == _products.end()) return 0;
	else return i->second;
}

void Reaction::setStoichiometry(MetabolitePtr m, double value) {
	assert(m->getOwner() == getOwner());
	if(value > -EPSILON && value < EPSILON) {
		// value == 0, so delete entry
		_stoichiometries.erase(m);
		_reactants.erase(m);
		_products.erase(m);
	}
	else {
		_stoichiometries[m] = value;
		if(value > 0) {
			_products[m] = value;
			_reactants.erase(m);
		}
		else { // value < 0
			_products.erase(m);
			_reactants[m] = -value;
		}
	}
	assert(getStoichiometry(m) == getProduct(m) - getReactant(m));
}

void Reaction::setReactant(MetabolitePtr m, double value) {
	assert(m->getOwner() == getOwner());
	if(value > -EPSILON && value < EPSILON) {
		// value == 0, so delete entry
		_reactants.erase(m);
		if(_products.find(m) == _products.end()) {
			_stoichiometries.erase(m);
		}
	}
	else {
		assert(value > 0);
		_reactants[m] = value;
		_stoichiometries[m] = getProduct(m) - value;
	}
	assert(getStoichiometry(m) == getProduct(m) - getReactant(m));
}

void Reaction::setProduct(MetabolitePtr m, double value) {
	assert(m->getOwner() == getOwner());
	if(value > -EPSILON && value < EPSILON) {
		// value == 0, so delete entry
		_products.erase(m);
		if(_reactants.find(m) == _reactants.end()) {
			_stoichiometries.erase(m);
		}
	}
	else {
		assert(value > 0);
		_products[m] = value;
		_stoichiometries[m] = value - getReactant(m);
	}
	assert(getStoichiometry(m) == getProduct(m) - getReactant(m));
}

const boost::shared_ptr<Model> Reaction::getOwner() const {
	if(_model.expired()) {
		BOOST_THROW_EXCEPTION(ModelOwnershipError() << reaction_name(getName()));
	}
	return _model.lock();
}


} /* namespace metaopt */
