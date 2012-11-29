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
 * Metabolite.cpp
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#include <string>
#include <boost/throw_exception.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

#include <iostream>

#include "../Properties.h"
#include "Reaction.h"
#include "Metabolite.h"
#include "Model.h"

using namespace std;
using namespace boost;

namespace metaopt {

Metabolite::Metabolite(weak_ptr<Model> model, std::string name) {
	_model = model;
	_name = name;
	_lb = -INFINITY;
	_ub = INFINITY;
	_obj = 0;
	_boundary = false;
}

Metabolite::~Metabolite() {
	// TODO Auto-generated destructor stub
}


void Metabolite::notifyChange() {
	if(_model.expired()) {
		BOOST_THROW_EXCEPTION(ModelOwnershipError() << metabolite_name(getName()));
	}
	_model.lock()->updateMetabolite(shared_from_this());
}

bool Metabolite::isIsolated() const {
	if(_model.expired()) {
		BOOST_THROW_EXCEPTION(ModelOwnershipError() << metabolite_name(getName()));
	}
	// check if any of the reactions contains this metabolite in the stoichiometries (which are the union of reactants and products)
	const unordered_set<ReactionPtr> & reactions = _model.lock()->getReactions();
	foreach( ReactionPtr r , reactions )
	{
		typedef const unordered_map<MetabolitePtr, double> map;
		foreach( map::value_type v, r->getStoichiometries()) {
			if(v.first == shared_from_this()) return false;
		}
	}
	return true;
}

const boost::shared_ptr<Model> Metabolite::getOwner() const {
	if(_model.expired()) {
		BOOST_THROW_EXCEPTION(ModelOwnershipError() << metabolite_name(getName()));
	}
	return _model.lock();
}

boost::shared_ptr<std::vector<ReactionPtr> > Metabolite::getProducers() const {
	boost::shared_ptr<std::vector<ReactionPtr> > out(new vector<ReactionPtr>());
	foreach(weak_ptr<Reaction> r, _producers) {
		ReactionPtr sr = r.lock();
		out->push_back(sr);
	}
	return out;
}
boost::shared_ptr<std::vector<ReactionPtr> > Metabolite::getConsumers() const {
	boost::shared_ptr<std::vector<ReactionPtr> > out(new vector<ReactionPtr>());
	foreach(weak_ptr<Reaction> r, _consumers) {
		ReactionPtr sr = r.lock();
		out->push_back(sr);
	}
	return out;
}

struct WeakFind {
	const weak_ptr<Reaction> wr;

	bool operator()(const weak_ptr<Reaction> o) {
		bool out = !(wr < o) && !(o < wr);
		return out;
	}

	WeakFind(ReactionPtr r) : wr(r) {}
};

void Metabolite::removeProducer(ReactionPtr r) {
	//std::cout << _producers.size() << std::endl;
	_producers.erase(remove_if(_producers.begin(), _producers.end(), WeakFind(r)), _producers.end());
	//std::cout << _producers.size() << std::endl;
}
void Metabolite::removeConsumer(ReactionPtr r) {
	//std::cout << _consumers.size() << std::endl;
	_consumers.erase(remove_if(_consumers.begin(), _consumers.end(), WeakFind(r)), _consumers.end());
	//std::cout << _consumers.size() << std::endl;
}
void Metabolite::addProducer(ReactionPtr r) {
	_producers.push_back(r);
}
void Metabolite::addConsumer(ReactionPtr r) {
	_consumers.push_back(r);
}

void Metabolite::setPotPrecision(PrecisionPtr precision) {
	_potPrecision = precision;
}

std::size_t hash_value(MetabolitePtr const & met ) {
	   return reinterpret_cast<std::size_t>(met.get());
}


} /* namespace metaopt */
