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

#include "../Properties.h"
#include "Reaction.h"
#include "Metabolite.h"
#include "Model.h"

using namespace std;
using namespace boost;

namespace metaopt {

Metabolite::Metabolite() {
	// TODO Auto-generated constructor stub

}

Metabolite::~Metabolite() {
	// TODO Auto-generated destructor stub
}

inline void Metabolite::notifyChange() {
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


} /* namespace metaopt */
