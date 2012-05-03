/*
 * ModelAddOn.cpp
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#include "ModelAddOn.h"

using namespace std;

namespace metaopt {

ModelAddOn::ModelAddOn(ScipModelPtr model, string name) {
	_model = model;
	_name = name;
	_destroyed = false;
}

ModelAddOn::~ModelAddOn() {
	assert(_destroyed);
	// nothing to do
}

void ModelAddOn::destroy(ScipModel* model) {
	_destroyed = true;
}

void ModelAddOn::setDirection(SCIP_NODE* node, ReactionPtr rxn, bool fwd) {
	// don't do anything by default.
}

boost::shared_ptr<const boost::unordered_set<ReactionPtr> > ModelAddOn::getFixedDirections() {
	boost::shared_ptr<const boost::unordered_set<ReactionPtr> > result(new boost::unordered_set<ReactionPtr>());
	return result;
}

} /* namespace metaopt */
