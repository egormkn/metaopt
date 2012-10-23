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
