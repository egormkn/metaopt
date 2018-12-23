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
 * FullModel.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: arne
 */

#include "FullModel.h"

using namespace std;
using namespace boost;

namespace metaopt {

FullModel::FullModel() {
}

FullModel::~FullModel() {
	// nothing to do
}

const ReactionPtr FullModel::createReaction(string name) {
	boost::weak_ptr<Model> ptr = shared_from_this();
	ReactionPtr r = boost::shared_ptr<Reaction>(new Reaction(ptr, name));
	_reactions.insert(r);
	updateReaction(r);
	return r;
}

void FullModel::deleteReaction(ReactionPtr reaction) {
	assert(_reactions.find(reaction) != _reactions.end());
	// first clear the reaction by removing all reference to metabolites
	vector<Stoichiometry> mets(reaction->getStoichiometries().begin(), reaction->getStoichiometries().end()); // copy list to get rid of concurrency issues
	foreach(Stoichiometry s, mets) {
		reaction->setStoichiometry(s.first,0);
	}
	_reactions.erase(reaction);
	_fluxforcing.erase(reaction);
	_objective.erase(reaction);
	_internal.erase(reaction);
}

void FullModel::deleteMetabolite(MetabolitePtr metabolite) {
	assert(_metabolites.find(metabolite) != _metabolites.end());
	assert(metabolite->isIsolated());
	_metabolites.erase(metabolite);
}


void FullModel::updateReaction(ReactionPtr reaction) {
	Model::updateReaction(reaction);
}

const MetabolitePtr FullModel::createMetabolite(string name) {
	boost::weak_ptr<Model> ptr = shared_from_this();
	MetabolitePtr m = boost::shared_ptr<Metabolite>(new Metabolite(ptr, name));
	_metabolites.insert(m);
	updateMetabolite(m);
	return m;
}

void FullModel::updateMetabolite(MetabolitePtr metabolite) {
	Model::updateMetabolite(metabolite);
}

const unordered_set<ReactionPtr>& FullModel::getReactions() const {
	return _reactions;
}

const unordered_set<MetabolitePtr>& FullModel::getMetabolites() const {
	return _metabolites;
}




} /* namespace metaopt */
