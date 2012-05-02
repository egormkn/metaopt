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
	weak_ptr<Model> ptr = shared_from_this();
	ReactionPtr r = shared_ptr<Reaction>(new Reaction(ptr, name));
	_reactions.insert(r);
	updateReaction(r);
	return r;
}

void FullModel::deleteReaction(ReactionPtr reaction) {
	assert(_reactions.find(reaction) != _reactions.end());
	_reactions.erase(reaction);
	_fluxforcing.erase(reaction);
	_objective.erase(reaction);
	_internal.erase(reaction);
}

void FullModel::deleteMetabolite(MetabolitePtr metabolite) {
	assert(_metabolites.find(metabolite) != _metabolites.end());
	assert(metabolite->isIsolated());
}


void FullModel::updateReaction(ReactionPtr reaction) {
	Model::updateReaction(reaction);
}

const MetabolitePtr FullModel::createMetabolite(string name) {
	weak_ptr<Model> ptr = shared_from_this();
	MetabolitePtr m = shared_ptr<Metabolite>(new Metabolite(ptr, name));
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
