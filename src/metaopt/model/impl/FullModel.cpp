/*
 * FullModel.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: arne
 */

#include "FullModel.h"

using namespace boost;

namespace metaopt {

FullModel::FullModel(shared_ptr<libsbml::Model> model) {
	_model = model;
}

FullModel::~FullModel() {
	// nothing to do
}

const ReactionPtr FullModel::createReaction() {
	weak_ptr<Model> ptr = shared_from_this();
	ReactionPtr r = shared_ptr<Reaction>(new Reaction(ptr));
	_reactions.insert(r);
	return r;
}

void FullModel::deleteReaction(ReactionPtr reaction) {
	assert(_reactions.find(reaction) != _reactions.end());
	_reactions.erase(reaction);
	_fluxforcing.erase(reaction);
	_objective.erase(reaction);
}

void FullModel::deleteMetabolite(MetabolitePtr metabolite) {
	assert(_metabolites.find(metabolite) != _metabolites.end());
	assert(metabolite->isIsolated());
}

const MetabolitePtr FullModel::createMetabolite() {
	return MetabolitePtr();
}

const unordered_set<ReactionPtr>& FullModel::getReactions() const {
	return _reactions;
}

const unordered_set<MetabolitePtr>& FullModel::getMetabolites() const {
	return _metabolites;
}




} /* namespace metaopt */
