/*
 * Model.cpp
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#include <vector>
#include <algorithm>

#include "../Properties.h"
#include "Model.h"

using namespace std;
using namespace boost;

namespace metaopt {

Model::Model() {
	// nothing to do yet
}

Model::~Model() {
	// nothing to do yet
}

const unordered_set<ReactionPtr>& Model::getFluxForcingReactions() const {
	return _fluxforcing;
}

const unordered_set<ReactionPtr>& Model::getObjectiveReactions() const {
	return _objective;
}

const unordered_set<ReactionPtr>& Model::getInternalReactions() const {
	return _internal;
}


void Model::updateReaction(ReactionPtr reaction) {
	assert(getReactions().find(reaction) != getReactions().end());
	bool was_fluxforcing = _fluxforcing.find(reaction) != _fluxforcing.end();
	bool was_objective = _objective.find(reaction) != _objective.end();
	bool was_internal = _internal.find(reaction) != _internal.end();
	if(was_fluxforcing) {
		if(reaction->getLb() < EPSILON && reaction->getUb() > -EPSILON) _fluxforcing.erase(reaction);
	}
	else {
		if(reaction->getLb() > EPSILON || reaction->getUb() < -EPSILON) _fluxforcing.insert(reaction);
	}
	if(was_objective) {
		if(reaction->getObj() < EPSILON && reaction->getObj() > -EPSILON) _objective.erase(reaction);
	}
	else {
		if(reaction->getObj() > EPSILON || reaction->getObj() < -EPSILON) _objective.insert(reaction);
	}
	if(was_internal) {
		if(reaction->isExchange()) _internal.erase(reaction);
	}
	else {
		if(!reaction->isExchange()) _internal.insert(reaction);
	}
}

void Model::updateMetabolite(MetabolitePtr metabolite) {
	assert(getMetabolites().find(metabolite) != getMetabolites().end());
	// nothing to do yet
}

MetabolitePtr Model::getMetabolite(string name) {
	foreach(MetabolitePtr met, getMetabolites()) {
		if(met->getName() == name) return met;
	}
	assert(false);
	return(MetabolitePtr());
}

ReactionPtr Model::getReaction(string name) {
	foreach(ReactionPtr rxn, getReactions()) {
		if(rxn->getName() == name) return rxn;
	}
	assert(false); // Reaction does not exist!
	return(ReactionPtr());
}

} /* namespace metaopt */
