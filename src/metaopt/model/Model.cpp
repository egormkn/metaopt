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

void Model::deleteReaction(ReactionPtr reaction) {
	assert(_reactions.find(reaction) != _reactions.end());
	_reactions.erase(reaction);
	_fluxforcing.erase(reaction);
	_objective.erase(reaction);
}

void Model::deleteMetabolite(MetabolitePtr metabolite) {
	assert(_metabolites.find(metabolite) != _metabolites.end());
	assert(metabolite->isIsolated());
}

const unordered_set<ReactionPtr>& Model::getFluxForcingReactions() const {
	return _fluxforcing;
}

const unordered_set<ReactionPtr>& Model::getObjectiveReactions() const {
	return _objective;
}

const unordered_set<ReactionPtr>& Model::getReactions() const {
	return _reactions;
}

const unordered_set<MetabolitePtr>& Model::getMetabolites() const {
	return _metabolites;
}

void Model::updateReaction(ReactionPtr reaction) {
	assert(_reactions.find(reaction) != _reactions.end());
	bool was_fluxforcing = _fluxforcing.find(reaction) != _fluxforcing.end();
	bool was_objective = _objective.find(reaction) != _objective.end();
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
}

void Model::updateMetabolite(MetabolitePtr metabolite) {
	assert(_metabolites.find(metabolite) != _metabolites.end());
	// nothing to do yet
}

} /* namespace metaopt */
