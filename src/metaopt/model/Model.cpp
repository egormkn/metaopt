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
 * Model.cpp
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#include <vector>
#include <algorithm>

#include "Precision.h"
#include "../Properties.h"
#include "Model.h"

using namespace std;
using namespace boost;

namespace metaopt {

#define DEFAULT_FLUX_PRECISION 1e-6
#define DEFAULT_POT_PRECISION 1e-6
#define DEFAULT_COEF_PRECISION 1e-10

Model::Model() :
	_flux_precision(new Precision(DEFAULT_FLUX_PRECISION)),
	_pot_precision(new Precision(DEFAULT_POT_PRECISION)),
	_coef_precision(new Precision(DEFAULT_COEF_PRECISION))
{
	// nothing else to do
}

Model::~Model() {
	// nothing to do yet
}

const unordered_set<ReactionPtr>& Model::getFluxForcingReactions() const {
	return _fluxforcing;
}

const unordered_set<ReactionPtr>& Model::getProblematicReactions() const {
	return _problematic;
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
	bool was_problematic = _problematic.find(reaction) != _problematic.end();
	bool was_objective = _objective.find(reaction) != _objective.end();
	bool was_internal = _internal.find(reaction) != _internal.end();
	if(was_fluxforcing) {
		if(!reaction->isFluxForcing()) _fluxforcing.erase(reaction);
	}
	else {
		if(reaction->isFluxForcing()) _fluxforcing.insert(reaction);
	}
	if(was_problematic) {
		if(!reaction->isProblematic()) _problematic.erase(reaction);
	}
	else {
		if(reaction->isProblematic()) _problematic.insert(reaction);
	}
	if(was_objective) {
		if(!reaction->isObjective()) _objective.erase(reaction);
	}
	else {
		if(reaction->isObjective()) _objective.insert(reaction);
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
		if(rxn->getName().compare(name) == 0) return rxn;
	}
	assert(false); // Reaction does not exist!
	return(ReactionPtr());
}

bool Model::hasMetabolite(string name) {
	foreach(MetabolitePtr met, getMetabolites()) {
		if(met->getName() == name) return true;
	}
	return false;
}

bool Model::hasReaction(string name) {
	foreach(ReactionPtr rxn, getReactions()) {
		if(rxn->getName() == name) return true;
	}
	return false;
}

void Model::setFluxPrecision(PrecisionPtr precision) {
	// directly update reactions and update lists of flux forcing, objective reactions.
	_flux_precision = precision;
	foreach(ReactionPtr rxn, getReactions()) {
		rxn->setFluxPrecision(precision);
		updateReaction(rxn);
	}
}

void Model::setPotPrecision(PrecisionPtr precision) {
	_pot_precision = precision;
	// also update prcision values of metabolites
	foreach(MetabolitePtr met, getMetabolites()) {
		met->setPotPrecision(precision);
	}
}

void Model::setCoefPrecision(PrecisionPtr precision) {
	_coef_precision = precision;
}


} /* namespace metaopt */
