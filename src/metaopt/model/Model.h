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
 * Model.h
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>
#include "Reaction.h"
#include "Metabolite.h"
#include "metaopt/Uncopyable.h"
#include "metaopt/Properties.h"

namespace metaopt {

struct ModelOwnershipError : virtual boost::exception, virtual std::exception {
public:
	ModelOwnershipError() {
		int foo = 0;
	}
};


/** This class represents metabolic network models.
 *
 * A metabolic network consists of a list of reactions and metabolites.
 * Depending on the properties of the reactions, some reactions are flux forcing, respectively contain non-zero objective coefficient.
 *
 * A Model stores a list of all reactions. A Reaction must not exist without a Model.
 * For performance reasons, Model also maintains lists of flux forcing reactions, i.e. reactions that have positive lower bound or negative upper bound.
 * Also, a list of reactions with objective coefficient unequal to zero is maintained for the same reason.
 * Hence, it is important to call updateReaction(ReactionPtr) if a reaction was modified (usually the reaction implementation does this automatically).
 */
class Model  : public boost::enable_shared_from_this<Model>, Uncopyable {
public:
	Model();
	virtual ~Model();

	/**
	 * Returns the list of reactions in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	virtual const boost::unordered_set<ReactionPtr>& getReactions() const = 0;

	/**
	 * Returns the list of flux forcing reactions in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	const boost::unordered_set<ReactionPtr>& getFluxForcingReactions() const;

	/**
	 * Returns the list of problematic reactions in this model.
	 * A reaction is problematic, if it is involved in additional constraints to steady-state and thermodynamics.
	 * Problematic reactions have similar effects on the constraint handler and heuristic as flux forcing reactions.
	 * However, we distinguish them explicitly from flux forcing reactions, because flux forcing is an implicit property,
	 * where as problematic reactions have to be marked as such explicitly by the user.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	const boost::unordered_set<ReactionPtr>& getProblematicReactions() const;

	/**
	 * Returns the list of objective reactions in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	const boost::unordered_set<ReactionPtr>& getObjectiveReactions() const;

	/**
	 * Returns the list of all internal reactions in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	const boost::unordered_set<ReactionPtr>& getInternalReactions() const;

	/**
	 * Updates a reaction in this Model.
	 *
	 * This method must be called, if bounds or objective coefficients have changed.
	 * In these cases it may happen that the reaction becomes flux forcing etc. and then the corresponding lists have to be updated.
	 *
	 * @param reaction the Reaction to update
	 */
	virtual void updateReaction(ReactionPtr reaction);

	/**
	 * Returns the list of metabolites in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	virtual const boost::unordered_set<MetabolitePtr>& getMetabolites() const = 0;

	/**
	 * Updates a metabolite in this Model.
	 *
	 * This method must be called, if bounds etc. have changed.
	 *
	 * @param metabolite the Metabolite to update
	 */
	virtual void updateMetabolite(MetabolitePtr metabolite);

	/**
	 * Tests, if there exists a metabolite with this name in this model
	 */
	bool hasMetabolite(std::string name);

	/**
	 * Tests, if there exists a reaction with this name in this model
	 */
	bool hasReaction(std::string name);

	/**
	 * gets a metabolite by name.
	 * Attention: Slow implementation!
	 * Attention: Names need not be unique!
	 */
	MetabolitePtr getMetabolite(std::string name);

	/**
	 * gets a reaction by name
	 * Attention: Slow implementation!
	 * Attention: Names need not be unique!
	 */
	ReactionPtr getReaction(std::string name);


protected:
	boost::unordered_set<ReactionPtr> _problematic; /** list of all reactions that are marked as problematic */
	boost::unordered_set<ReactionPtr> _fluxforcing; /** list of all flux forcing reactions */
	boost::unordered_set<ReactionPtr> _objective; /** list of all objective reactions */
	boost::unordered_set<ReactionPtr> _internal; /** list of all internal reactions */

};

/** Please use ModelPtr to handle pointers to Model instances */
typedef boost::shared_ptr<Model> ModelPtr;


} /* namespace metaopt */
#endif /* MODEL_H_ */
