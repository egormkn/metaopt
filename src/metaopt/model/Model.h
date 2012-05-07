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
	 * gets a metabolite by name.
	 * Attention: Slow implementation!
	 * Only use this for debugging purposes!
	 */
	MetabolitePtr getMetabolite(std::string name);

	/**
	 * gets a reaction by name
	 * Attention: Slow implementation!
	 * Only use this for debugging purposes!
	 */
	ReactionPtr getReaction(std::string name);


protected:
	boost::unordered_set<ReactionPtr> _fluxforcing; /** list of all flux forcing reactions */
	boost::unordered_set<ReactionPtr> _objective; /** list of all objective reactions */
	boost::unordered_set<ReactionPtr> _internal; /** list of all internal reactions */

};

/** Please use ModelPtr to handle pointers to Model instances */
typedef boost::shared_ptr<Model> ModelPtr;


} /* namespace metaopt */
#endif /* MODEL_H_ */
