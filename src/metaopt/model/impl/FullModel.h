/*
 * FullModel.h
 *
 *  Created on: Mar 28, 2012
 *      Author: arne
 */

#ifndef FULLMODEL_H_
#define FULLMODEL_H_

#include "metaopt/Properties.h"
#include "metaopt/model/Model.h"

namespace metaopt {

/**
 * This class represents fully initialized models
 */
class FullModel: public metaopt::Model {
public:
	FullModel();
	virtual ~FullModel();

	//////////////////////////////////////////////////
	// implemented methods from Model
	//////////////////////////////////////////////////

	/**
	 * Returns the list of reactions in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	const boost::unordered_set<ReactionPtr>& getReactions() const;

	/**
	 * Returns the list of metabolites in this model.
	 *
	 * Important: This returns a reference to the internally stored list. Do not modify it.
	 */
	const boost::unordered_set<MetabolitePtr>& getMetabolites() const;

	/**
	 * Creates a new Reaction in this model. The created Reaction is automatically added to this Model.
	 *
	 * @returns the created Reaction for setting bounds, objective coefficient and stoichiometries.
	 */
	const ReactionPtr createReaction(std::string name);

	/**
	 * Deletes a Reaction from this Model.
	 *
	 * @param reaction the Reaction to delete
	 */
	void deleteReaction(ReactionPtr reaction);

	/**
	 * Updates a reaction in this Model.
	 *
	 * This method must be called, if bounds or objective coefficients have changed.
	 * In these cases it may happen that the reaction becomes flux forcing etc. and then the corresponding lists have to be updated.
	 *
	 * @param reaction the Reaction to update
	 */
	void updateReaction(ReactionPtr reaction);

	/**
	 * Creates a new Metabolite in this model. The created Metabolite is automatically added to this Model.
	 *
	 * @returns the created Metabolite for setting bounds, stoichiometries, etc.
	 */
	const MetabolitePtr createMetabolite(std::string name);

	/**
	 * Deletes a Metabolite from this Model.
	 *
	 * @param metabolite the Metabolite to delete
	 */
	void deleteMetabolite(MetabolitePtr metabolite);

	/**
	 * Updates a metabolite in this Model.
	 *
	 * This method must be called, if bounds etc. have changed.
	 *
	 * @param metabolite the Metabolite to update
	 */
	void updateMetabolite(MetabolitePtr metabolite);

private:
	boost::unordered_set<ReactionPtr> _reactions; /** list of all reactions */
	boost::unordered_set<MetabolitePtr> _metabolites; /** list of all metabolites */
};

} /* namespace metaopt */
#endif /* FULLMODEL_H_ */
