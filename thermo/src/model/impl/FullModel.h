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
 * FullModel.h
 *
 *  Created on: Mar 28, 2012
 *      Author: arne
 */

#ifndef FULLMODEL_H_
#define FULLMODEL_H_

#include "Properties.h"
#include "model/Model.h"
#include "Uncopyable.h"

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

typedef boost::shared_ptr<FullModel> FullModelPtr;

} /* namespace metaopt */
#endif /* FULLMODEL_H_ */
