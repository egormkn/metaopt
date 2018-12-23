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
 * ReducedScipFluxModel.h
 *
 *  Created on: 04.07.2012
 *      Author: arnem
 */

#ifndef REDUCEDSCIPFLUXMODEL_H_
#define REDUCEDSCIPFLUXMODEL_H_

#include "boost/unordered_map.hpp"

#include "model/scip/AbstractScipFluxModel.h"
#include "model/scip/ScipModel.h"
#include "Properties.h"

namespace metaopt {

class ReducedScipFluxModel : public AbstractScipFluxModel {
public:
	ReducedScipFluxModel(ScipModelPtr origin);
	virtual ~ReducedScipFluxModel();

	/** returns true if objective sense is set to maximize.
	 * If objective sense is minimize, false is returned.
	 */
	virtual bool isMaximize() const;

	/**
	 * fetches the scip variable representing the flux.
	 * If the variable does not exist yet, it is created.
	 * The default bounds are the bounds stored by the reaction.
	 *
	 * The pointer will be valid as long as this ScipModel lives.
	 */
	virtual SCIP_VAR* getFlux(ReactionPtr rxn);

	/**
	 * Associates a (transformed) variable to the reaction
	 */
	void setFlux(ReactionPtr rxn, SCIP_VAR* var);

	/**
	 * fetches the current upper bound of the flux.
	 * If the corresponding variable exists, its value of the current node is returned.
	 * Else, the default bound is returned.
	 */
	virtual double getCurrentFluxUb(ReactionPtr rxn);
	/**
	 * fetches the current lower bound of the flux.
	 * If the corresponding variable exists, its value of the current node is returned.
	 * Else, the default bound is returned.
	 */
	virtual double getCurrentFluxLb(ReactionPtr rxn);

	/**
	 * Checks, if there exists a current flux solution (of the current node in search tree).
	 * This is the case, if the problem has been solved, or if the LP relaxation of the current node has been solved to optimality.
	 * It can also happen that no flux variables have been initialized at all. In this case, this method will always return false.
	 */
	virtual bool hasCurrentFlux();

	/**
	 * Returns the flux value through the specified reaction.
	 * If the problem has been solved, the value of an optimal solution is returned.
	 * Else, the value of the current LP-relaxation is returned. The LP-relaxation, however has to be solved to optimality.
	 *
	 * Precondition: hasFluxVar(rxn),  hasCurrentFlux()
	 *
	 * If the model has no current flux, a PreconditionViolatedException is thrown
	 */
	virtual double getCurrentFlux(ReactionPtr rxn);


	/**
	 * Returns the flux value of specified solution through the specified reaction.
	 */
	virtual double getFlux(SolutionPtr sol, ReactionPtr rxn);

	/**
	 * checks if the variable representing the flux has already been created
	 */
	virtual bool hasFluxVar(ReactionPtr rxn);

	/**
	 * Tests if this model has any flux variables initialized
	 */
	virtual bool hasFlux();

private:
	boost::weak_ptr<ScipModel> _origin;

	inline ScipModelPtr getScip() const;

	boost::unordered_map<ReactionPtr, SCIP_VAR*> _reactions;
};

inline ScipModelPtr ReducedScipFluxModel::getScip() const {
	assert(!_origin.expired());
	return _origin.lock();
}

typedef boost::shared_ptr<ReducedScipFluxModel> ReducedScipFluxModelPtr;

} /* namespace metaopt */
#endif /* REDUCEDSCIPFLUXMODEL_H_ */
