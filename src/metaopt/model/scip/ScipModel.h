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
 * ScipModel.h
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#ifndef SCIPMODEL_H_
#define SCIPMODEL_H_

#include <boost/unordered_map.hpp>
#include <vector>
#include "objscip/objscip.h"

#include "metaopt/model/Model.h"
#include "metaopt/model/scip/AbstractScipFluxModel.h"
#include "Solution.h"
#include "metaopt/scip/ScipError.h"
#include "metaopt/Uncopyable.h"
#include "metaopt/Properties.h"

namespace metaopt {

struct PreconditionViolatedException : virtual boost::exception, virtual std::exception { };

class ModelAddOn;
typedef boost::shared_ptr<ModelAddOn> ModelAddOnPtr;

class ScipModel : public AbstractScipFluxModel, Uncopyable {
public:
	/** do not use, use factory method of Model instead! */
	ScipModel(ModelPtr model);
	virtual ~ScipModel();

	ModelPtr getModel() const;

	/** set if we want to maximize or minimize */
	inline void setObjectiveSense(bool maximize);

	/** returns true if objective sense is set to maximize.
	 * If objective sense is minimize, false is returned.
	 */
	inline bool isMaximize() const;

#if 0
	/**
	 * sets the solver to solve with the desired precision
	 */
	void setPrecision(double precision);

	/**
	 * returns the precision with which the problem will be solved.
	 */
	double getPrecision();
#endif

	/**
	 * fetches the scip variable representing the flux.
	 * If the variable does not exist yet, it is created.
	 * The default bounds are the bounds stored by the reaction.
	 *
	 * The pointer will be valid as long as this ScipModel lives.
	 */
	SCIP_VAR* getFlux(ReactionPtr rxn);
	/**
	 * fetches the scip variable representing the metabolite concentration.
	 * If the variable does not exist yet, it is created.
	 * The default bounds are the bounds stored by the metabolite.
	 *
	 * The pointer will be valid as long as this ScipModel lives.
	 */
	SCIP_VAR* getPotential(MetabolitePtr met);

	/**
	 * fetches the current upper bound of the flux.
	 * If the corresponding variable exists, its value of the current node is returned.
	 * Else, the default bound is returned.
	 */
	double getCurrentFluxUb(ReactionPtr rxn);
	/**
	 * fetches the current lower bound of the flux.
	 * If the corresponding variable exists, its value of the current node is returned.
	 * Else, the default bound is returned.
	 */
	double getCurrentFluxLb(ReactionPtr rxn);

	/**
	 * fetches the current upper bound of the potential.
	 * If the corresponding variable exists, its value of the current node is returned.
	 * Else, the default bound is returned.
	 */
	double getCurrentPotentialUb(MetabolitePtr met);
	/**
	 * fetches the current lower bound of the potential.
	 * If the corresponding variable exists, its value of the current node is returned.
	 * Else, the default bound is returned.
	 */
	double getCurrentPotentialLb(MetabolitePtr met);

	/**
	 * Checks, if there exists a current flux solution (of the current node in search tree).
	 * This is the case, if the problem has been solved, or if the LP relaxation of the current node has been solved to optimality.
	 * It can also happen that no flux variables have been initialized at all. In this case, this method will always return false.
	 */
	bool hasCurrentFlux();

	/**
	 * Checks, if there exists a current solution for the potentials (of the current node in search tree).
	 * This is the case, if the problem has been solved, or if the LP relaxation of the current node has been solved to optimality.
	 * It can also happen that no potential variables have been initialized at all. In this case, this method will always return false.
	 */
	bool hasCurrentPotentials();

	/**
	 * Returns the flux value through the specified reaction.
	 * If the problem has been solved, the value of an optimal solution is returned.
	 * Else, the value of the current LP-relaxation is returned. The LP-relaxation, however has to be solved to optimality or the primal ray proving unboundedness has been computed.
	 *
	 * Precondition: hasFluxVar(rxn),  hasCurrentFlux()
	 *
	 * If the model has no current flux, a PreconditionViolatedException is thrown
	 */
	double getCurrentFlux(ReactionPtr rxn);

	/**
	 * Returns the potential of the specified metabolite.
	 * If the problem has been solved, the value of an optimal solution is returned.
	 * Else, the value of the current LP-relaxation is returned. The LP-relaxation, however has to be solved to optimality.
	 *
	 * Precondition: hasPotentialVar(met), hasCurrentPotential()
	 *
	 * If the model has no current potentials, a PreconditionViolatedException is thrown
	 */
	double getCurrentPotential(MetabolitePtr met);


	/**
	 * Returns the flux value of specified solution through the specified reaction.
	 * If the solution is unbounded, the values of the ray are returned instead.
	 */
	double getFlux(SolutionPtr sol, ReactionPtr rxn);

	/**
	 * Returns the potential of the specified metabolite in the specified solution.
	 */
	double getPotential(SolutionPtr sol, MetabolitePtr met);

	/**
	 * checks if the variable representing the flux has already been created
	 */
	bool hasFluxVar(ReactionPtr rxn);
	/**
	 * checks if the variable representing the metabolite potential has already been created.
	 */
	bool hasPotentialVar(MetabolitePtr met);

	/**
	 * Tests if this model has any flux variables initialized
	 */
	bool hasFlux();

	/**
	 * Tests if this model has any potential variables initialized
	 */
	bool hasPotentials();

	/**
	 * Solves the problem that has been setup on this ScipModel.
	 *
	 * By default a scip model has no constraints etc.
	 * Constraints have to be added using separate constraint implementations.
	 */
	void solve();

	/**
	 * If the problem has been solved, the best objective value can be fetched using this method.
	 */
	double getObjectiveValue();

	/**
	 * If the problem has been solved, checks if the computed solution is optimal
	 */
	bool isOptimal();

	/**
	 * If the problem has been solved, checks if the computed solution is unbounded.
	 */
	bool isUnbounded();

	/**
	 * If the problem is unbounded, returns truw if also a primal ray is stored.
	 */
	bool hasPrimalRay();

	/**
	 * If the problem has been solved, returns true if the problem is infeasible
	 */
	bool isInfeasible();

	/**
	 * Returns the internal scip instance.
	 * The pointer is valid as long as this ScipModel lives.
	 */
	SCIP* getScip();

	/**
	 * register an addon.
	 *
	 * All the addons this addon depends on, must have been added before this one.
	 * It is recommended to add addons right after construction.
	 */
	void addAddOn(ModelAddOnPtr addon);

	/**
	 * try to compute values for variables introduced by the addons.
	 *
	 * @return true if successful, false otherwise
	 */
	bool computeAddOnValues(SolutionPtr sol);

	/**
	 * Set the direction of a reaction for the specified node of the search tree.
	 *
	 * This method is used for setting branching decisions.
	 * In particular, the specified direction does not only manipulate the direction of a reaction,
	 * but also the sign of the potential difference.
	 * This is not the same, because a reaction may be irreversible (allowing only forward flux),
	 * but can have negative potential difference (if it does not carry any flux).
	 *
	 * If fwd is set to true, the potential difference of the reaction is forced to be non-positive.
	 * This implies that backward flux is impossible.
	 * If fwd is set to false, the potential difference of the reaction is forced to be non-negative.
	 * This implies that forward flux is impossible.
	 *
	 * The pointer node must not be stored anywhere. It may only live as long as the function call goes.
	 *
	 * This method will call the setDirection method of the addons allowing addons to implement additional fixings.
	 */
	void setDirection(SCIP_NODE* node, ReactionPtr rxn, bool fwd);

	/**
	 * Block one flux in one direction of the specified reaction for the specified node of the search tree.
	 *
	 * This method is different to setDirection in the sense, that it does not say anything about the potential difference of the reaction.
	 * If you only know that this reaction cannot carry flux (e.g. by steady-state assumption etc.), use this method to propagate the information.
	 *
	 * If the information comes from a branching decision you may want to use setDirection.
	 *
	 * This method does not notify the addons!
	 *
	 * @param fwd set to true, to block forward directions
	 */
	void setBlockedFlux(SCIP_NODE* node, ReactionPtr rxn, bool fwd);

	/**
	 * returns all reactions that have fixed directions.
	 */
	boost::shared_ptr<boost::unordered_set<ReactionPtr> > getFixedDirections();

private:
	ModelPtr _model;
	SCIP* _scip;

	SCIP_RETCODE init_scip();
	SCIP_RETCODE free_scip();

	boost::unordered_map<MetabolitePtr, SCIP_VAR*> _metabolites; // map of initialized potential vars
	boost::unordered_map<ReactionPtr, SCIP_VAR*> _reactions; // map of initialized reaction vars

	/**
	 * A topologically sorted list of addons.
	 * For evaluation of the method computeSolutionVals, some of the addons require that
	 * this method was already executed for other addons.
	 * Cyclic dependencies are not allowed, hence we obtain a partial ordering.
	 * This vector contains the addons in a topologically sorted way w.r.t. this relation.
	 * Hence, all dependencies get satisfied, if the computeSolutionVals method is called from beginning (index 0) to end.
	 */
	std::vector<ModelAddOnPtr> _addons;

};

inline SCIP* ScipModel::getScip() {
	return _scip;
}

inline ModelPtr ScipModel::getModel() const {
	return _model;
}

inline void ScipModel::setObjectiveSense(bool maximize) {
	if(maximize) {
		BOOST_SCIP_CALL( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );
	}
	else {
		BOOST_SCIP_CALL( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MINIMIZE) );
	}
}

inline bool ScipModel::isMaximize() const {
	return SCIPgetObjsense(_scip) == SCIP_OBJSENSE_MAXIMIZE;
}

typedef boost::shared_ptr<ScipModel> ScipModelPtr;

/** If an error is thrown by a method of one of the getCurrentFlux or getCurrentPotential methods, this error_info states the reason */
typedef boost::error_info<struct tag_state,std::string> var_state;

} /* namespace metaopt */
#endif /* SCIPMODEL_H_ */
