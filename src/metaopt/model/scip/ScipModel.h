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
#include "Solution.h"

namespace metaopt {

struct PreconditionViolatedException : virtual boost::exception, virtual std::exception { };

class ModelAddOn;
typedef boost::shared_ptr<ModelAddOn> ModelAddOnPtr;

class ScipModel {
public:
	/** do not use, use factory method of Model instead! */
	ScipModel(ModelPtr model);
	virtual ~ScipModel();

	ModelPtr getModel() const;

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
	 * Else, the value of the current LP-relaxation is returned. The LP-relaxation, however has to be solved to optimality.
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
	 * checks if the variable representing the flux has already been created
	 */
	bool hasFluxVar(ReactionPtr rxn);
	/**
	 * checks if the variable representing the metabolite potential has already been created.
	 */
	bool hasPotentialVar(MetabolitePtr met);

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

typedef boost::shared_ptr<ScipModel> ScipModelPtr;

/** If an error is thrown by a method of one of the getCurrentFlux or getCurrentPotential methods, this error_info states the reason */
typedef boost::error_info<struct tag_state,std::string> var_state;

} /* namespace metaopt */
#endif /* SCIPMODEL_H_ */
