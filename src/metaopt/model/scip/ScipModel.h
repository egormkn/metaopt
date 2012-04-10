/*
 * ScipModel.h
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#ifndef SCIPMODEL_H_
#define SCIPMODEL_H_

#include "metaopt/model/Model.h"
#include "objscip/objscip.h"
#include <boost/unordered_map.hpp>

namespace metaopt {

struct PreconditionViolatedException : virtual boost::exception, virtual std::exception { };

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
	 */
	SCIP_VAR* getFlux(ReactionPtr rxn);
	/**
	 * fetches the scip variable representing the metabolite concentration.
	 * If the variable does not exist yet, it is created.
	 * The default bounds are the bounds stored by the metabolite.
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
	 * This solution is usually the value of the LP relaxation of the current node.
	 * Sometimes, the LP relaxation is not yet solved to feasibility and then there exists no current flux.
	 * It can also happen that no flux variables have been initialized at all. In this case, this method will always return false.
	 */
	bool hasCurrentFlux();

	/**
	 * Checks, if there exists a current solution for the potentials (of the current node in search tree).
	 * This solution is usually the value of the LP relaxation of the current node.
	 * Sometimes, the LP relaxation is not yet solved to feasibility and then there exists no current potentials.
	 * It can also happen that no potential variables have been initialized at all. In this case, this method will always return false.
	 */
	bool hasCurrentPotentials();

	/**
	 * Returns the flux value through the specified reaction.
	 *
	 * Precondition: hasCurrentFlux
	 *
	 * If the model has no current flux, a PreconditionViolatedException is thrown
	 */
	double getCurrentFlux(ReactionPtr rxn);

	/**
	 * Returns the potential of the specified metabolite.
	 *
	 * Precondition: hasCurrentPotential
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

private:
	ModelPtr _model;
	SCIP* _scip;

	SCIP_RETCODE init_scip();
	SCIP_RETCODE free_scip();

	boost::unordered_map<MetabolitePtr, SCIP_VAR*> _metabolites;
	boost::unordered_map<ReactionPtr, SCIP_VAR*> _reactions;
};

typedef boost::shared_ptr<ScipModel> ScipModelPtr;

/** If an error is thrown by a method of one of the getCurrentFlux or getCurrentPotential methods, this error_info states the reason */
typedef boost::error_info<struct tag_state,std::string> var_state;

} /* namespace metaopt */
#endif /* SCIPMODEL_H_ */
