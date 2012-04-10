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

} /* namespace metaopt */
#endif /* SCIPMODEL_H_ */
