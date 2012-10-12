/*
 * FVA.h
 *
 *  Created on: 14.07.2012
 *      Author: arne
 */

#ifndef FVA_H_
#define FVA_H_

#include <boost/unordered_map.hpp>
#include <utility>

#include "metaopt/model/Model.h"
#include "metaopt/model/scip/LPFlux.h"

#include "metaopt/algorithms/ModelFactory.h"
#include "metaopt/Properties.h"

namespace metaopt {

struct FVASettings {
	double timeout;
	boost::unordered_set<ReactionPtr> reactions;

	FVASettings() : timeout(-1), reactions() {};
};

typedef boost::shared_ptr<FVASettings> FVASettingsPtr;

/**
 * Runs ordinary flux variablity analysis on the given model. The objective coefficients are ignored, the whole flux space is analyzed.
 * Result is stored in the maps min and max. min contains the minimal possible flux, max contains the maximal possible flux.
 * If min,max are not empty, existing values may be overridden.
 */
void fva(ModelPtr model, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max );

/**
 * Runs ordinary flux variability analysis on the given LPFlux model.
 * Objective coefficients of model must be zero initially and will be zero afterwards.
 * The objective sense will be set to maximize.
 * Result is stored in the maps min and max. min contains the minimal possible flux, max contains the maximal possible flux.
 * If min,max are not empty, existing values may be overridden.
 */
void fva(LPFluxPtr model, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max );

/**
 * Runs flux variability on the given model.
 * Additional constraints, like thermodynamic constraints can be given via the model factory.
 * This solves a CIP for each reaction!
 */
void fva(ModelPtr model, ModelFactory& factory, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max );

/**
 * Runs thermodynamic FVA on the given model.
 * @param model the model to solve thermodynamic FVA on
 * @param reactions the reactions for which to do FVA
 * @param min minimal flux values
 * @param max maximal flux values
 */
void tfva(ModelPtr model, FVASettingsPtr settings, boost::unordered_map<ReactionPtr,double >& min , boost::unordered_map<ReactionPtr,double >& max );

/** Used if an reaction is not found */
struct TimeoutError : virtual boost::exception, virtual std::exception {};


} /* namespace metaopt */
#endif /* FVA_H_ */
