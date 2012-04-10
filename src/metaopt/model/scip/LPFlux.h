/*
 * LPFlux.h
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#ifndef LPFLUX_H_
#define LPFLUX_H_

#include "scip/lpi.h"

#include "metaopt/model/Model.h"
#include "metaopt/model/scip/ScipModel.h"

namespace metaopt {

class LPFlux;
typedef boost::shared_ptr<LPFlux> LPFluxPtr;

class LPFlux {
public:
	/**
	 * creates a new LPFlux.
	 * If exchange is specified as true, the LP model also contains exchange reactions.
	 * Else, it only contains internal reactions.
	 */
	LPFlux(ModelPtr model, bool exchange);
	virtual ~LPFlux();

	/** set the lower bound of the specified variable to the specified value */
	void setLb(ReactionPtr r, double lb);
	/** set the upper bound of the specified variable to the specified value */
	void setUb(ReactionPtr r, double ub);
	/** set the objective coefficient of the specified variable to the specified value */
	void setObj(ReactionPtr r, double obj);
	/** set all objective coefficients to zero */
	void setZeroObj();

	/**
	 * Set the bounds to the same values as in the specified flux.
	 * Every reaction that appears in this flux, must also appear in the other flux.
	 */
	void setBounds(LPFluxPtr other);
	/**
	 * Set the bounds to the same values as in the specified flux.
	 * Every reaction that appears in this flux, must also appear in the other flux.
	 */
	void setBounds(ScipModelPtr other);

	/**
	 * Sets flux bounds according to the current solution of specified LPFlux.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to 1.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -1.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirections(LPFluxPtr flux);
	/**
	 * Sets flux bounds according to the current solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to 1.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -1.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirections(ScipModelPtr flux);

private:
	ModelPtr _model;
	boost::unordered_map<ReactionPtr, int> _reactions;
	SCIP_LPI* _lpi;

	SCIP_RETCODE init_lp(bool exchange);
	SCIP_RETCODE free_lp();
};


} /* namespace metaopt */
#endif /* LPFLUX_H_ */
