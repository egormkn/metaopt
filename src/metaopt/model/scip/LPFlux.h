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
#include <boost/unordered_map.hpp>

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
	void setDirectionBounds(LPFluxPtr flux);
	/**
	 * Sets flux bounds according to the current solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to 1.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -1.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBounds(ScipModelPtr flux);

	/**
	 * Sets flux bounds according to the specified solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to 1.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -1.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBounds(SolutionPtr sol, ScipModelPtr flux);

	/**
	 * Sets flux bounds according to the current solution of specified LPFlux.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to +infinity.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -infinity.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBoundsInfty(LPFluxPtr flux);
	/**
	 * Sets flux bounds according to the current solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to +infinity.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -infinity.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBoundsInfty(ScipModelPtr flux);

	/**
	 * Sets flux bounds according to the specified solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to +infinity.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -infinity.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBoundsInfty(SolutionPtr sol, ScipModelPtr flux);

	/**
	 * sets the objective value such that reactions with positive flux are maximized,
	 * reactions with negative flux are minimized and reactions without flux have objective coef of 0.
	 */
	void setDirectionObj(LPFluxPtr flux);

	/**
	 * sets the objective value such that reactions with positive flux are maximized,
	 * reactions with negative flux are minimized and reactions without flux have objective coef of 0.
	 */
	void setDirectionObj(ScipModelPtr flux);

	/**
	 * sets the objective value such that reactions with positive flux (in the specified solution) are maximized,
	 * reactions with negative flux are minimized and reactions without flux have objective coef of 0.
	 */
	void setDirectionObj(SolutionPtr sol, ScipModelPtr flux);

	/**
	 * specifies, if we want to maximize (true) or minimize (false)
	 */
	void setObjSense(bool maximize);

	/**
	 * Computes optimal flux using primal simplex.
	 * Use this, if you only changed objective coefficients (or the objective sense).
	 */
	void solvePrimal();

	/**
	 * Computes optimal flux using dual simplex.
	 * Use this, if you only changed bounds.
	 */
	void solveDual();

	/**
	 * Computes optimal flux from scratch.
	 * Use this, if you modified many things.
	 */
	void solve();

	double getFlux(ReactionPtr rxn);

	double getDual(MetabolitePtr met);

	/**
	 * Is the current LP solution optimal ?
	 */
	bool isOptimal();

	/**
	 * Is the current LP solution feasible ?
	 */
	bool isFeasible();

	/**
	 * Returns current objective value.
	 * Attention: Does not start a solve, if the solution is not optimal.
	 */
	double getObjVal();

	/**
	 * print a human readable representation. For debugging.
	 */
	void print();

	/**
	 * copy the flux of the current solution of the ScipModel into this solution
	 */
	void set(ScipModelPtr);

	/**
	 * copy the flux of the specified solution of the ScipModel into this solution
	 */
	void set(SolutionPtr sol, ScipModelPtr);

	/**
	 * subtracts current solution of flux from this solution.
	 * Every Reaction in flux must also exist in this LPFlux.
	 * The converse needs not be true. If a reaction does not exist, a flux of 0 is assumed.
	 */
	void subtract(LPFluxPtr flux, double scale);

	/**
	 * compute maximal subtraction factor such that when source is subtracted from this flux, no flux direction changes.
	 * It is assumed that if a reaction carries nonzero flux in source, it also carries non-zero flux in this LPFlux in the same direction.
	 *
	 * Hence, the returned scale is always positive at successful termination.
	 * If no scale could be found -1 is returned.
	 */
	double getSubScale(LPFluxPtr source);

private:
	ModelPtr _model;
	boost::unordered_map<ReactionPtr, int> _reactions; // in the internal LP problem, columns are only identified by indices, so we have to map reactions to indices
	boost::unordered_map<MetabolitePtr, int> _metabolites; // in the internal LP problem, rows are only identified by indices, so we have to map metabolites to indices
	SCIP_LPI* _lpi; //< internal LP problem
	int _num_metabolites; //< number of metabolites in the LP model (depends if we include exchange fluxes or not)
	int _num_reactions; //< number of reactions in the LP model (depends if we include exchange fluxes or not)

	std::vector<double> _primsol; // use vector to store primal solution to circumvent deallocation hassle

	SCIP_RETCODE init_lp(bool exchange);
	SCIP_RETCODE free_lp();
};

typedef boost::shared_ptr<LPFlux> LPFluxPtr;

} /* namespace metaopt */
#endif /* LPFLUX_H_ */
