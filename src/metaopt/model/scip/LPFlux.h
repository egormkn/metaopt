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
#include "metaopt/Uncopyable.h"
#include "metaopt/model/scip/ISSupply.h"
#include "metaopt/model/scip/PotSpaceConstraint.h"
#include "metaopt/Properties.h"

namespace metaopt {

class LPFlux;
typedef boost::shared_ptr<LPFlux> LPFluxPtr;

class LPFlux : Uncopyable, public ISSupply {
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

	/** gets the lower bound on the specified reaction */
	double getLb(ReactionPtr r);

	/** gets the upper bound on the specified reaction */
	double getUb(ReactionPtr r);

	/** gets the objective on the specified reaction */
	double getObj(ReactionPtr r);

	/**
	 * Set the bounds to the same values as in the specified flux.
	 * Every reaction that appears in this flux, must also appear in the other flux.
	 */
	void setBounds(LPFluxPtr other);
	/**
	 * Set the bounds to the same values as in the specified flux.
	 * Every reaction that appears in this flux, must also appear in the other flux.
	 */
	void setBounds(AbstractScipFluxModelPtr other);

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
	void setDirectionBounds(AbstractScipFluxModelPtr flux);

	/**
	 * Sets flux bounds according to the specified solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to 1.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -1.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBounds(SolutionPtr sol, AbstractScipFluxModelPtr flux);

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
	void setDirectionBoundsInfty(AbstractScipFluxModelPtr flux);

	/**
	 * Sets flux bounds according to the specified solution of specified ScipModel.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to +infinity.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -infinity.
	 * If the reaction carries no flux, the lower and upper bound are set to 0.
	 */
	void setDirectionBoundsInfty(SolutionPtr sol, AbstractScipFluxModelPtr flux);

	/**
	 * sets the objective value such that reactions with positive flux are maximized,
	 * reactions with negative flux are minimized and reactions without flux have objective coef of 0.
	 */
	void setDirectionObj(LPFluxPtr flux);

	/**
	 * sets the objective value such that reactions with positive flux are maximized,
	 * reactions with negative flux are minimized and reactions without flux have objective coef of 0.
	 */
	void setDirectionObj(AbstractScipFluxModelPtr flux);

	/**
	 * sets the objective value such that reactions with positive flux (in the specified solution) are maximized,
	 * reactions with negative flux are minimized and reactions without flux have objective coef of 0.
	 */
	void setDirectionObj(SolutionPtr sol, AbstractScipFluxModelPtr flux);

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

	/**
	 * resets this LPFlux. The next solve will be from scratch.
	 */
	void resetState();

	double getFlux(ReactionPtr rxn);

	/**
	 * If LPFlux is used to compute infeasible sets, the flux on the reactions gives the corresponding coefficients.
	 * This means, this method gives the same results as getFlux.
	 */
	double getAlpha(ReactionPtr rxn);

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
	 * writes the current LP into a file
	 */
	void write(const char* filename);

	/**
	 * writes the current state of the LP into a file (i.e. basis information)
	 */
	void writeState(const char* filename);

	/**
	 * copy the flux of the current solution of the ScipModel into this solution
	 */
	void set(AbstractScipFluxModelPtr);

	/**
	 * copy the flux of the specified solution of the ScipModel into this solution
	 */
	void set(SolutionPtr sol, AbstractScipFluxModelPtr);

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

	/**
	 * fetches the model that generated this LPFlux
	 */
	inline ModelPtr getModel();

	/**
	 * updates list of extra constraints
	 */
	void setExtraPotConstraints(boost::unordered_set<PotSpaceConstraintPtr>& psc);

#ifndef NDEBUG
	// these methods are for debugging only! A state of the LP can be stored and fetched later on
	std::vector<int> cstat;
	std::vector<int> rstat;

	void loadState();

	void setOldState();

	void initStateInfo();

	SCIP_LPI* getLPI();
	int getIndex(ReactionPtr rxn);
#endif

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

	boost::unordered_map<PotSpaceConstraintPtr, int> _extraConstraints;

};

inline ModelPtr LPFlux::getModel() {
	return _model;
}

typedef boost::shared_ptr<LPFlux> LPFluxPtr;

} /* namespace metaopt */
#endif /* LPFLUX_H_ */
