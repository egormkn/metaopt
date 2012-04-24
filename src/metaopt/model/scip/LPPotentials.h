/*
 * LPPotentials.h
 *
 *  Created on: 23.04.2012
 *      Author: arnem
 */

#ifndef LPPOTENTIALS_H_
#define LPPOTENTIALS_H_

#include "ScipModel.h"
#include "LPFlux.h"

namespace metaopt {

class LPPotentials {
public:
	LPPotentials(ModelPtr model);
	virtual ~LPPotentials();

	/**
	 * Sets reaction directions (bounds on potential differences) according to the current solution of specified LPFlux.
	 * If the reaction carries positive flux, the lower bound is set to 0 and the upper bound is set to Infinity.
	 * If the reaction carries negative flux, the upper bound is set to 0 and the lower bound is set to -Infinity.
	 * If the reaction carries no flux, the lower and upper bound are unconstrained.
	 */
	void setDirections(LPFluxPtr flux);

	void setDirection(ReactionPtr rxn, bool fwd);

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
	 * fetch computed potential
	 */
	double getPotential(MetabolitePtr met);

	/**
	 * checks, if the current solution is feasible
	 * does no computation.
	 * TODO: maybe change this!
	 */
	bool isFeasible();

private:
	ModelPtr _model;
	boost::unordered_map<ReactionPtr, int> _reactions; // in the internal LP problem, rows are only identified by indices, so we have to map reactions to indices
	boost::unordered_map<MetabolitePtr, int> _metabolites; // in the internal LP problem, columns are only identified by indices, so we have to map metabolites to indices
	SCIP_LPI* _lpi; //< internal LP problem
	int _num_metabolites; //< number of metabolites in the LP model
	int _num_reactions; //< number of reactions in the LP model

	std::vector<double> _primsol; // use vector to store primal solution to circumvent deallocation hassle

	SCIP_RETCODE init_lp();
	SCIP_RETCODE free_lp();

};

typedef boost::shared_ptr<LPPotentials> LPPotentialsPtr;

} /* namespace metaopt */
#endif /* LPPOTENTIALS_H_ */
