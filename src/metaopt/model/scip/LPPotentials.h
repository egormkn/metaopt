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

/**
 * LPPotentials models the space of feasible potentials.
 * I.e. we work on solutions that satisfy
 *
 * l <= mu <= u
 *   Dmu_i <  0 	i in fwd
 *   Dmu_i >  0		i in bwd
 *
 * Since this space is usually not closed, feasibility and so on are nontrivial matters.
 *
 * To test strict feasibility of Ax <= a && Bx < b, we test
 *
 * max z :
 *       Ax <= a
 * 	z1 + Bx <= b
 * 	z       <= 1
 *
 * In out application, A models the lower and upper bounds on mu.
 * It is assumed that l <= u (with no epsilon-tolerance!!); hence, Ax <= a is always feasible.
 * Thus, the strict-feasibility testing problem is always feasible.
 * Since z is bounded at 1, it always has a finite optimal solution.
 *
 * An usage cycle consists of the following steps:
 *
 * 0. Change bounds and objective
 *
 * 1. test for strict feasibility by solving the feasibility subproblem
 * 		a) load optimal basic solution of last feasibility test solve
 * 			We only changed bounds and one column
 * 		b) use dual simplex to compute the new solution
 * 		c) the computed solution is stored for the next iteration
 *
 * 2. compute the optimal solution
 * 		a) change to original objective function
 * 		b) use primal simplex to compute the solution
 *
 *
 * Possible problems:
 * 1. the column of the feastest variable is basic and is changed
 * 		can happen, if the directions are changed
 * 		basic solution may not be basic anymore
 *
 */
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
