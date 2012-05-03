/*
 * DualPotentials.h
 *
 *  Created on: 02.05.2012
 *      Author: arnem
 */

#ifndef DUALPOTENTIALS_H_
#define DUALPOTENTIALS_H_

#include "ScipModel.h"
#include "LPFlux.h"

namespace metaopt {

/**
 * This class distinguishes from LPPotentials that it is designed to find small infeasible sets.
 *
 * Although, the LP formulation is derived from the dual of LPPotentials, it has some modifications.
 * These modifications would have created a mess, if implemented as an additional feature in LPPotentials.
 *
 * The problem that we will solve looks as follows:
 *
 * min \Eins \alpha_(N unfixed) - \Eins \alpha_(P unfixed)
 * 		u \beta^+ - \ell \beta^- + \gamma \leq 0
 * 		S \alpha + beta^+ - \beta^- = 0
 * 		\gamma + \Eins \alpha_N - \Eins \alpha_P = 1 \\
 * 		\alpha_N \geq 0
 * 		\alpha_P \leq 0
 * 		\alpha_(not N nor P) = 0
 * 		\beta^+, \beta^-, \gamma \geq 0
 *
 *
 */
class DualPotentials {
public:
	DualPotentials(ModelPtr model);
	virtual ~DualPotentials();

	/**
	 * set the directions of the reactions.
	 * The directions are set according to the specified flux.
	 *
	 * The problem will be initialized in such a way that a small infeasible set w.r.t. the specified configuration is found.
	 *
	 * Reactions that are already fixed (according to the parameter) do not count into the size of the infeasible set.
	 */
	void setDirections(LPFluxPtr flux, boost::shared_ptr<boost::unordered_set<ReactionPtr> > fixed_rxns);

	/**
	 * run the optuimization step
	 */
	void optimize();

	/**
	 * check, if a feasible solution (a violating set) has been found
	 */
	bool isFeasible();

	/**
	 * compute a small infeasible set
	 */
	boost::shared_ptr<boost::unordered_set<ReactionPtr> > getIS();

private:
	std::vector<double> _primsol; // use vector to store primal solution to circumvent deallocation hassle

	ModelPtr _model;
	boost::unordered_map<ReactionPtr, int> _reactions; // in the internal LP problem, columns are only identified by indices, so we have to map reactions to indices
	boost::unordered_map<MetabolitePtr, int> _metabolites; // in the internal LP problem, rows are only identified by indices, so we have to map metabolites to indices
	SCIP_LPI* _lpi; //< internal LP problem
	int _num_metabolites; //< number of metabolites in the LP model
	int _num_reactions; //< number of reactions in the LP model

	int _num_beta_vars;

	SCIP_RETCODE init_lp();
	SCIP_RETCODE free_lp();
};

typedef boost::shared_ptr<DualPotentials> DualPotentialsPtr;

} /* namespace metaopt */
#endif /* DUALPOTENTIALS_H_ */
