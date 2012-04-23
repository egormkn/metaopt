/*
 * CycleDeletionHeur.h
 *
 *  Created on: 23.04.2012
 *      Author: arnem
 */

#ifndef CYCLEDELETIONHEUR_H_
#define CYCLEDELETIONHEUR_H_

#include "objscip/objscip.h"
#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/model/scip/LPFlux.h"

namespace metaopt {

class CycleDeletionHeur : public scip::ObjHeur{
public:
	CycleDeletionHeur(ScipModelPtr _scip);
	virtual ~CycleDeletionHeur();

	/**
	 * interface method to scip
	 */
	virtual SCIP_RETCODE scip_exec(
		SCIP*              	scip,               /**< SCIP data structure */
		SCIP_HEUR*   	 	heur,
		SCIP_HEURTIMING    	heurtiming,
		SCIP_RESULT *    	result
		);

private:
	ScipModelPtr _scip;
	LPFluxPtr _difficultyTestFlux; // flux for checking if objective reactions are contained in internal circuits, contains only internal reactions
	LPFluxPtr _tflux; // flux for storing intermediate flux results, is not used for solving LPs, contains exchange reactions
	LPFluxPtr _cycle; // flux for finding cycles in _tflux, contains no exchange reactions

	/**
	 * checks if it will be hard to compute a non-trivial feasible solution
	 */
	bool isDifficult();

	/**
	 * computes a flux solution. Computed solution may not be feasible.
	 * This may be caused by flux forcing reactions.
	 * Additional effects produced by metabolite concentration information etc. are also not considered.
	 */
	bool computeFluxSolution(SolutionPtr sol);

	/**
	 * try to compute potentials that fit to the already computed flux
	 */
	bool computePotentials(SolutionPtr sol);
};

typedef boost::shared_ptr<CycleDeletionHeur> CycleDeletionHeurPtr;

} /* namespace metaopt */
#endif /* CYCLEDELETIONHEUR_H_ */
