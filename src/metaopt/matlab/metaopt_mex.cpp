/*
 * metaopt_mex.cpp
 *
 *  Created on: 07.08.2012
 *      Author: arnem
 */

#include "metaopt_mex.h"

#include <boost/unordered_map.hpp>

#include "metaopt/Properties.h"
#include "metaopt/model/matlab/MatlabLoader.h"
#include "metaopt/algorithms/FVA.h"
#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/model/scip/LPFlux.h"
#include "metaopt/matlab/ExitEventHandler.h"
#include "metaopt/scip/constraints/SteadyStateConstraint.h"
#include "metaopt/scip/constraints/ThermoConstraintHandler.h"
#include "metaopt/scip/heur/CycleDeletionHeur.h"

namespace metaopt_mex {

using namespace boost;
using namespace std;
using namespace metaopt;

/**
 * Matlab wrapper to plain old FBA
 *
 * Second parameter must encode a metabolic model loadable by metaopt::Matlabloader
 */
void fba(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	assert(nrhs >= 1);

	if(nrhs == 1) {
		mexErrMsgTxt("fba: To run FBA you must specify a model as the second parameter!");
		return;
	}

	MatlabLoader loader;
	loader.load(prhs[1]);
	ModelPtr model = loader.getModel();
	LPFluxPtr flux(new LPFlux(model, true));
	flux->solve();

	if(nlhs >= 1) {
		plhs[0] = mxCreateDoubleScalar(flux->getObjVal());
	}
	if(nlhs >= 2) {
		plhs[1] = mxCreateDoubleMatrix(model->getReactions().size(), 1, mxREAL);
		double* fluxSolData = mxGetPr(plhs[1]);
		int numRxns = model->getReactions().size();
		for(int i = 0; i < numRxns; i++) {
			ReactionPtr rxn = loader.getReaction(i);
			fluxSolData[i] = flux->getFlux(rxn);
		}
	}
}

/**
 * Matlab wrapper to thermodynamically constrained FBA
 */
void tfba(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	assert(nrhs >= 1);

	if(nrhs == 1) {
		mexErrMsgTxt("tfba: To run FBA you must specify a model as the second parameter!");
		return;
	}

	MatlabLoader loader;
	loader.load(prhs[1]);
	ModelPtr model = loader.getModel();
	ScipModelPtr scip(new ScipModel(model));
	createSteadyStateConstraint(scip);
	createThermoConstraint(scip);
	createCycleDeletionHeur(scip);
	registerExitEventHandler(scip);

	scip->solve();

	if(nlhs >= 1) {
		plhs[0] = mxCreateDoubleScalar(scip->getObjectiveValue());
	}
	if(nlhs >= 2) {
		plhs[1] = mxCreateDoubleMatrix(model->getReactions().size(), 1, mxREAL);
		double* fluxSolData = mxGetPr(plhs[1]);
		int numRxns = model->getReactions().size();
		for(int i = 0; i < numRxns; i++) {
			ReactionPtr rxn = loader.getReaction(i);
			fluxSolData[i] = scip->getCurrentFlux(rxn);
		}
	}
}

/**
 * Stores fva results in the matlab structure (as a matrix with two columns)
 */
void convert_fva_result(mxArray** m, metaopt::ModelPtr& model, metaopt::MatlabLoader& loader, boost::unordered_map<metaopt::ReactionPtr, double>& min, boost::unordered_map<metaopt::ReactionPtr, double>& max) {
	*m = mxCreateDoubleMatrix(model->getReactions().size(), 2, mxREAL);
	double* outData = mxGetPr(*m);
	mwIndex subs[2];
	bool warned = false;
	int numRxns = model->getReactions().size();
	for(int i = 0; i < numRxns; i++) {
		ReactionPtr rxn = loader.getReaction(i);
		unordered_map<ReactionPtr, double >::iterator iter_min = min.find(rxn);
		unordered_map<ReactionPtr, double >::iterator iter_max = max.find(rxn);
		subs[0] = i;
		subs[1] = 0;
		if(iter_min == min.end()) {
			if(!warned) {
				mexWarnMsgTxt("Error: Did not run FVA for a reaction - computation aborted?\n");
				//mexWarnMsgTxt("setting lower and upper bounds to +-Infinity\n");
				warned = true;
			}
			outData[mxCalcSingleSubscript(*m, 2, subs)] = -INFINITY;
		}
		else {
			outData[mxCalcSingleSubscript(*m, 2, subs)] = iter_min->second;
		}
		subs[1] = 1;
		if(iter_min == min.end()) {
			if(!warned) {
				mexWarnMsgTxt("Error: Did not run FVA for a reaction - computation aborted?\n");
				//mexWarnMsgTxt("setting lower and upper bounds to +-Infinity\n");
				warned = true;
			}
			outData[mxCalcSingleSubscript(*m, 2, subs)] = INFINITY;
		}
		else {
			outData[mxCalcSingleSubscript(*m, 2, subs)] = iter_max->second;
		}
	}
}

/**
 * Matlab wrapper to plain old FVA
 */
void fva(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	assert(nrhs >= 1);

	if(nrhs == 1) {
		mexErrMsgTxt("fva: To run FVA you must specify a model as the second parameter!");
		return;
	}

	if(nlhs == 0) {
		mexErrMsgTxt("fva: For FVA you must at least give one output parameter");
		return;
	}

	MatlabLoader loader;
	loader.load(prhs[1]);
	ModelPtr model = loader.getModel();

	unordered_map<ReactionPtr, double> min, max;

	fva(model, min, max);

	convert_fva_result(&(plhs[0]), model, loader, min, max);
}

class ThermoModelFactory : public ModelFactory {
public:
	ScipModelPtr build(ModelPtr m) {
		ScipModelPtr scip(new ScipModel(m));
		createSteadyStateConstraint(scip);
		createThermoConstraint(scip);
		createCycleDeletionHeur(scip);
		registerExitEventHandler(scip);

		return scip;
	}
};

/**
 * Matlab wrapper to thermodynamically constrained FVA
 */
void tfva(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	assert(nrhs >= 1);

	if(nrhs == 1) {
		mexErrMsgTxt("fva: To run FVA you must specify a model as the second parameter!");
		return;
	}

	if(nlhs == 0) {
		mexErrMsgTxt("fva: For FVA you must at least give one output parameter");
		return;
	}

	MatlabLoader loader;
	loader.load(prhs[1]);
	ModelPtr model = loader.getModel();

	unordered_map<metaopt::ReactionPtr, double> min, max;

	ThermoModelFactory factory;

	fva(model, factory, min, max);

	convert_fva_result(&(plhs[0]), model, loader, min, max);
}

};

using namespace metaopt_mex;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if(nrhs == 0) {
		mexErrMsgTxt("You must give at least the name of the function you want to call!\n");
		return;
	}

	if(!mxIsChar(prhs[0])) {
		mexErrMsgTxt("The first parameter must be the name of the function you want to call!\n");
		return;
	}

	char* fname = mxArrayToString(prhs[0]);

	if(strcmp(fname, "fba") == 0) {
		fba(nlhs, plhs, nrhs, prhs);
	}
	else if(strcmp(fname, "tfba") == 0) {
		tfba(nlhs, plhs, nrhs, prhs);
	}
	else if(strcmp(fname, "fva") == 0) {
		fva(nlhs, plhs, nrhs, prhs);
	}
	else if(strcmp(fname, "tfva") == 0) {
		tfva(nlhs, plhs, nrhs, prhs);
	}
}

