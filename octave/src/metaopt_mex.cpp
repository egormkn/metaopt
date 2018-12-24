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
 * metaopt_mex.cpp
 *
 *  Created on: 07.08.2012
 *      Author: arnem
 */

#include "metaopt_mex.h"

#include <iostream>
#include <inttypes.h>
#include <boost/unordered_map.hpp>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "Properties.h"
#include "model/octave/OctaveLoader.h"
#include "algorithms/FVA.h"
#include "algorithms/BlockingSet.h"
#include "model/scip/ScipModel.h"
#include "model/scip/LPFlux.h"
#include "scip/constraints/SteadyStateConstraint.h"
#include "scip/constraints/ThermoConstraintHandler.h"
#include "scip/heur/CycleDeletionHeur.h"
#include "model/Precision.h"

//#include "ExitEventHandler.h"

namespace metaopt_mex {

using namespace boost;
using namespace std;
using namespace metaopt;

/**
 * Creates the extra linear constraints A v <= a
 */
void createExtraConstraint(ScipModelPtr scip, OctaveLoader& loader, mxArray* A, mxArray* a) {
	assert(A != NULL); // A does not exist
	assert(a != NULL); // a does not exist
	assert(!mxIsChar(A) && !mxIsComplex(A) && mxIsNumeric(A) && mxIsSparse(A)); // A must be sparse and contain reals
	assert(!mxIsChar(a) && !mxIsComplex(a) && mxIsNumeric(a) && !mxIsSparse(a)); // a must be full and contain reals

	double* data = mxGetPr(A);
	mwIndex* rows = mxGetIr(A);
	mwIndex* pos = mxGetJc(A);

	int numRows = mxGetM(A);
	int numCols = mxGetN(A);

	assert((int) mxGetNumberOfElements(a) == numRows); // a does not have as many elements as A has rows

	vector<vector<SCIP_VAR*> > vars(numRows);
	vector<vector<double> > coef(numRows);

	ModelPtr model = scip->getModel();
	PrecisionPtr prec = model->getCoefPrecision();

	for(int i = 0; i < numCols; i++) {
		ReactionPtr rxn = loader.getReaction(i);
		SCIP_VAR* var = scip->getFlux(rxn);
		for(mwIndex j = pos[i]; j < pos[i+1]; j++) {
			if(data[j] < -prec->getCheckTol() || data[j] > prec->getCheckTol()) {
				vars.at(rows[j]).push_back(var);
				coef.at(rows[j]).push_back(data[j]);
				rxn->setProblematic(true);
			}
		}
	}

	double* aData = mxGetPr(a);
	double lhs = -INFINITY;

	for(int i = 0; i < numRows; i++) {
		vector<SCIP_VAR*>& cvars = vars[i];
		vector<double>& ccoef = coef[i];
		SCIP_CONS* cons = NULL;
		BOOST_SCIP_CALL( SCIPcreateConsLinear(scip->getScip(), &cons, "additional cons", cvars.size(), cvars.data(), ccoef.data(), lhs, aData[i], true, true, true,true, true, false, false, false, false, false) );
		BOOST_SCIP_CALL( SCIPaddCons(scip->getScip(), cons) );
		BOOST_SCIP_CALL( SCIPreleaseCons(scip->getScip(), &cons) );
	}
}

/**
 * Octave wrapper to plain old FBA
 *
 * Second parameter must encode a metabolic model loadable by metaopt::Octaveloader
 */
int fba(mxArray *argout, mxArray *argin) {
	int nargin = mxGetNumberOfElements(argin);

	if(nargin == 0) {
		cout << "fba: To run FBA you must specify a model as the second parameter!" << endl;
		return 16;
	}

	OctaveLoader loader;
	loader.load(mxGetCell(argin, 0));
	ModelPtr model = loader.getModel();

	LPFluxPtr flux(new LPFlux(model, true));
	flux->solve();

	int nargout = mxGetNumberOfElements(argout);

	if(nargout >= 1) {
		mxArray* val = mxCreateDoubleScalar(flux->getObjVal());
		mxSetCell(argout, 0, val);
		//mxDestroyArray(val);
	}
	if(nargout >= 2) {
		mxArray* fluxA = mxCreateDoubleMatrix(model->getReactions().size(), 1, mxREAL);
		double* fluxSolData = mxGetPr(fluxA);
		int numRxns = model->getReactions().size();
		for(int i = 0; i < numRxns; i++) {
			ReactionPtr rxn = loader.getReaction(i);
			fluxSolData[i] = flux->getFlux(rxn);
		}
		mxSetCell(argout, 1, fluxA);
		//mxDestroyArray(fluxA);
	}

	return 0;
}

/**
 * Octave wrapper to thermodynamically constrained FBA
 */
int tfba(mxArray *argout, mxArray *argin) {
	int nargin = mxGetNumberOfElements(argin);

	if(nargin == 0) {
		cout << "fba: To run FBA you must specify a model as the second parameter!" << endl;
		return 17;
	}

	OctaveLoader loader;
	mxArray* sys = mxGetCell(argin, 0);
	loader.load(sys);
	ModelPtr model = loader.getModel();
	ScipModelPtr scip(new ScipModel(model));
	createSteadyStateConstraint(scip);
	mxArray* extra_A = mxGetField(sys, 0, "A");
	mxArray* extra_a = mxGetField(sys, 0, "a");
	if(extra_A != NULL) {
		assert(extra_a != NULL);
		createExtraConstraint(scip, loader, extra_A, extra_a);
	}
	createThermoConstraint(scip);
	createCycleDeletionHeur(scip);
	//registerExitEventHandler(scip);

	scip->solve();

	bool solFound = scip->isOptimal();

	int nargout = mxGetNumberOfElements(argout);

	if(nargout >= 1) {
		if(solFound) {
			mxSetCell(argout, 0, mxCreateDoubleScalar(scip->getObjectiveValue()));
		}
		else if(scip->isUnbounded()) {
			if(scip->isMaximize()) {
				mxSetCell(argout, 0, mxCreateDoubleScalar(INFINITY));
			}
			else {
				mxSetCell(argout, 0, mxCreateDoubleScalar(-INFINITY));
			}
		}
		else {
			if(!scip->isInfeasible()) {
				cout << "Warning: Did not compute a solution, but not infeasible nor unbounded" << endl;
			}
			mxSetCell(argout, 0, mxCreateDoubleScalar(NAN));
		}
	}
	if(nargout >= 2) {
		if(solFound) {
			mxArray* fluxA = mxCreateDoubleMatrix(model->getReactions().size(), 1, mxREAL);
			double* fluxSolData = mxGetPr(fluxA);
			int numRxns = model->getReactions().size();
			for(int i = 0; i < numRxns; i++) {
				ReactionPtr rxn = loader.getReaction(i);
				fluxSolData[i] = scip->getCurrentFlux(rxn);
			}
			mxSetCell(argout, 1, fluxA);
		}
		else {
			mxArray* fluxA = mxCreateDoubleMatrix(0, 0, mxREAL);
			mxSetCell(argout, 1, fluxA);
		}
	}
	if(nargout >= 3 && (solFound || scip->hasPrimalRay())) {
		LPFluxPtr flux(new LPFlux(model, true));
		//LPFluxPtr m(new LPFlux(model, true));
		flux->set(scip); // copy optimal solution
		boost::shared_ptr<std::vector<DirectedReaction> > block = findBlockingSet(flux);
		int num_fwd = 0;
		int num_bwd = 0;
		foreach(DirectedReaction& d, *block) {
			if(d._fwd) num_fwd++;
			else num_bwd++;
		}

		cout << "num fwd reactions in blocking set: " << num_fwd << endl;
		cout << "num bwd reactions in blocking set: " << num_bwd << endl;

		mxArray* blockFwd = mxCreateNumericMatrix(num_fwd, 1, mxUINT32_CLASS, mxREAL);
		mxArray* blockBwd = mxCreateNumericMatrix(num_bwd, 1, mxUINT32_CLASS, mxREAL);

		// it is important to get the same type, else the indices of the array don't match
		uint32_t* fwd = (uint32_t*) mxGetData(blockFwd);
		uint32_t* bwd = (uint32_t*) mxGetData(blockBwd);

		unordered_map<ReactionPtr, int> rxnIndex;
		int numRxns = model->getReactions().size();
		for(int i = 0; i < numRxns; i++) {
			rxnIndex[loader.getReaction(i)] = i;
		}

		int i_fwd = 0;
		int i_bwd = 0;
		foreach(DirectedReaction& d, *block) {
			if(d._fwd) {
				assert(rxnIndex.find(d._rxn) != rxnIndex.end());
				fwd[i_fwd] = rxnIndex[d._rxn] + 1; // octave indexing starts at 1
				cout << "blocked fwd reaction: " << d._rxn->getName() << " (" << rxnIndex[d._rxn] << ")" << endl;
				i_fwd++;
			}
			else {
				assert(rxnIndex.find(d._rxn) != rxnIndex.end());
				bwd[i_bwd] = rxnIndex[d._rxn] + 1; // octave indexing starts at 1
				cout << "blocked bwd reaction: " << d._rxn->getName() << " (" << rxnIndex[d._rxn] << ")" << endl;
				i_bwd++;
			}
		}
		assert(i_fwd == num_fwd);
		assert(i_bwd == num_bwd);

		const char* fwd_name = "fwd";
		const char* bwd_name = "bwd";
		const char* fieldnames[2] = {fwd_name, bwd_name};

		mxArray* block_m = mxCreateStructMatrix(1,1,2,fieldnames);
		mxSetField(block_m, 0, fwd_name, blockFwd);
		mxSetField(block_m, 0, bwd_name, blockBwd);

		mxSetCell(argout, 2, block_m);
	}

	return 0;
}

/**
 * Stores fva results in the octave structure (as a matrix with two columns)
 * The order of the elements is the same as they were loaded by the octave loader
 *
 * @param m the result matrix that is created
 * @param model the model on which FVA had been run
 * @param loader the loader with which the model had been loaded
 * @param min the minimal flux values
 * @param max the maximal flux values
 */
void convert_fva_result(mxArray** m, metaopt::ModelPtr& model, metaopt::OctaveLoader& loader, boost::unordered_map<metaopt::ReactionPtr, double>& min, boost::unordered_map<metaopt::ReactionPtr, double>& max) {
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
				cout << "Error: Did not run FVA for a reaction - computation aborted?" << endl;
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
				cout << "Error: Did not run FVA for a reaction - computation aborted?" << endl;
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
 * Stores fva results in the octave structure (as a matrix with two columns)
 * The order of the elements is the same as given in the reaction list
 *
 * @param m the result matrix that is created
 * @param model the model on which FVA had been run
 * @param loader the loader with which the model had been loaded
 * @param rxnList a octave array containing the octave-indices (starting from 1!) of the reactions for which FVA was performed
 * @param min the minimal flux values
 * @param max the maximal flux values
 */
void convert_fva_result(mxArray** m, metaopt::ModelPtr& model, metaopt::OctaveLoader& loader, mxArray* rxnList, boost::unordered_map<metaopt::ReactionPtr, double>& min, boost::unordered_map<metaopt::ReactionPtr, double>& max) {
	int num = mxGetNumberOfElements(rxnList);
	assert( mxIsDouble(rxnList) );
	double* rxnListData = mxGetPr(rxnList);

	*m = mxCreateDoubleMatrix(num, 2, mxREAL);
	double* outData = mxGetPr(*m);
	mwIndex subs[2];
	bool warned = false;
	for(int i = 0; i < num; i++) {
		int rxnId = (int) (rxnListData[i]+0.5);
		assert(rxnId >= 1);

		ReactionPtr rxn = loader.getReaction(rxnId-1);
		unordered_map<ReactionPtr, double >::iterator iter_min = min.find(rxn);
		unordered_map<ReactionPtr, double >::iterator iter_max = max.find(rxn);
		subs[0] = i;
		subs[1] = 0;
		if(iter_min == min.end()) {
			if(!warned) {
				cout << "Error: Did not run FVA for a reaction - computation aborted?" << endl;
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
				cout << "Error: Did not run FVA for a reaction - computation aborted?" << endl;
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
 * Octave wrapper to plain old FVA
 */
int fva(mxArray *argout, mxArray *argin) {
	int nargin = mxGetNumberOfElements(argin);

	if(nargin == 0) {
		cout << "fva: To run FVA you must specify a model as the second parameter!" << endl;
		return 18;
	}

	OctaveLoader loader;
	loader.load(mxGetCell(argin, 0));
	ModelPtr model = loader.getModel();

	unordered_map<ReactionPtr, double> min, max;

	fva(model, min, max);

	mxArray* res;

	convert_fva_result(&res, model, loader, min, max);
	mxSetCell(argout, 0, res);

	return 0;
}

/*class ThermoModelFactory : public ModelFactory {
public:
	ScipModelPtr build(ModelPtr m) {
		ScipModelPtr scip(new ScipModel(m));
		createSteadyStateConstraint(scip);
		createThermoConstraint(scip);
		createCycleDeletionHeur(scip);
		//registerExitEventHandler(scip);

		return scip;
	}
};*/

/**
 * Octave wrapper to thermodynamically constrained FVA
 */
int tfva(mxArray *argout, mxArray *argin) {
	int nargin = mxGetNumberOfElements(argin);

	if(nargin == 0) {
		cout << "tfva: To run FVA you must specify a model as the second parameter!" << endl;
		return 19;
	}
	OctaveLoader loader;
	loader.load(mxGetCell(argin, 0));
	ModelPtr model = loader.getModel();

	// FVA settings
	FVASettingsPtr settings(new FVASettings());

	mxArray* rxnList = NULL;

	if(nargin == 1) {
		// default settings, if no setting parameter is given
		settings->reactions = model->getReactions();
	}
	else {
		// load settings from second parameter
		mxArray* set = mxGetCell(argin, 1);

		rxnList = mxGetField(set, 0, "reactions");
		mxArray* timeout = mxGetField(set, 0, "timeout");

		if(rxnList != NULL) {
			// load reaction list from second parameter
			if( !mxIsNumeric(rxnList) ) {
				cout << "reactions must contain indices of reactions" << endl;
				return 19;
			}
			int num = mxGetNumberOfElements(rxnList);
			double* rxnListData = mxGetPr(rxnList);
			for(int i = 0; i < num; i++) {
				int rxnId = (int) (rxnListData[i]+0.5);
				if(rxnId <= 0) {
					cout << "reaction indices must be greater then 0!" << endl;
					return 19;
				}
				settings->reactions.insert( loader.getReaction(rxnId-1) ); // subtract 1 to get C++ index
			}
		}
		else {
			settings->reactions = model->getReactions();
		}
		if(timeout != NULL) {
			if( !mxIsNumeric(timeout) ) {
				cout << "timeout parameter must be a number" << endl;
				return 19;
			}
			settings->timeout = mxGetScalar(timeout);
		}
	}

	unordered_map<metaopt::ReactionPtr, double> min, max;

	//ThermoModelFactory factory;

	try {
		metaopt::tfva(model, settings, min, max);
	} catch(std::exception& ex) {
		std::cout << diagnostic_information(ex) << std::endl;
		return 19;
	}

	mxArray* res;

	if(rxnList == NULL) {
		// use order of the loader
		convert_fva_result(&res, model, loader, min, max);
	}
	else {
		// use given order
		convert_fva_result(&res, model, loader, rxnList, min, max);
	}

	mxSetCell(argout, 0, res);

	return 0;
}

int help(mxArray *argout, mxArray *argin) {
	cout << "Metaopt Version " << VERSION << endl;
	cout << "Copyright (C) 2012 Arne Müller <arne.mueller@fu-berlin.de> " << endl;
	cout << endl;
	cout << "This program is free software: you can redistribute it and/or modify" << endl;
	cout << "it under the terms of the GNU General Public License as published by" << endl;
	cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
	cout << "(at your option) any later version." << endl;
	cout << endl;
	cout << "This program is distributed in the hope that it will be useful," << endl;
	cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
	cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
	cout << "GNU General Public License for more details." << endl;
	cout << endl;
	cout << "You should have received a copy of the GNU General Public License" << endl;
	cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
	cout << endl;
	cout << endl;
	cout << "Usage:" << endl;
	cout << "  <output parameters> = metaopt(<function name>, <input parameters>)" << endl;
	cout << "  <output parameters> = metaopt_time(<function name>, <timeout>, <input parameters>)" << endl;
	cout << "- function name has to be one of the functions listed below." << endl;
	cout << "- timout has has to be a number followed by either s, m, h, or d for seconds, minutes, hours, or days respectively." << endl;
	cout << "- input parameters depend on the function called." << endl;
	cout << "- output parameters depend on the function called." << endl;
	cout << endl;
	cout << "Available functions for metaopt:" << endl;
	cout << endl;
	cout << "fba: for ordinary fba (without thermodynamic constraints)." << endl;
	cout << "  Input Parameters:" << endl;
	cout << "   - a metabolic network model. See below for a precise specification." << endl;
	cout << "  Output Parameters:" << endl;
	cout << "   - optimal objective value (scalar)" << endl;
	cout << "   - an optimal flux distribution (column-vector)" << endl;
	cout << endl;
	cout << "tfba: for tfba (with thermodynamic constraints)." << endl;
	cout << "  Input Parameters:" << endl;
	cout << "   - a metabolic network model. See below for a precise specification." << endl;
	cout << "  Output Parameters:" << endl;
	cout << "   - optimal objective value (scalar)" << endl;
	cout << "   - an optimal flux distribution (column-vector)" << endl;
	cout << "   - blocking set. This is a struct with two fields 'fwd', 'bwd'." << endl;
	cout << "      - fwd is a list of reactions to be blocked in forward direction" << endl;
	cout << "      - bwd is a list of reactions to be blocked in backward direction" << endl;
	cout << "     If this blocking set is enforced, the tfba objective value will equal the fba objective value (only works without potential bounds)." << endl;
	cout << endl;
	cout << "fva: for fva (without thermodynamic constraints). This ignores any specified objective function." << endl;
	cout << "     If you want to run fva only on the optimal flux space, you have to make a separate call to fba and constrain the flux space accordingly." << endl;
	cout << "  Input Parameters:" << endl;
	cout << "   - a metabolic network model. See below for a precise specification." << endl;
	cout << "  Output Parameters:" << endl;
	cout << "   - variabilities of fluxes (matrix). " << endl;
	cout << "     The first column gives the minimal flux for each reaction." << endl;
	cout << "     The second column gives maximal flux for each reaction." << endl;
	cout << endl;
	cout << "tfva: for tfva (with thermodynamic constraints). This ignores any specified objective function." << endl;
	cout << "      If you want to run fva only on the optimal flux space, you have to make a separate call to tfba and constrain the flux space accordingly." << endl;
	cout << "  Input Parameters:" << endl;
	cout << "   - a metabolic network model. See below for a precise specification." << endl;
	cout << "   - a struct settings wich can have the following parameters:" << endl;
	cout << "      - reactions: Specifies a list (indices) of reactions for which fva should be performed." << endl;
	cout << "      - timeout: This feature is deprecated and does not work reliably. Use metaopt_timeout instead." << endl;
	cout << "  Output Parameters:" << endl;
	cout << "   - variabilities of fluxes (matrix). " << endl;
	cout << "     The first column gives the minimal flux for each reaction." << endl;
	cout << "     The second column gives maximal flux for each reaction." << endl;
	cout << endl;
	cout << "help: Prints this message." << endl;
	cout << endl;
	cout << "Specification of metabolic network model:" << endl;
	cout << "  The metabolic network model is struct which is basically a COBRA model with additional fields:" << endl;
	cout << "   - S         stoichiometric matrix (must be a sparse matrix)" << endl;
	cout << "   - ub        flux upper bounds" << endl;
	cout << "   - lb        flux lower bounds" << endl;
	cout << "   - c         objective function on fluxes" << endl;
	cout << "   - rxns      cell array of reaction names" << endl;
	cout << "   - mets      cell array of metabolite names" << endl;
	cout << "   - boundary  list of indices specifying metabolites with boundary condition, i.e. steady-state assumption does not need to hold for this metabolite (experimental)." << endl;
	cout << "   - p_ub      potential upper bounds (experimental)" << endl;
	cout << "   - p_lb      potential lower bounds (experimental)" << endl;
	cout << "   - p_c       potential objective function (not yet used)" << endl;
	cout << "   - int_rxns  cell array of internal reaction names" << endl;
	cout << "   - precision precision with which the model should be solved" << endl;
	cout << "   - A         coefficient matrix (sparse) for extra constraints A v <= a (only for tfba)" << endl;
	cout << "   - a         left hand sides of extra constraints A v <= a (only for tfba)" << endl;
	cout << "  Remark: The right-hand side vector b (of S v = b) is assumed to be 0." << endl;
	cout << endl;
	cout << "Please report bugs, feature requests etc. to me: arne.mueller@fu-berlin.de" << endl;

	if(argout != NULL) {
		mxArray* ans = mxCreateString("If you did not see a help message, output to stdout of metaopt is not shown to you. To see the help message, simply call metaopt from the terminal without parameters.");
		mxSetCell(argout, 0, ans);
	}

	return 0;
}

};

using namespace metaopt_mex;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if(nrhs < 2) {
		mexPrintf ("Usage: \n");
		mexPrintf ("metaopt(funcname, inarg)\n");
		mexPrintf ("infilename: name of the file containing input arguments\n");
		mexPrintf ("outfilename: name of the file containing where the output should be written to\n");
		mexPrintf ("funcname: name of the function to be called\n");
		mexPrintf ("nargout: number of output arguments to be produced\n");
		mexPrintf ("\n");
		// help(NULL, NULL);
		return;
	}

	char* fname = mxArrayToString(prhs[0]);
	mxArray* argin = mxDuplicateArray(prhs[1]);


	mxArray* argout = mxCreateCellMatrix(1, nlhs);

	int code = 0;

	try {

		if(strcmp(fname, "fba") == 0) {
			code = fba(argout, argin);
		}
		else if(strcmp(fname, "tfba") == 0) {
			code = tfba(argout, argin);
		}
		else if(strcmp(fname, "fva") == 0) {
			code = fva(argout, argin);
		}
		else if(strcmp(fname, "tfva") == 0) {
			code = tfva(argout, argin);
		}
		else {
			code = help(argout, argin);
		}

	}
	catch(std::exception& ex) {
		mexPrintf ("crashed because of exception\n");
		mexPrintf (diagnostic_information(ex).c_str());
		mexPrintf ("\n");
		return;
	}

	if(code != 0) {
		return;
	}

	plhs[0] = argout;

	return;
}
