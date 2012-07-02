/*
 * MatlabLoader.cpp
 *
 *  Created on: 02.07.2012
 *      Author: arnem
 */

#include "MatlabLoader.h"
#include <boost/throw_exception.hpp>
#include <boost/exception/all.hpp>
#include <iostream>

#include "metaopt/model/Reaction.h"
#include "metaopt/model/Metabolite.h"
#include "metaopt/model/impl/FullModel.h"

using namespace boost;
using namespace std;

#define F_S "S"
#define F_UB "ub"
#define F_LB "lb"
#define F_C "c"
#define F_RXNS "rxns"
#define F_METS "mets"
#define F_BOUNDARY "boundary"
#define F_POT_UB "p_ub"
#define F_POT_LB "p_lb"
#define F_POT_C "p_c"
#define F_INT_RXNS "int_rxns"
#define F_INT_METS "int_mets"

#define ERROR_MSG(field,msg) BOOST_THROW_EXCEPTION( MatlabLoaderError() << field_name(field) << matlab_error_message(msg))

#define STRICT

namespace metaopt {

MatlabLoader::MatlabLoader() {
	// TODO Auto-generated constructor stub

}

MatlabLoader::~MatlabLoader() {
	// TODO Auto-generated destructor stub
}

void MatlabLoader::load(const mxArray* m) {
	_model = shared_ptr<FullModel>(new FullModel());
	// model given by stoichiometric matrix etc.
	_S = mxGetField(m, 0, F_S);
	_ub = mxGetField(m, 0, F_UB);
	_lb = mxGetField(m, 0, F_LB);
	_c = mxGetField(m, 0, F_C);
	_rxns = mxGetField(m, 0, F_RXNS);
	_mets = mxGetField(m, 0, F_METS);
	_boundary = mxGetField(m, 0, F_BOUNDARY);
	_pot_ub = mxGetField(m,0, F_POT_UB);
	_pot_lb = mxGetField(m,0, F_POT_LB);
	_pot_c = mxGetField(m,0, F_POT_C);
	_int_rxns = mxGetField(m,0, F_INT_RXNS);
	_int_mets = mxGetField(m,0, F_INT_METS);

	// test if data is has proper format
	testFields();
	testTypes();
	testDimensions();

	// create metabolites
	int num_species = mxGetM(_S); // rows = # species
	int num_reactions = mxGetN(_S); // cols = # reactions
	for(int i = 0; i < num_species; i++) {
		mxArray* name = mxGetCell(_mets, i);
		if(!mxIsChar(name)) ERROR_MSG(F_METS,"Metabolite names must be strings");
		char* cname = mxArrayToString(name);
		assert(cname != NULL);
		string sname(cname);
		mxFree(cname);
		_metabolites.push_back(_model->createMetabolite(sname));
	}

	// set boundary condition
	if(_boundary != NULL) {
		double* boundData = mxGetPr(_boundary);
		int num = mxGetNumberOfElements(_boundary);
		for(int i = 0; i < num; i++) {
			int j = (int) boundData[i];
			if(0 <= j && j < num_species) {
				_metabolites[j]->setBoundaryCondition(true);
			}
			else {
				cout << "WARNING: Boundary condition index " << j << " out of bounds (no such metabolite)" << endl;
			}
		}
	}

	// create reactions
	double* data = mxGetPr(_S);
	mwIndex* rows = mxGetIr(_S);
	mwIndex* pos = mxGetJc(_S);

	double* ubData = mxGetPr(_ub);
	double* lbData = mxGetPr(_lb);
	double* cData = mxGetPr(_c);

	for(int i = 0; i < num_reactions; i++) {
		// create reaction (and read the name)
		mxArray* name = mxGetCell(_rxns, i);
		if(!mxIsChar(name)) ERROR_MSG(F_RXNS,"The metabolite names must be strings");
		char* cname = mxArrayToString(name);
		assert(cname != NULL);
		string sname(cname);
		mxFree(cname);
		ReactionPtr reaction = _model->createReaction(sname);

		// initialize stoichiometries
		for(mwIndex j = pos[i]; j < pos[i+1]; j++) {
			//cout << i << " " << rows[j] << " " << data[j] << endl;
			MetabolitePtr met = _metabolites.at(rows[j]);
			double stoich = data[j];
			reaction->setStoichiometry(met, stoich);
		}
		reaction->setLb(lbData[i]);
		reaction->setUb(ubData[i]);
		reaction->setObj(cData[i]);
		reaction->setExchange(true); // default is exchange, may be switched back later
	}

	// set internal reactions
	if(_int_rxns != NULL) {
		// use the explicit information given
		int num_internal_rxns = mxGetNumberOfElements(_int_rxns);
		for(int i = 0; i < num_internal_rxns; i++) {
			mxArray* name = mxGetCell(_int_rxns, i);
			if(!mxIsChar(name)) ERROR_MSG(F_METS,"Internal reaction names must be strings");
			char* cname = mxArrayToString(name);
			assert(cname != NULL);
			string sname(cname);
			mxFree(cname);
			if(!_model->hasReaction(sname)) ERROR_MSG(F_INT_RXNS,"There exists no reaction with name " + sname);
			ReactionPtr rxn = _model->getReaction(sname); // not efficient, but this is only init.
			rxn->setExchange(false);
		}
	}
	else {
		// use the implicit rule that every reaction with only one metabolite is an exchange reaction.
		foreach(ReactionPtr rxn, _model->getReactions()) {
			rxn->setExchange(rxn->getStoichiometries().size() <= 1);
		}
	}

	// set potential bounds
	double* pubData = NULL;
	double* plbData = NULL;
	double* pcData = NULL;

	if(_pot_ub != NULL) pubData = mxGetPr(_pot_ub);
	if(_pot_lb != NULL) plbData = mxGetPr(_pot_lb);
	if(_pot_c != NULL) pcData = mxGetPr(_pot_c);

	int num_int_mets = mxGetNumberOfElements(_int_mets);
	for(int i = 0; i < num_int_mets; i++) {
		mxArray* name = mxGetCell(_int_mets, i);
		if(!mxIsChar(name)) ERROR_MSG(F_METS,"Internal metabolite names must be strings");
		char* cname = mxArrayToString(name);
		assert(cname != NULL);
		string sname(cname);
		mxFree(cname);
		if(!_model->hasMetabolite(sname)) ERROR_MSG(F_INT_METS,"There exists no metabolite with name " + sname);
		MetabolitePtr met = _model->getMetabolite(sname);
		if(pubData != NULL) met->setPotUb(pubData[i]);
		if(plbData != NULL) met->setPotLb(plbData[i]);
		if(pcData != NULL) met->setPotObj(pcData[i]);
	}
}

void MatlabLoader::testFields() {
	// check, that all the necessary fields are specified
	if(_S == NULL) ERROR_MSG(F_S, "The model must contain a stoichiometric matrix S");
	if(_ub == NULL) ERROR_MSG(F_UB, "The model must contain flux upper bounds ub");
	if(_lb == NULL) ERROR_MSG(F_LB, "The model must contain flux lower bounds lb");
	if(_c == NULL) ERROR_MSG(F_C, "The model must contain an objective function c");
	if(_rxns == NULL) ERROR_MSG(F_RXNS, "The model must contain the names of reactions (rxns)");
	if(_mets == NULL) ERROR_MSG(F_METS, "The model must contain the names of metabolites (mets)");
}

void MatlabLoader::testTypes() {
	// check, if all fields have correct type
	if(mxIsChar(_S) || mxIsComplex(_S) || !mxIsNumeric(_S) || !mxIsSparse(_S)) ERROR_MSG(F_S, "The stoichiometric matrix must be sparse and contain reals");
	if(mxIsChar(_ub) || mxIsComplex(_ub) || !mxIsNumeric(_ub) || mxIsSparse(_ub)) ERROR_MSG(F_UB, "The upper bounds must be full and contain reals");
	if(mxIsChar(_lb) || mxIsComplex(_lb) || !mxIsNumeric(_lb) || mxIsSparse(_lb)) ERROR_MSG(F_LB, "The lower bounds matrix must be full and contain reals");
	if(mxIsChar(_c) || mxIsComplex(_c) || !mxIsNumeric(_c) || mxIsSparse(_c)) ERROR_MSG(F_C, "The objective vector must be full and contain reals");
	if(!mxIsCell(_rxns)) ERROR_MSG(F_RXNS, "The list of reactions must be given as a cell array of strings");
	if(!mxIsCell(_mets)) ERROR_MSG(F_METS, "The list of metabolites must be given as a cell array of strings");
	// if boundary is specified, check if it has correct type
	if(_boundary != NULL && (!mxIsNumeric(_boundary) || mxIsComplex(_boundary) || mxIsDouble(_boundary))) {
#ifdef STRICT
		ERROR_MSG(F_BOUNDARY, "A list of boundary metabolites was specified, but not of correct type. Must be a list of indices of metabolites.");
#else
		cout << "A list of boundary metabolites was specified, but not of correct type. Must be a list of indices of metabolites. Ignoring boundary metabolites." << endl;
		_boundary = NULL;
#endif
	}
	if(_pot_ub != NULL && (mxIsChar(_pot_ub) || mxIsComplex(_pot_ub) || !mxIsNumeric(_pot_ub) || mxIsSparse(_pot_ub))) {
#ifdef STRICT
		ERROR_MSG(F_POT_UB, "The upper bounds must be full and contain reals");
#else
		cout << "The upper bounds must be full and contain reals. Ignoring potential upper bounds." << endl;
		_pot_ub = NULL;
#endif
	}
	if(_pot_lb != NULL && (mxIsChar(_pot_lb) || mxIsComplex(_pot_lb) || !mxIsNumeric(_pot_lb) || mxIsSparse(_pot_lb))) {
#ifdef STRICT
		ERROR_MSG(F_POT_LB, "The lower bounds must be full and contain reals");
#else
		cout << "The lower bounds must be full and contain reals. Ignoring potential lower bounds." << endl;
		_pot_lb = NULL;
#endif
	}
	if(_pot_c != NULL && (mxIsChar(_pot_c) || mxIsComplex(_pot_c) || !mxIsNumeric(_pot_c) || mxIsSparse(_pot_c))) {
#ifdef STRICT
		ERROR_MSG(F_POT_C, "The objective vector must be full and contain reals");
#else
		cout << "The objective vector must be full and contain reals. Ignoring potential objectives." << endl;
		_pot_c = NULL;
#endif
	}
	if(_int_rxns != NULL && (!mxIsCell(_int_rxns))) {
#ifdef STRICT
		ERROR_MSG(F_INT_RXNS, "The list of internal reactions must be given as a cell array of strings");
#else
		cout << "The list of internal reactions must be given as a cell array of strings. Ignoring specified internal reactions." << endl;
		_int_rxns = NULL;
#endif
	}
	if(_int_mets != NULL && (!mxIsCell(_int_mets))) {
#ifdef STRICT
		ERROR_MSG(F_INT_RXNS, "The list of internal metabolites must be given as a cell array of strings");
#else
		cout << "The list of internal metabolites must be given as a cell array of strings. Assuming all metabolites to be internal." << endl;
		_int_mets = _mets;
#endif
	}
	else if(_int_mets == NULL) {
		_int_mets = _mets;
	}
}

void MatlabLoader::testDimensions() {
	unsigned int num_species = mxGetM(_S); // rows = # species
	unsigned int num_reactions = mxGetN(_S); // cols = # reactions

	// check dimensions
	if(mxGetNumberOfElements(_ub) != num_reactions) ERROR_MSG(F_UB,"You must give exactly as many upper bounds as the stoichiometric matrix has columns");
	if(mxGetNumberOfElements(_lb) != num_reactions) ERROR_MSG(F_LB,"You must give exactly as many lower bounds as the stoichiometric matrix has columns");
	if(mxGetNumberOfElements(_c) != num_reactions) ERROR_MSG(F_C,"The objective vector must have as many elements as the stoichiometric matrix has columns");
	if(mxGetNumberOfElements(_rxns) != num_reactions) ERROR_MSG(F_RXNS,"You must give a name for every reaction and no more");
	if(mxGetNumberOfElements(_mets) != num_species) ERROR_MSG(F_METS,"You must give a name for every metabolite and no more");

	unsigned int num_int_species = mxGetNumberOfElements(_int_mets);
	if(_pot_ub != NULL && mxGetNumberOfElements(_pot_ub) != num_int_species) ERROR_MSG(F_POT_UB, "You must give exactly as many upper bounds on potentials as there are (internal) metabolites");
	if(_pot_c != NULL && mxGetNumberOfElements(_pot_lb) != num_int_species) ERROR_MSG(F_POT_LB, "You must give exactly as many lower bounds on potentials as there are (internal) metabolites");
	if(_pot_lb != NULL && mxGetNumberOfElements(_pot_c) != num_int_species) ERROR_MSG(F_POT_C, "You must give exactly as many objective coefficients on potentials as there are (internal) metabolites");
}

ModelPtr MatlabLoader::getModel() const {
  return _model;
}

void MatlabLoader::clear() {
	// clear loaded model
	_model.reset();
	_metabolites.clear();
}

ModelPtr loadMatlabModel(const mxArray* m) {
	MatlabLoader loader;
	loader.load(m);
	ModelPtr model = loader.getModel();
	return model;
}

} /* namespace metaopt */
