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
 * OctaveLoader.h
 *
 *  Created on: 02.07.2012
 *      Author: arnem
 */

#ifndef OCTAVELOADER_H_
#define OCTAVELOADER_H_

#include <string>
#include <boost/exception/all.hpp>

#include "Properties.h"
#include "model/Model.h"
#include "model/Metabolite.h"
#include <boost/unordered_map.hpp>
#include "model/impl/FullModel.h"
#include "Uncopyable.h"
#include "mex.h"


namespace metaopt {

struct OctaveLoaderError  : virtual boost::exception, virtual std::exception {
public:
	OctaveLoaderError() {
		//assert(false);
	}
};

class OctaveLoader : Uncopyable {
public:
	OctaveLoader();
	virtual ~OctaveLoader();

	/** load the model */
	void load(const mxArray* m);

	/** fetch the loaded model */
	ModelPtr getModel() const;

	/**
	 * fetches the metabolite with the specified index.
	 * Requires that a model has been loaded.
	 */
	MetabolitePtr getMetabolite(int index);

	/**
	 * fetches the reaction with the specified index.
	 * Requires that a model has been loaded.
	 */
	ReactionPtr getReaction(int index);

	/**
	 * loads coupling information for the model
	 */
	//CouplingPtr loadCoupling(mxArray* fctable, mxArray* blocked);

	/** clear loader for loading another model */
	void clear();

private:
	boost::shared_ptr<FullModel> _model;

	mxArray* _S; 		///< stoichiometric matrix
	mxArray* _lb; 		///< flux lower bounds
	mxArray* _ub;		///< flux upper bounds
	mxArray* _c;			///< flux cost function
	mxArray* _rxns;		///< reaction names
	mxArray* _mets;		///< metabolite names
	mxArray* _boundary;	///< indicator if boundary metabolite (i.e. if met has to satisfy mass-balance)
	mxArray* _pot_ub;	///< upper bound on potentials (for internal metabolites)
	mxArray* _pot_lb;	///< lower bound on potentials (for internal metabolites)
	mxArray* _pot_c;		///< cost function on potentials (for internal metabolites)
	mxArray* _int_rxns;  ///< list of internal reactions (that have to satisfy thermodynamics)
	mxArray* _int_mets;	///< list of internal metabolites (that have well defined potential)
	mxArray* _precision; /// scalar giving the precision with which the model should be analyzed
	mxArray* _dprecision; /// scalar giving the dual precision with which the model should be analyzed

	/** map of metabolite indices to created objects */
	std::vector<MetabolitePtr> _metabolites;
	/** map of reaction indices to created objects */
	std::vector<ReactionPtr> _reactions;

	/**
	 * test if all necessary fields are supplied.
	 * If not, a OctaveLoaderError is thrown.
	 */
	void testFields();

	/**
	 * test if all fields have their proper type.
	 * If not, a OctaveLoaderError is thrown.
	 */
	void testTypes();

	/**
	 * test if the dimensions of all specified fields match.
	 * If not, a OctaveLoaderError is thrown.
	 */
	void testDimensions();


};

/** If something went wrong with a certain field by loading the model, this gives the name of the field */
typedef boost::error_info<struct tag_field_name,std::string> field_name;
typedef boost::error_info<struct tag_message,std::string> octave_error_message;



ModelPtr loadOctaveModel(const mxArray* m);


} /* namespace metaopt */
#endif /* OCTAVELOADER_H_ */
