/*
 * MatlabLoader.h
 *
 *  Created on: 02.07.2012
 *      Author: arnem
 */

#ifndef MATLABLOADER_H_
#define MATLABLOADER_H_

#include <string>
#include <boost/exception/all.hpp>

#include "metaopt/Properties.h"
#include "metaopt/model/Model.h"
#include "metaopt/model/Metabolite.h"
#include <boost/unordered_map.hpp>
#include "metaopt/model/impl/FullModel.h"
#include "metaopt/Uncopyable.h"
#include "matrix.h"
#include "mat.h"


namespace metaopt {

struct MatlabLoaderError  : virtual boost::exception, virtual std::exception {
public:
	MatlabLoaderError() {
		//assert(false);
	}
};

class MatlabLoader : Uncopyable {
public:
	MatlabLoader();
	virtual ~MatlabLoader();

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

	/** map of metabolite indices to created objects */
	std::vector<MetabolitePtr> _metabolites;
	/** map of reaction indices to created objects */
	std::vector<ReactionPtr> _reactions;

	/**
	 * test if all necessary fields are supplied.
	 * If not, a MatlabLoaderError is thrown.
	 */
	void testFields();

	/**
	 * test if all fields have their proper type.
	 * If not, a MatlabLoaderError is thrown.
	 */
	void testTypes();

	/**
	 * test if the dimensions of all specified fields match.
	 * If not, a MatlabLoaderError is thrown.
	 */
	void testDimensions();


};

/** If something went wrong with a certain field by loading the model, this gives the name of the field */
typedef boost::error_info<struct tag_field_name,std::string> field_name;
typedef boost::error_info<struct tag_message,std::string> matlab_error_message;



ModelPtr loadMatlabModel(const mxArray* m);


} /* namespace metaopt */
#endif /* MATLABLOADER_H_ */
