/*
 * SBMLLoader.h
 *
 *  Created on: 03.04.2012
 *      Author: arnem
 */

#ifndef SBMLLOADER_H_
#define SBMLLOADER_H_

#include "metaopt/Properties.h"
#include "metaopt/model/Model.h"
#include "metaopt/model/Metabolite.h"
#include "sbml/SBMLTypes.h"
#include <boost/unordered_map.hpp>
#include "metaopt/model/impl/FullModel.h"
#include "metaopt/Uncopyable.h"

namespace metaopt {

class SBMLLoader : Uncopyable {
public:
	SBMLLoader();
	virtual ~SBMLLoader();

	/** load the model */
	void load(const libsbml::Model* m);

	/** fetch the loaded model */
	ModelPtr getModel() const;

	/** clear loader for loading another model */
	void clear();

private:
	/** temporary map of metabolite names to created objects */
	boost::unordered_map<std::string, MetabolitePtr> _metabolites;
	boost::shared_ptr<FullModel> _model;

	void buildReaction(const libsbml::Reaction* r);
};

ModelPtr loadSBMLModel(const libsbml::Model* m);


}

#endif /* SBMLLOADER_H_ */
