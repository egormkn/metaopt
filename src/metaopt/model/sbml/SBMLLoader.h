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

/**
 * The SBMLLoader can be used to load networks in form of SBML files.
 */
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
