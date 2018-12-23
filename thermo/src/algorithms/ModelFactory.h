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
 * ModelFactory.h
 *
 *  Created on: 23.07.2012
 *      Author: arnem
 */

#ifndef MODELFACTORY_H_
#define MODELFACTORY_H_

#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/Properties.h"

namespace metaopt {

class ModelFactory {
public:
	ModelFactory();
	virtual ~ModelFactory();

	/**
	 * build a brand new ScipModel.
	 * Implement to build ScipModel with your desired features.
	 * @model the model to build the ScipModel from.
	 */
	virtual ScipModelPtr build(ModelPtr model) = 0;
};

} /* namespace metaopt */
#endif /* MODELFACTORY_H_ */
