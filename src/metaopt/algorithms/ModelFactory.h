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
