/*
 * ModelAddOn.cpp
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#include "ModelAddOn.h"

namespace metaopt {

ModelAddOn::ModelAddOn(ScipModelPtr model) {
	_model = model;
}

ModelAddOn::~ModelAddOn() {
	// nothing to do
}

} /* namespace metaopt */
