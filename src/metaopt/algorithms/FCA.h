/*
 * FCA.h
 *
 *  Created on: 14.07.2012
 *      Author: arne
 */

#ifndef FCA_H_
#define FCA_H_

#include "metaopt/model/Model.h"
#include "metaopt/model/Coupling.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Runs flux coupling analysis and inserts the computed couplings into coupling.
 * This only adds coupled reaction pairs.
 * If you want to use the data, you must first call coupling->computeClosure();
 */
void fca(ModelPtr model, CouplingPtr coupling);


} /* namespace metaopt */
#endif /* FCA_H_ */
