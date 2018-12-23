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
