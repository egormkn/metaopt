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
 * SteadyStateConstraint.h
 *
 *  Created on: 11.04.2012
 *      Author: arnem
 */

#ifndef STEADYSTATECONSTRAINT_H_
#define STEADYSTATECONSTRAINT_H_

#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Creates linear constraints that enforce the steady state assumption on the flux variables.
 */
void createSteadyStateConstraint(ScipModelPtr smodel);

} /* namespace metaopt */
#endif /* STEADYSTATECONSTRAINT_H_ */
