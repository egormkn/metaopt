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
 * WeakNaiveThermoConstraint.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef RELAXEDNAIVETHERMOCONSTRAINT_H_
#define RELAXEDNAIVETHERMOCONSTRAINT_H_

#include "model/scip/addon/PotentialDifferences.h"
#include "model/scip/addon/ReactionDirections.h"
#include "Properties.h"

namespace metaopt {

/**
 * Creates a naive implementation of thermodynamic constraints where all strict inequalities are relaxed to weak inequalities.
 *
 * The relaxed constraint simply states that flux and potential difference must have different signs, or one of them has to be zero.
 *
 * Note, that this formulation is significantly weaker than even the weak formulation of thermodynamic constraints.
 *
 * This constraint requires ReactionDirections as addon.
 * It basically simply initializes the constraints of the reaction direction vars of all internal reactions
 */
void createRelaxedNaiveThermoConstraint(ScipModelPtr model, ReactionDirectionsPtr dirs);

} /* namespace metaopt */
#endif /* RELAXEDNAIVETHERMOCONSTRAINT_H_ */
