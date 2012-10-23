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
 * PotSpaceConstraint.h
 *
 *  Created on: 16.07.2012
 *      Author: arnem
 */

#ifndef POTSPACECONSTRAINT_H_
#define POTSPACECONSTRAINT_H_

#include <boost/unordered_map.hpp>
#include "scip/scip.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Use this constraint to formulate side constraints on the potential space variables, even if we do not create the potential variables explicitly.
 * Formulates _coef * mu >= 0
 */
struct PotSpaceConstraint {
	boost::unordered_map<MetabolitePtr, double> _coef;
};

typedef boost::shared_ptr<PotSpaceConstraint> PotSpaceConstraintPtr;


} /* namespace metaopt */
#endif /* POTSPACECONSTRAINT_H_ */
