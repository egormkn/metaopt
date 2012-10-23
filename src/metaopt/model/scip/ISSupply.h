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
 * ISSupply.h
 *
 *  Created on: 10.07.2012
 *      Author: arne
 */

#ifndef ISSUPPLY_H_
#define ISSUPPLY_H_

#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>
#include <vector>
#include "metaopt/model/Reaction.h"
#include "metaopt/model/scip/PotSpaceConstraint.h"
#include "metaopt/Properties.h"

namespace metaopt {

class ISSupply {
public:
	virtual ~ISSupply();

	/**
	 * For rxn in the computed infeasible set, gives the coefficient.
	 */
	virtual double getAlpha(ReactionPtr rxn) = 0;

	/**
	 * updates list of extra constraints
	 */
	virtual void setExtraPotConstraints(boost::unordered_set<PotSpaceConstraintPtr>& psc) = 0;
};

typedef boost::shared_ptr<ISSupply> ISSupplyPtr;

} /* namespace metaopt */
#endif /* ISSUPPLY_H_ */
