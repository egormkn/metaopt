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
 * ThermoInfeasibleSetPool.h
 *
 *  Created on: 12.11.2012
 *      Author: arnem
 */

#ifndef THERMOINFEASIBLESETPOOL_H_
#define THERMOINFEASIBLESETPOOL_H_

#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/model/Coupling.h"
#include "metaopt/model/DirectedReaction.h"
#include "metaopt/model/scip/ISSupply.h"

#include "metaopt/Properties.h"

namespace metaopt {

/**
 * Since we do not want to recompute the infeasible sets every time, store the good ones.
 */
struct ThermoInfeasibleSet {
	boost::unordered_set<DirectedReaction> set;
	double priority;

	ThermoInfeasibleSet() : priority(0) {}
};

typedef boost::shared_ptr<ThermoInfeasibleSet> ThermoInfeasibleSetPtr;


class ThermoInfeasibleSetPool {
public:
	ThermoInfeasibleSetPool();
	virtual ~ThermoInfeasibleSetPool();

	void add(ThermoInfeasibleSetPtr set);

	/**
	 * returns a list of the best (smallest priority) infeasible sets.
	 *
	 * Attention: Returns a reference to its internally stored list. Do not modify!
	 */
	std::vector<ThermoInfeasibleSetPtr>&  getInfeasibleSets();

private:
	// we want a data structure with which we can easily just iterate through the best elements,
	// and also update elements easily

	// for now: simply maintain a list of a fixed number of elements, and only store the best ones (smallest priority value)
	std::vector<ThermoInfeasibleSetPtr> infeasibleSets;
	// value of the worst stored item
	double worst;

};

} /* namespace metaopt */
#endif /* THERMOINFEASIBLESETPOOL_H_ */
