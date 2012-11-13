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
 * ThermoInfeasibleSetPool.cpp
 *
 *  Created on: 12.11.2012
 *      Author: arnem
 */

#include "ThermoInfeasibleSetPool.h"

#define NUMBER_INFEASIBLE_SETS 10


namespace metaopt {

ThermoInfeasibleSetPool::ThermoInfeasibleSetPool() {
	worst = INFINITY;
}

ThermoInfeasibleSetPool::~ThermoInfeasibleSetPool() {
	// nothing to do
}

void ThermoInfeasibleSetPool::add(ThermoInfeasibleSetPtr tis) {
	if(infeasibleSets.size() < NUMBER_INFEASIBLE_SETS || worst >= tis->priority) {
		// we add the new element
		ThermoInfeasibleSetPtr toInsert = tis;
		for(int i = 0; i < infeasibleSets.size(); i++) {
			if(infeasibleSets[i]->priority > toInsert->priority) {
				ThermoInfeasibleSetPtr old = infeasibleSets[i];
				infeasibleSets[i] = toInsert;
				toInsert = old;
				if(toInsert->set == tis->set) {
					// we don't need the same element twice, so just let it drop out
					// the check is performed here, to allow an improving update of priority
					return;
				}
			}
		}
		if(infeasibleSets.size() < NUMBER_INFEASIBLE_SETS) {
			infeasibleSets.push_back(toInsert);
		}
		worst = infeasibleSets[infeasibleSets.size()-1]->priority;
	}

}

std::vector<ThermoInfeasibleSetPtr>&  ThermoInfeasibleSetPool::getInfeasibleSets() {
	return infeasibleSets;
}

} /* namespace metaopt */
