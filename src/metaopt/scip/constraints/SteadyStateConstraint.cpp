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
 * SteadyStateConstraint.cpp
 *
 *  Created on: 11.04.2012
 *      Author: arnem
 */

#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

#include "SteadyStateConstraint.h"
#include "metaopt/model/Metabolite.h"
#include "metaopt/model/Reaction.h"
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "metaopt/Properties.h"
#include "metaopt/scip/ScipError.h"

using namespace boost;
using namespace std;

namespace metaopt {

void createSteadyStateConstraint(ScipModelPtr smodel) {
	ModelPtr m = smodel->getModel();

	// simply create a bunch of linear constraints
	// we have to satisfy flow conservation for every species
	unordered_map<MetabolitePtr,SCIP_CONS*> met_cons;

	// create empty linear constraints for every species and then add the vars
	foreach(MetabolitePtr met, m->getMetabolites()) {
		// only create constraints for internal metabolites,
		// since all other metabolites do not need to satisfy steady state
		if(!met->hasBoundaryCondition()) {
			SCIP_CONS* cons;

			double lhs = 0;
			double rhs = 0;

			BOOST_SCIP_CALL( SCIPcreateConsLinear(smodel->getScip(),
					&cons,
					met->getCName(),   // name
					0,      // num vars
					NULL,   // no vars
					NULL,   // no coefficients
					lhs,    // lhs
					rhs,    // rhs
					true,   // initial
					true,   // separate
					true,   // enforce
					true,   // check
					true,   // propagate
					false, false, false, false, false ));

			BOOST_SCIP_CALL( SCIPaddCons(smodel->getScip(), cons));
			met_cons[met] = cons;
		}
	}

	foreach(ReactionPtr rxn, m->getReactions()) {

		SCIP_VAR* var = smodel->getFlux(rxn);
		foreach(Stoichiometry s, rxn->getStoichiometries()) {
			if(!s.first->hasBoundaryCondition()) {
				BOOST_SCIP_CALL( SCIPaddCoefLinear(smodel->getScip(), met_cons.at(s.first), var , s.second));
			}
		}
	}

	// free species constraints (flow conservation)
	typedef pair<MetabolitePtr,SCIP_CONS*> MetCons;
	foreach(MetCons met, met_cons) {
		SCIP_CONS* con = met.second;
		BOOST_SCIP_CALL( SCIPreleaseCons(smodel->getScip(), &con) );
	}
}

} /* namespace metaopt */
