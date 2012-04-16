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
				BOOST_SCIP_CALL( SCIPaddCoefLinear(smodel->getScip(), met_cons[s.first], var , s.second));
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
