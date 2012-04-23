/*
 * ReactionDirections.cpp
 *
 *  Created on: 19.04.2012
 *      Author: arnem
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "metaopt/scip/ScipError.h"

#include "ReactionDirections.h"

using namespace boost;

namespace metaopt {

// For the strict inequalities, we have to introduce an epsilon for now
#define REACTION_DIRECTIONS_EPSILON 1

ReactionDirections::ReactionDirections(ScipModelPtr model, PotentialDifferencesPtr potDiff)
		: ModelAddOn(model, "reaction directions")
{
	_potDiff = potDiff;
}

ReactionDirections::~ReactionDirections() {
	// already happened in destroy
}

void ReactionDirections::destroy(ScipModel* model) {
	// free vars
	typedef std::pair<ReactionPtr, SCIP_VAR*> DirVar;
	foreach(DirVar r, _dirs) {
		SCIP_VAR* var = r.second;
		int code = SCIPreleaseVar(model->getScip(), &var);
		assert(code == SCIP_OKAY);
	}
	_dirs.clear(); // make sure there are no invalid pointers hanging around
	_potDiff.reset();
	ModelAddOn::destroy(model);
}

SCIP_VAR* ReactionDirections::getDirection(ReactionPtr rxn) {
	assert(!isDestroyed());
	unordered_map<ReactionPtr, SCIP_VAR*>::iterator iter = _dirs.find(rxn);
	if(iter != _dirs.end()) {
		return iter->second;
	}
	else {
		// Create Variable
		ScipModelPtr model = getModel();
		SCIP* scip = model->getScip();

		std::string name = rxn->getName() + "_dir";

		SCIP_VAR* var = NULL;
		BOOST_SCIP_CALL( SCIPcreateVar(scip,
				&var,
				name.c_str(),
				0,
				1,
				0,
				SCIP_VARTYPE_BINARY,
				true,
				false,
				0,0,0,0,0) );

		BOOST_SCIP_CALL( SCIPaddVar(scip, var) );

		BOOST_SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) ); //TODO: Is this really necessary?

		// Now we also have to create the constraints that define the variable's value
		// we will use indicator constraints which enforce ax <= b if z = 1, where z is the binary var.

		SCIP_VAR* neg_var = NULL;
		BOOST_SCIP_CALL( SCIPgetNegatedVar(scip, var, &neg_var) );

		SCIP_VAR* flux_var = model->getFlux(rxn);
		SCIP_VAR* pot_var = _potDiff->getPotDiff(rxn);

		// if var = 1, flux >= 0
		SCIP_CONS* cons = NULL;
		std::string name_flux_fwd = name+"_flux_fwd";
		double val = -1; // coefficient of flux var
		BOOST_SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name_flux_fwd.c_str(), var, 1, &flux_var, &val, 0, true, true, true, true, true, false, false, false, false) );
		BOOST_SCIP_CALL( SCIPaddCons(scip, cons));
		BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons));

		// if var = 0, flux <= 0
		std::string name_flux_bwd = name+"_flux_bwd";
		val = 1; // coefficient of flux var
		BOOST_SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name_flux_bwd.c_str(), neg_var, 1, &flux_var, &val, 0, true, true, true, true, true, false, false, false, false) );
		BOOST_SCIP_CALL( SCIPaddCons(scip, cons));
		BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons));

		// if var = 1, potDiff <= -EPSILON
		std::string name_pot_fwd = name+"_pot_fwd";
		val = 1; // coefficient of pot var
		BOOST_SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name_pot_fwd.c_str(), var, 1, &pot_var, &val, -REACTION_DIRECTIONS_EPSILON, true, true, true, true, true, false, false, false, false) );
		BOOST_SCIP_CALL( SCIPaddCons(scip, cons));
		BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons));

		// if var = 0, potDiff >= EPSILON
		std::string name_pot_bwd = name+"_pot_bwd";
		val = -1; // coefficient of pot var
		BOOST_SCIP_CALL( SCIPcreateConsIndicator(scip, &cons, name_pot_bwd.c_str(), neg_var, 1, &pot_var, &val, -REACTION_DIRECTIONS_EPSILON, true, true, true, true, true, false, false, false, false) );
		BOOST_SCIP_CALL( SCIPaddCons(scip, cons));
		BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons));

		_dirs[rxn] = var;
		return var;
	}
}

bool ReactionDirections::computeSolutionVals(SolutionPtr sol) const {
	assert(!isDestroyed());
	typedef std::pair<ReactionPtr, SCIP_VAR*> DirVar;

	foreach(DirVar v, _dirs) {
		ReactionPtr rxn = v.first;
		ScipModelPtr model = getModel();
		// fetch the value of the potential difference
		double val = SCIPgetSolVal(model->getScip(), sol.get(), _potDiff->getPotDiff(rxn));

		if(val < 0) {
			// if it is negative, set direction to 1 (fwd flux)
			BOOST_SCIP_CALL( SCIPsetSolVal(model->getScip(), sol.get(), v.second, 1) );
		}
		else {
			// set to 0 (bwd flux)
			BOOST_SCIP_CALL( SCIPsetSolVal(model->getScip(), sol.get(), v.second, 0) );
		}
	}
	return true;
}


double ReactionDirections::getCurrentDirection(ReactionPtr rxn) {
	assert(!isDestroyed());
	ScipModelPtr model = getModel();
	SCIP* scip = model->getScip();
	assert( hasDirection(rxn) );
	if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING) {
		assert( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL); // LP not solved to optimality
		return SCIPgetSolVal(scip, NULL, getDirection(rxn));
	}
	else if(SCIPgetStage(scip) == SCIP_STAGE_SOLVED) {
		SCIP_SOL* sol = SCIPgetBestSol(scip);
		return SCIPgetSolVal(scip, sol, getDirection(rxn));
	}
	else {
		assert( false ); // Solving has not yet started!
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("Solving has not yet started!") );
		return 0; // never called
	}
}

bool ReactionDirections::hasDirection(ReactionPtr rxn) {
	assert(!isDestroyed());
	return ( _dirs.find(rxn) != _dirs.end() );
}


} /* namespace metaopt */
