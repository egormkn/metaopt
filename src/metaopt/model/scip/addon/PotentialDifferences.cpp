/*
 * PotentialDifferences.cpp
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "PotentialDifferences.h"
#include "metaopt/scip/ScipError.h"

using namespace boost;

namespace metaopt {

PotentialDifferences::PotentialDifferences(ScipModelPtr model) : ModelAddOn(model, "potential differences") {
	// nothing to do
}

PotentialDifferences::~PotentialDifferences() {
	// already happened in destroy
}

void PotentialDifferences::destroy(ScipModel* model) {
	// free vars
	typedef std::pair<ReactionPtr, SCIP_VAR*> PotDiffVar;
	foreach(PotDiffVar r, _potDiff) {
		SCIP_VAR* var = r.second;
		int code = SCIPreleaseVar(model->getScip(), &var);
		assert(code == SCIP_OKAY);
	}
	_potDiff.clear(); // make sure there are no invalid pointers hanging around
	ModelAddOn::destroy(model);
}

SCIP_VAR* PotentialDifferences::getPotDiff(ReactionPtr rxn) {
	assert(!isDestroyed());
	unordered_map<ReactionPtr, SCIP_VAR*>::iterator iter = _potDiff.find(rxn);
	if(iter != _potDiff.end()) {
		return iter->second;
	}
	else {
		// Create Variable
		ScipModelPtr model = getModel(); // make sure the model stays living
		SCIP* scip = model->getScip();

		std::string name = rxn->getName() + "_potDiff";

		SCIP_VAR* var = NULL;
		BOOST_SCIP_CALL( SCIPcreateVar(scip,
				&var,
				name.c_str(),
				-SCIPinfinity(scip),
				SCIPinfinity(scip),
				0,
				SCIP_VARTYPE_CONTINUOUS,
				true,
				false,
				0,0,0,0,0) );

		BOOST_SCIP_CALL( SCIPaddVar(scip, var) );

		BOOST_SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, var) ); //TODO: Is this really necessary?

		// Now we also have to create the constraint that defines the variables value
		int nvars = rxn->getStoichiometries().size() + 1; // do not forget the potential difference var!
		SCIP_VAR* vars[nvars];
		double vals[nvars];
		// add potential difference var
		vars[0] = var;
		vals[0] = -1;
		int i = 1;
		foreach(Stoichiometry s, rxn->getStoichiometries()) {
			vars[i] = model->getPotential(s.first);
			vals[i] = s.second;
			i++;
		}
		assert(i == nvars);

		SCIP_CONS* cons = NULL;
		BOOST_SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name.c_str(), nvars, vars, vals, 0,0, true, true, true, true, true, false, false, false, false, false));
		BOOST_SCIP_CALL( SCIPaddCons(scip, cons));
		BOOST_SCIP_CALL( SCIPreleaseCons(scip, &cons));

		_potDiff[rxn] = var;
		return var;
	}
}

bool PotentialDifferences::computeSolutionVals(SolutionPtr sol) const {
	assert(!isDestroyed());
	typedef std::pair<ReactionPtr, SCIP_VAR*> PotDiffVar;

	ScipModelPtr model = getModel();

	foreach(PotDiffVar v, _potDiff) {
		ReactionPtr rxn = v.first;
		double val = 0;
		foreach(Stoichiometry s, rxn->getStoichiometries()) {
			SCIP_VAR* met = model->getPotential(s.first);
			val += s.second * SCIPgetSolVal(model->getScip(), sol.get(), met);
		}
		BOOST_SCIP_CALL( SCIPsetSolVal(model->getScip(), sol.get(), v.second, val) );
	}
	return true;
}

double PotentialDifferences::getCurrentPotDiff(ReactionPtr rxn) {
	assert(!isDestroyed());
	ScipModelPtr model = getModel();
	SCIP* scip = model->getScip();
	assert( hasPotDiff(rxn) );
	if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING) {
		assert( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL); // LP not solved to optimality
		return SCIPgetSolVal(scip, NULL, getPotDiff(rxn));
	}
	else if(SCIPgetStage(scip) == SCIP_STAGE_SOLVED) {
		SCIP_SOL* sol = SCIPgetBestSol(scip);
		return SCIPgetSolVal(scip, sol, getPotDiff(rxn));
	}
	else {
		assert( false ); // Solving has not yet started!
		BOOST_THROW_EXCEPTION( PreconditionViolatedException() << var_state("Solving has not yet started!") );
		return 0; // never called
	}
}

bool PotentialDifferences::hasPotDiff(ReactionPtr rxn) {
	assert(!isDestroyed());
	return ( _potDiff.find(rxn) != _potDiff.end() );
}

} /* namespace metaopt */
