/*
 * Solution.h
 *
 * Defines a shared_ptr wrapper for SCIP_SOL* Objekte
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_

#include "scip/scip.h"
#include <boost/shared_ptr.hpp>
#include "metaopt/model/scip/ScipModel.h"

namespace metaopt {

typedef boost::shared_ptr<SCIP_SOL> SolutionPtr;

struct SolutionDestructor {
	ScipModelPtr _scip; // storing an instance of the ScipModel makes sure that the scip pointer stays live

	SolutionDestructor(ScipModelPtr scip) {
		_scip = scip;
	}

	SolutionDestructor(const SolutionDestructor& d) {
		_scip = d._scip;
	}

	void operator()(SCIP_SOL* sol) {
		SCIPfreeSol(_scip->getScip(), &sol);
	}
};

inline SolutionPtr wrap(SCIP_SOL* sol, ScipModelPtr scip) {
	return SolutionPtr(sol, SolutionDestructor(scip));
}

}

#endif /* SOLUTION_H_ */
