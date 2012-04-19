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

namespace metaopt {

class ScipModel;
typedef boost::shared_ptr<SCIP_SOL> SolutionPtr;

struct SolutionDestructor {
	boost::shared_ptr<ScipModel> _scip; // storing an instance of the ScipModel makes sure that the scip pointer stays live

	SolutionDestructor(boost::shared_ptr<ScipModel> scip);

	SolutionDestructor(const SolutionDestructor& d);

	void operator()(SCIP_SOL* sol);
};

inline SolutionPtr wrap(SCIP_SOL* sol, boost::shared_ptr<ScipModel> scip) {
	return SolutionPtr(sol, SolutionDestructor(scip));
}

}

#endif /* SOLUTION_H_ */
