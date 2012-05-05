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

#include "metaopt/Uncopyable.h"

namespace metaopt {

class ScipModel;
typedef boost::shared_ptr<SCIP_SOL> SolutionPtr;

struct SolutionDestructor {
	boost::shared_ptr<ScipModel> _scip; // storing an instance of the ScipModel makes sure that the scip pointer stays live

	SolutionDestructor(boost::shared_ptr<ScipModel> scip);

	SolutionDestructor(const SolutionDestructor& d);

	void operator()(SCIP_SOL* sol);
};

/**
 * Use for newly allocated solutions where the user needs to manually free the solution.
 */
inline SolutionPtr wrap(SCIP_SOL* sol, boost::shared_ptr<ScipModel> scip) {
	return SolutionPtr(sol, SolutionDestructor(scip));
}

inline void noOpSolutionDestructor(SCIP_SOL* sol) {} // does nothing

/**
 * Use for solutions that are supplied by scip and that must not be deleted manually.
 */
inline SolutionPtr wrap_weak(SCIP_SOL* sol) {

	return SolutionPtr(sol, noOpSolutionDestructor);
}

}

#endif /* SOLUTION_H_ */
