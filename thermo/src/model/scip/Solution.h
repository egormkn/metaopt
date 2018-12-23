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
 * Solution.h
 *
 * Defines a boost::shared_ptr wrapper for SCIP_SOL* Objekte
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_

#include "scip/scip.h"
#include <boost/shared_ptr.hpp>

#include "Uncopyable.h"
#include "Properties.h"

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
