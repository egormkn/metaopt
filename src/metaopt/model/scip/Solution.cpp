/*
 * Solution.cpp
 *
 *  Created on: 19.04.2012
 *      Author: arnem
 */

#include "Solution.h"
#include "ScipModel.h"

namespace metaopt {

SolutionDestructor::SolutionDestructor(ScipModelPtr scip) {
	_scip = scip;
}

SolutionDestructor::SolutionDestructor(const SolutionDestructor& d) {
	_scip = d._scip;
}

void SolutionDestructor::operator()(SCIP_SOL* sol) {
	SCIPfreeSol(_scip->getScip(), &sol);
}

}
