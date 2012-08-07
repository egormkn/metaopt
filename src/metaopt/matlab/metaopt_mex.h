/*
 * metaopt_mex.h
 *
 *  Created on: 07.08.2012
 *      Author: arnem
 */

#ifndef METAOPT_MEX_H_
#define METAOPT_MEX_H_

#include "mex.h"

/**
 * General interface function to matlab.
 * Every matlab call of the metaopt lib calls this function.
 *
 * The first parameter of the matlab function call must be a string.
 * This string will specify which method will actually be executed.
 *
 * @param nrhs must be at least one
 * @param prhs must contain at least one entry
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* METAOPT_MEX_H_ */
