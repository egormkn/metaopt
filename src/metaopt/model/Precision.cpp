/*
 * Precision.cpp
 *
 *  Created on: 27.11.2012
 *      Author: arnem
 */

#include "Precision.h"

namespace metaopt {

// check factor is bigger than 1, since it should be tolerant to accumulating errors
#define CHECK_FACTOR_DEFAULT 1e1
// slave factor is less than 1, since the slave should solve with more precision, so that accumulation of errors is not fatal for the master
#define SLAVE_FACTOR_DEFAULT 1e-2

Precision::Precision(double tol) {
	_tol = tol;
	_slave_factor = SLAVE_FACTOR_DEFAULT;
	_check_factor = CHECK_FACTOR_DEFAULT;
	_check_tol = _tol*_check_factor;
}

Precision::Precision(double tol, double slave_factor, double check_factor) {
	_tol = tol;
	_slave_factor = slave_factor;
	_check_factor = check_factor;
	_check_tol = _tol*_check_factor;
}

Precision::~Precision() {
	// nothing to destruct
}

PrecisionPtr Precision::getSlavePrecision() {
	return PrecisionPtr(new Precision(_slave_factor*_tol, _slave_factor, _check_factor));
}

} /* namespace metaopt */
