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
	_primal_tol = tol;
	_dual_tol = tol;
	_slave_factor = SLAVE_FACTOR_DEFAULT;
	_check_factor = CHECK_FACTOR_DEFAULT;
	_check_tol = _primal_tol*_check_factor;
}

Precision::Precision(double primal_tol, double dual_tol) {
	_primal_tol = primal_tol;
	_dual_tol = dual_tol;
	_slave_factor = SLAVE_FACTOR_DEFAULT;
	_check_factor = CHECK_FACTOR_DEFAULT;
	_check_tol = _primal_tol*_check_factor;
}

Precision::Precision(double primal_tol, double dual_tol, double slave_factor, double check_factor) {
	_primal_tol = primal_tol;
	_dual_tol = dual_tol;
	_slave_factor = slave_factor;
	_check_factor = check_factor;
	_check_tol = _primal_tol*_check_factor;
}

Precision::~Precision() {
	// nothing to destruct
}

PrecisionPtr Precision::getPrimalSlavePrecision() {
	return PrecisionPtr(new Precision(_slave_factor*_primal_tol, _dual_tol, _slave_factor, _check_factor));
}

PrecisionPtr Precision::getDualSlavePrecision() {
	return PrecisionPtr(new Precision(_primal_tol, _slave_factor*_dual_tol, _slave_factor, _check_factor));
}

} /* namespace metaopt */
