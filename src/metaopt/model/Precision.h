/*
 * Precision.h
 *
 *  Created on: 27.11.2012
 *      Author: arnem
 */

#ifndef PRECISION_H_
#define PRECISION_H_

#include <boost/shared_ptr.hpp>

namespace metaopt {

class Precision {
public:
	Precision(double tol);

	Precision(double primal_tol, double dual_tol);

	Precision(double primal_tol, double dual_tol, double slave_factor, double check_factor);

	virtual ~Precision();

	inline double getPrimalFeasTol();

	inline double getDualFeasTol();

	inline double getCheckTol();

	/**
	 * Gets precision with tighter primal feasibility tolerance.
	 * Use for computing primal feasible solutions.
	 */
	boost::shared_ptr<Precision> getPrimalSlavePrecision();

	/**
	 * Gets precision with tighter dual feasibility tolerance.
	 * Use for computing dual feasible solutions.
	 */
	boost::shared_ptr<Precision> getDualSlavePrecision();

private:
	double _primal_tol;
	double _dual_tol;
	double _slave_factor;
	double _check_factor;
	double _check_tol;
};

typedef boost::shared_ptr<Precision> PrecisionPtr;

inline double Precision::getPrimalFeasTol() {
	return _primal_tol;
}

inline double Precision::getDualFeasTol() {
	return _dual_tol;
}

inline double Precision::getCheckTol() {
	return _check_tol;
}



} /* namespace metaopt */
#endif /* PRECISION_H_ */
