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

	Precision(double tol, double slave_factor, double check_factor);

	virtual ~Precision();

	inline double getPrimalFeasTol();

	inline double getDualFeasTol();

	inline double getCheckTol();

	boost::shared_ptr<Precision> getSlavePrecision();

private:
	double _tol;
	double _slave_factor;
	double _check_factor;
	double _check_tol;
};

typedef boost::shared_ptr<Precision> PrecisionPtr;

inline double Precision::getPrimalFeasTol() {
	return _tol;
}

inline double Precision::getDualFeasTol() {
	return _tol;
}

inline double Precision::getCheckTol() {
	return _check_tol;
}



} /* namespace metaopt */
#endif /* PRECISION_H_ */
