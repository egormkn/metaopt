/*
 * ISSupply.h
 *
 *  Created on: 10.07.2012
 *      Author: arne
 */

#ifndef ISSUPPLY_H_
#define ISSUPPLY_H_

#include <boost/shared_ptr.hpp>
#include "metaopt/model/Reaction.h"

namespace metaopt {

class ISSupply {
public:
	virtual ~ISSupply();

	/**
	 * For rxn in the computed infeasible set, gives the coefficient.
	 */
	virtual double getAlpha(ReactionPtr rxn) = 0;
};

typedef boost::shared_ptr<ISSupply> ISSupplyPtr;

} /* namespace metaopt */
#endif /* ISSUPPLY_H_ */
