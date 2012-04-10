/*
 * LPFlux.cpp
 *
 *  Created on: 10.04.2012
 *      Author: arnem
 */

#include "LPFlux.h"
#include "objscip/objscip.h"
#include "metaopt/Properties.h"

namespace metaopt {

LPFlux::LPFlux(ModelPtr model, bool exchange) {
	_model = model;
	init_lp(exchange);
}

SCIP_RETCODE LPFlux::init_lp(bool exchange) {
	_lpi = NULL;
	SCIP_CALL( SCIPlpiCreate(&_lpi, "LPFlux", SCIP_OBJSEN_MAXIMIZE) );

	foreach(ReactionPtr r, _model->getReactions()) {
		if(exchange || !r->isExchange()) {
			SCIPlpiAddRows(_lpi, )
		}
	}
}

LPFlux::~LPFlux() {
	// TODO Auto-generated destructor stub
}

} /* namespace metaopt */
