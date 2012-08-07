/*
 * ExitEventHandler.h
 *
 *  Created on: Jan 17, 2012
 *      Author: arne
 */

#ifndef EXITEVENTHANDLER_H_
#define EXITEVENTHANDLER_H_

#include "objscip/objscip.h"

#include "metaopt/model/scip/ScipModel.h"
#include "metaopt/scip/ScipError.h"

using namespace scip;

namespace metaopt {

class ExitEventHandler : public ObjEventhdlr {
public:
	ExitEventHandler(SCIP* scip);
	SCIP_RETCODE scip_exec(SCIP* scip, SCIP_EVENTHDLR *eventhdlr, SCIP_EVENT* event, SCIP_EVENTDATA *eventdata);
	virtual ~ExitEventHandler();
};

void registerExitEventHandler(ScipModelPtr scip);

};

#endif /* EXITEVENTHANDLER_H_ */
