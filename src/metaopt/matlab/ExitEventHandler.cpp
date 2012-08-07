/*
 * ExitEventHandler.cpp
 *
 *  Created on: Jan 17, 2012
 *      Author: arne
 */

#include "ExitEventHandler.h"
#include "mex.h"

#include <iostream>

using namespace std;
using namespace scip;

#define EXIT_EVENT_HANDLER_NAME "ExitHandler"
#define EXIT_EVENT_HANDLER_DESC "on user request, aborts computation"

extern "C" bool utIsInterruptPending();
//extern "C" void utSetInterruptPending();

namespace metaopt {

ExitEventHandler::ExitEventHandler( SCIP*  scip )
	: ObjEventhdlr(scip, EXIT_EVENT_HANDLER_NAME, EXIT_EVENT_HANDLER_DESC) {
}

ExitEventHandler::~ExitEventHandler() {
	// nothing to do
}

SCIP_RETCODE ExitEventHandler::scip_exec(SCIP* scip, SCIP_EVENTHDLR *eventhdlr, SCIP_EVENT* event, SCIP_EVENTDATA *eventdata) {
	//cout << "executing event handler" << endl;
	if(utIsInterruptPending()) {
		mexPrintf("\n Ctrl-C detected. Aborting SCIP computation\n");
		SCIP_CALL( SCIPinterruptSolve(scip) );
	}
	return SCIP_OKAY;
}

void registerExitEventHandler(ScipModelPtr scip) {
	SCIP_EVENTHDLR* evthdlr = SCIPfindEventhdlr(scip->getScip(), EXIT_EVENT_HANDLER_NAME);
	//cout << "found eventhandler? "<<evthdlr;
	if(evthdlr == NULL) {
	      SCIPerrorMessage(EXIT_EVENT_HANDLER_NAME);
	      SCIPerrorMessage(" not found\n");
	      BOOST_SCIP_CALL( SCIP_PLUGINNOTFOUND );
	}
	BOOST_SCIP_CALL( SCIPcatchEvent(scip->getScip(), SCIP_EVENTTYPE_NODESOLVED, evthdlr, NULL, NULL) );
}

};
