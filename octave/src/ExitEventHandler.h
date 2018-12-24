/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    fast-tfva - efficient thermodynamic constrained flux variability analysis.
    Copyright (C) 2012  Arne Müller, arne.mueller@fu-berlin.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * ExitEventHandler.h
 *
 *  Created on: Jan 17, 2012
 *      Author: arne
 */

#ifndef EXITEVENTHANDLER_H_
#define EXITEVENTHANDLER_H_

#include "objscip/objscip.h"

#include "model/scip/ScipModel.h"
#include "scip/ScipError.h"

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
