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

#ifndef SCIPERROR_H_
#define SCIPERROR_H_

#include <boost/throw_exception.hpp>
#include <boost/exception/errinfo_at_line.hpp>
#include "scip/scip.h"
#include "Properties.h"

#define ERROR_DEBUG

namespace metaopt {

struct ScipError : virtual boost::exception, virtual std::exception { };

typedef boost::error_info<struct tag_error_code,int> error_code;

#ifndef ERROR_DEBUG

#define BOOST_SCIP_CALL(x) { \
	int boost_scip_call_error_code = x; \
	if( boost_scip_call_error_code != SCIP_OKAY) { \
		BOOST_THROW_EXCEPTION( ScipError() << error_code(boost_scip_call_error_code)); \
	} \
}

#else

#define BOOST_SCIP_CALL(x) { \
	int boost_scip_call_error_code = x; \
	assert( boost_scip_call_error_code == SCIP_OKAY ); \
	if( boost_scip_call_error_code != SCIP_OKAY) { \
		BOOST_THROW_EXCEPTION( ScipError() << metaopt::error_code(boost_scip_call_error_code)); \
	} \
}

#endif

}

#endif
