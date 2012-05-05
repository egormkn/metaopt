#ifndef SCIPERROR_H_
#define SCIPERROR_H_

#include <boost/throw_exception.hpp>
#include <boost/exception/errinfo_at_line.hpp>
#include "scip/scip.h"
#include "metaopt/Properties.h"

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
		BOOST_THROW_EXCEPTION( ScipError() << error_code(boost_scip_call_error_code)); \
	} \
}

#endif

}

#endif