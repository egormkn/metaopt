/**
 * This header file stores global properties
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#define EPSILON 0.0001

#ifdef __CDT_PARSER__
	#define foreach(a, b) for(a : b)
	#define BOOST_THROW_EXCEPTION(x)
#else
    #define foreach BOOST_FOREACH
#endif

#endif
