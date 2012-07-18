/**
 * This header file stores global properties
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#define EPSILON 0.0001

#define LPSOLVER_SOPLEX

#ifdef __CDT_PARSER__
	#define foreach(a, b) for(a : b)
	#define BOOST_THROW_EXCEPTION(x)
#else
    #define foreach BOOST_FOREACH
#endif

//#define _GLIBCXX_DEBUG

#endif
