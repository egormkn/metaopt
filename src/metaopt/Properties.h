/**
 * This header file stores global properties
 */

#include <cstddef>
#include <boost/shared_ptr.hpp>

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#define EPSILON 0.0001

#define LPSOLVER_SOPLEX

#ifdef __CDT_PARSER__
	#define foreach(a, b) for(a : b)
	#define BOOST_THROW_EXCEPTION(x)
	#define CLOCKS_PER_SEC 0
#else
    #define foreach BOOST_FOREACH
#endif

namespace boost
{
   template <class T>
   std::size_t
   hash_value(boost::shared_ptr<T> const & _ptr)
   {
      return reinterpret_cast<std::size_t>( _ptr.get() );
   }
}

//#define _GLIBCXX_DEBUG

#endif
