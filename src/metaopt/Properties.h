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

/**
 * This header file stores global properties
 */

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

//#define EPSILON 0.0001

#define LPSOLVER_SOPLEX

#define VERSION "2.0.2"


#ifdef __CDT_PARSER__
	#define foreach(a, b) for(a : b)
	#define BOOST_THROW_EXCEPTION(x)
	#define CLOCKS_PER_SEC 0
#else
    #define foreach BOOST_FOREACH
#endif

/** Thrown id an endless loop has been detected */
struct EndlessLoopError : virtual boost::exception, virtual std::exception {
public:
	EndlessLoopError() {}
};

/** If an endless loop is detected, this error message indicates the number of iterations after the loop was aborted */
typedef boost::error_info<struct tag_iteration_count,int> iteration_count;


namespace metaopt
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
