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

#define VERSION "2.0.4"


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


#if 0
/**
 * Implement hashing support for MetabolitePtr and ReactionPtr
 */
namespace metaopt
{
#if 0 //BOOST_VERSION < 104700
   // somewhere between version 1.46 and 1.47 boost starts to have a definition for hash_value(shared_ptr<T>)
   // if it does not, we have to define it
   template <class T>
   std::size_t
   hash_value(boost::shared_ptr<T> const & _ptr)
   {
      return reinterpret_cast<std::size_t>( _ptr.get() );
   }

#else
   // otherwise the hash_function simply calls the native pointer hash function
   // so we have to define that
   std::size_t hash_value(MetabolitePtr const & met ) {
	   return reinterpret_cast<std::size_t>(met.get());
   }

   std::size_t hash_value(ReactionPtr const & rxn ) {
	   return reinterpret_cast<std::size_t>(rxn.get());
   }
#endif
}
#endif

//#define _GLIBCXX_DEBUG

#endif
