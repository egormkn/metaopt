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


#ifndef DIRECTEDREACTION_H_
#define DIRECTEDREACTION_H_

#include <string>

#include "metaopt/model/Reaction.h"
#include "metaopt/Properties.h"

using namespace std;
using namespace boost;

namespace metaopt {

/**
 * Some reactions can be operated in both directions.
 * However, we can sometimes say only something about one particular direction.
 * For those cases we use DirectedReactions, where we additionally store also the direction of the reaction.
 * The fwd direction of a reaction is the direction of positive flux values.
 */
struct DirectedReaction {
	ReactionPtr _rxn;
	bool _fwd;

	DirectedReaction(ReactionPtr rxn, bool fwd) : _rxn(rxn), _fwd(fwd) {}

	bool operator==(DirectedReaction const& other) const {
		return _rxn == other._rxn && _fwd == other._fwd;
	}

	friend std::size_t hash_value(DirectedReaction const& p) {
		std::size_t seed = 0;
		boost::hash_combine(seed, p._rxn);
		boost::hash_combine(seed, p._fwd);
		return seed;
	}
};



}

#endif
