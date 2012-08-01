
#ifndef DIRECTEDREACTION_H_
#define DIRECTEDREACTION_H_

#include <string>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "metaopt/model/Reaction.h"

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
