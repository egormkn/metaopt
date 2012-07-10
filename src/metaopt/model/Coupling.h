/*
 * Coupling.h
 *
 *  Created on: 10.07.2012
 *      Author: arnem
 */

#ifndef COUPLING_H_
#define COUPLING_H_

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include "metaopt/model/Reaction.h"


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

struct CoverReaction {
	DirectedReaction reaction;
	boost::shared_ptr<std::vector<DirectedReaction> > covered;
};

/**
 * This class manages directed flux coupling information on reactions.
 *
 * Flux coupling is a partial ordering a -> b with the semantic that if a carries (positive) flux, also b carries positive flux.
 * We will also speak about maximal elements, where a is bigger than b if a -> b.
 */
class Coupling {
public:
	Coupling();
	virtual ~Coupling();

	/**
	 * Stores that a is directionally coupled to b.
	 * This means that if a carries (positive) flux, also b carries (positive) flux.
	 */
	void setCoupled(DirectedReaction a, DirectedReaction b);

	/**
	 * Checks, if a is directionally coupled to b.
	 * This also uses the transitivity property of flux coupling.
	 * Hence, the entry does not need to have been added directly to the list of couplings.
	 */
	bool isCoupled(DirectedReaction a, DirectedReaction b);

	/**
	 * Given a set of reactions A, computes a subset B of A such that if every reaction in B carries (positive) flow,
	 * every reaction in A carries (positive) flow.
	 *
	 * The set B is a subset of the maximal elements of A (w.r.t. a >= b <=> a -> b)
	 *
	 * The result is the list of reactions B. For each b in B a list of covered reactions is stored, such that each reaction in A is either in B or covered.
	 */
	boost::shared_ptr<std::vector<CoverReaction> > computeCover(boost::unordered_set<DirectedReaction>& reactions);

private:
	/**
	 * This map stores for every directed reaction b the list of directed reactions a with a -> b.
	 */
	boost::unordered_multimap<DirectedReaction, DirectedReaction> couplings;

};

typedef boost::shared_ptr<Coupling> CouplingPtr;

} /* namespace metaopt */
#endif /* COUPLING_H_ */
