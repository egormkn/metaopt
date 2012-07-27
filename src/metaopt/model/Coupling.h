/*
 * Coupling.h
 *
 *  Created on: 10.07.2012
 *      Author: arnem
 */

#ifndef COUPLING_H_
#define COUPLING_H_

#include <string>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "metaopt/model/Reaction.h"
#include "metaopt/model/DirectedReaction.h"

#ifdef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
#error The transitive closure algorithm uses partial specialization.
#endif

namespace metaopt {

struct CoverReaction {
	DirectedReaction reaction;
	boost::shared_ptr<std::vector<DirectedReaction> > covered;

	CoverReaction(DirectedReaction& r) : reaction(r), covered(new std::vector<DirectedReaction>()) {}
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
	 * copies the coupling table
	 */
	Coupling(const Coupling& c);

	/**
	 * Stores that a is directionally coupled to b.
	 * This means that if a carries (positive) flux, also b carries (positive) flux.
	 * This is the fast implementation that does not maintain transitive closure.
	 */
	void addCoupled(DirectedReaction a, DirectedReaction b);

	/**
	 * Computes the transitive closure.
	 */
	void computeClosure();

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

	/**
	 * return some status information on the couplings
	 */
	std::string getStat();

	/**
	 * copy this Coupling information
	 */
	boost::shared_ptr<Coupling> copy() const;

private:

	// type for input graph
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS > graph_t;
	typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_t;
	graph_t _G;

	// type for transitive closure graph
	typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS> tc_graph_t;
	typedef boost::graph_traits < tc_graph_t >::vertex_descriptor tc_vertex_t;
	tc_graph_t _TC;

	std::vector<tc_vertex_t> _to_tc_vec;


	boost::unordered_map<DirectedReaction, vertex_t> _vertices;

	//boost::unordered_map<vertex_t, tc_vertex_t> _g_to_tc;

	bool israw;
};

typedef boost::shared_ptr<Coupling> CouplingPtr;

} /* namespace metaopt */
#endif /* COUPLING_H_ */
