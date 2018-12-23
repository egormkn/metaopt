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


#include "model/Reaction.h"
#include "model/DirectedReaction.h"
#include "Properties.h"

namespace metaopt {

struct CoverReaction {
	DirectedReaction reaction;
	boost::shared_ptr<std::vector<DirectedReaction> > covered;

	CoverReaction(const DirectedReaction& r) : reaction(r), covered(new std::vector<DirectedReaction>()) {}
};

typedef boost::shared_ptr<CoverReaction> CoverReactionPtr;

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

	// new implementation

	enum Color {
		WHITE,
		GREY,
		BLACK
	};

	struct Node {
		// adjacency list for a DirectedReaction
		const DirectedReaction _reaction;
		// use raw pointers, because the alternative would be weak pointers and thats too much hassle!
		boost::unordered_set<Node*> _to; // nodes n with this coupled to n
		boost::unordered_set<Node*> _from; // nodes n with n coupled to this
		Color _color;
		unsigned int _finishtime;

		Node(DirectedReaction& d);
	};

	std::size_t hash_value(Node* p );

	typedef boost::shared_ptr<Node> NodePtr;

	boost::unordered_map<DirectedReaction, NodePtr> _nodes;

	struct StrongComponent {
		// stores list of strong components this is coupled to
		boost::unordered_set<StrongComponent*> _coupledTo;
	};

	std::size_t hash_value(StrongComponent* p );

	typedef boost::shared_ptr<StrongComponent> StrongComponentPtr;

	boost::unordered_map<DirectedReaction, StrongComponentPtr> _components;

	int num_components; // for statistical purposes only

	// indicates if closure was already computed
	bool israw;

	/**
	 * Runs depth first search starting from the given node and stores finishing times in nodes.
	 * Only searches nodes reachable from node.
	 *
	 * This is the dfs run in the first round.
	 *
	 * @param node starting node
	 * @param next finsih time to be assigned (gets updated by algorithm)
	 */
	void dfs(Node* node, unsigned int &time);

	/**
	 * This is the dfs run in the second round.
	 *
	 * It operates on the reverse arcs and stores the found nodes in the given strong component.
	 */
	void dfs2(Node* node, const StrongComponentPtr & comp);
};

typedef boost::shared_ptr<Coupling> CouplingPtr;

} /* namespace metaopt */
#endif /* COUPLING_H_ */
