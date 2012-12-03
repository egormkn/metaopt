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
 * Coupling.cpp
 *
 *  Created on: 10.07.2012
 *      Author: arnem
 */

#include "metaopt/Properties.h"

#include "Coupling.h"
#include <iostream>
#include <utility>                   // for std::pair
#include <list>


namespace metaopt {

using namespace std;
using namespace boost;

Coupling::Coupling() {
	israw = true;
}

Coupling::Coupling(const Coupling &c) {
	israw = true;
	// only copy nodes, and force the other stuff to be recomputed
	typedef pair<DirectedReaction, NodePtr> NodeEntry;

	// create the nodes
	foreach(NodeEntry e, c._nodes) {
		NodePtr& n = _nodes[e.first];
		n.reset(new Node(e.first));
	}
	// add the links
	foreach(NodeEntry e, c._nodes) {
		const NodePtr& n = _nodes[e.first];
		foreach(Node* k, n->_to) {
			n->_to.insert(_nodes[k->_reaction].get());
		}
		foreach(Node* k, n->_from) {
			n->_from.insert(_nodes[k->_reaction].get());
		}
		// color and finish time we don't have to set, since they are set in the algorithm when needed
	}

	// components get computed, so don't copy
}

Coupling::~Coupling() {
	// nothing to do
}

//typedef unordered_multimap<DirectedReaction,DirectedReaction>::iterator citerator;

std::size_t Coupling::hash_value(Node* p ) {
	return reinterpret_cast<std::size_t>(p);
}

std::size_t Coupling::hash_value(StrongComponent* p ) {
	return reinterpret_cast<std::size_t>(p);
}

Coupling::Node::Node(DirectedReaction& d) : _reaction(d) {}

void Coupling::addCoupled(DirectedReaction a, DirectedReaction b) {

	// use the property that unordered_map allocates new objects if it doesn't find them
	NodePtr& na = _nodes[a];
	NodePtr& nb = _nodes[b];

	// create nodes if necessary
	if(na.use_count() == 0) {
		na.reset(new Node(a));
	}
	if(nb.use_count() == 0) {
		nb.reset(new Node(b));
	}

	na->_to.insert(nb.get());
	nb->_from.insert(na.get());

	israw = true;
}

void Coupling::dfs(Node* node, unsigned int &time) {
	node->_color = GREY;
	foreach(Node* child, node->_to) {
		if(child->_color == WHITE) {
			dfs(child, time);
		}
	}
	node->_color = BLACK;
	node->_finishtime = time++;
}

void Coupling::dfs2(Node* node, const StrongComponentPtr & comp) {
	node->_color = GREY;
	_components[node->_reaction] = comp; // store this connected component to belong to the node
	foreach(Node* child, node->_from) {
		if(child->_color == WHITE) {
			dfs2(child, comp);
		}
	}
	node->_color = BLACK;
}

void Coupling::computeClosure() {
#ifndef SILENT
	cout << "starting compute closure... ";
	cout.flush();
#endif

	/*
	 * We compute the strong components (correspond to sets of fully coupled reactions).
	 * As a side result, we also get a topological ordering which we use to compute the transitive couplings.
	 */

	// reset components
	_components.clear();

	typedef pair<DirectedReaction, NodePtr> NodeEntry;
	// run first dfs
	// first reset nodes (turn them all white)
	foreach(NodeEntry e, _nodes) {
		e.second->_color = WHITE;
	}
	unsigned int time = 0; // we start with zero, so that we can use the finish times as indicies
	foreach(NodeEntry e, _nodes) {
		if(e.second->_color == WHITE) {
			dfs(e.second.get(), time);
		}
	}
	assert(time == _nodes.size());

	// run second dfs
	// the second dfs has to process the nodes in order of decreasing finish time
	// so we create an ordered list of the nodes
	vector<Node*> ordered;
	ordered.resize(_nodes.size(), NULL);
	// reset nodes (turn them all white)
	// and put them in their correct position
	foreach(NodeEntry e, _nodes) {
		e.second->_color = WHITE;
		assert(e.second->_finishtime >= 0);
		assert(e.second->_finishtime < _nodes.size());
		ordered[e.second->_finishtime] = e.second.get();
	}
	vector<StrongComponentPtr> comps; // stores topologically ordered strong components
	for(int i = ordered.size()-1; i >= 0; i--) {
		Node* n = ordered[i];
		assert(n != NULL);
		if(n->_color == WHITE) {
			StrongComponentPtr comp(new StrongComponent());
			dfs2(n, comp);
			comps.push_back(comp);
		}
	}
	num_components = comps.size(); // store number of components for statistical purposes

	// iterate through all nodes and transfer directed couplings to strong components
	foreach(NodeEntry e, _nodes) {
		const StrongComponentPtr& s = _components[e.first];
		foreach(Node* n, e.second->_to) {
			const StrongComponentPtr& t = _components[n->_reaction];
			if(t != s) s->_coupledTo.insert(t.get()); // only add couplings to different components
		}
	}
	// now we have a topologically ordered list of components
	// we can now use a dynamic programming variant to add the partial couplings.
	// to do so, we iterate through the list in reverse order
	for(int i = comps.size()-1; i >= 0; i--) {
		const StrongComponentPtr& s = comps[i];
		unordered_set<StrongComponent*> to = s->_coupledTo; // make a copy, so we can modify the original while iterating
		// since we iterate through the components in topological order,
		// we already know the transitive closure for couplings from each t in to
		// hence, this simple update is sufficient
		foreach(StrongComponent* t, to) {
			s->_coupledTo.insert(t->_coupledTo.begin(), t->_coupledTo.end());
		}
	}

#ifndef SILENT
	cout << "finished" << endl;
#endif
	israw = false;
}


bool Coupling::isCoupled(DirectedReaction a, DirectedReaction b) {
	assert(!israw);

	const StrongComponentPtr& sa = _components[a];
	const StrongComponentPtr& sb = _components[b];

	return(sa == sb || sa->_coupledTo.find(sb.get()) != sa->_coupledTo.end());
}

shared_ptr<vector<CoverReaction> > Coupling::computeCover(boost::unordered_set<DirectedReaction>& reactions) {
	assert(!israw);

#define NOTCOMPUTE 0
#if NOTCOMPUTE
	shared_ptr<vector<CoverReaction> > cover(new vector<CoverReaction>());
	foreach(const DirectedReaction& d, reactions) {
		CoverReaction c(d);
		cover->push_back(c);
	}
#else

	shared_ptr<vector<CoverReaction> > cover(new vector<CoverReaction>());

	list<DirectedReaction> rxns(reactions.begin(), reactions.end());

	// we have a directed acyclic graph of the components;
	// hence, we just have to find components that point to no components of the set
	// use a hashmap to store the components, since it directly eliminates duplicates
	unordered_set<StrongComponent*> comps; // use pointers, because for those we have a hash function
	for(list<DirectedReaction>::iterator i = rxns.begin(); i != rxns.end();) {
		DirectedReaction d = *i;
		unordered_map<DirectedReaction, StrongComponentPtr>::iterator iter = _components.find(d);
		if(iter != _components.end()) {
			StrongComponentPtr s = iter->second;
			assert(s.use_count() > 0);
			comps.insert(s.get());
			i++;
		}
		else {
			// this really can happen
			// for example if the reaction is uncoupled to every other reaction
			cover->push_back(CoverReaction(d));
			i = rxns.erase(i); // returns iterator to element after erased element
		}
	}

	foreach(StrongComponent* c, comps) {
		bool isMax = true;
		foreach(StrongComponent* t, c->_coupledTo) {
			if(comps.find(t) != comps.end()) {
				// we found a component to which c is coupled, hence it is not maximal
				isMax = false;
				break;
			}
		}
		if(isMax) {
			shared_ptr<DirectedReaction> m;
			shared_ptr<vector<DirectedReaction> > covered(new vector<DirectedReaction>());
			for(list<DirectedReaction>::iterator i = rxns.begin(); i != rxns.end();) {
				DirectedReaction d = *i;
				const StrongComponentPtr& cd = _components[d];
				if(cd.get() == c) {
					if(m.use_count() == 0) {
						m.reset(new DirectedReaction(d));
						i = rxns.erase(i);
					}
					else {
						covered->push_back(d); // reaction is covered
						i = rxns.erase(i);
					}
				}
				else {
					// we have to check if the reaction is covered
					if(cd->_coupledTo.find(c) != cd->_coupledTo.end()) {
						// it is coupled
						covered->push_back(d);
						i = rxns.erase(i);
					}
					else {
						i++;
					}
				}
			}
			assert(m.use_count() != 0);
			CoverReaction cr(*m);
			cr.covered = covered;
			cover->push_back(cr);
		}
	}
#ifndef NDEBUG
	if(!rxns.empty()) {
		foreach(const DirectedReaction& d, rxns) {
			cout << "retained: " << d._rxn->getName() << endl;
		}

		foreach(const CoverReaction& cr, *cover) {
			cout << cr.reaction._rxn->getName() << " ( ";
			foreach(const DirectedReaction& d, *(cr.covered)) {
				cout << d._rxn->getName() << " ";
			}
			cout << ")" << endl;
		}
	}
#endif
	assert(rxns.empty());
#endif

	return cover;
}

string Coupling::getStat() {
	stringstream ss;
	int num_raw = 0;

	typedef pair<DirectedReaction, NodePtr> NodeEntry;
	// run first dfs
	// first reset nodes (turn them all white)
	foreach(NodeEntry e, _nodes) {
		num_raw += e.second->_to.size();
	}

	if(israw) {
		ss << string("number of raw couplings: ");
		return ss.str();
	}
	else {
		// don't compute number of aggregated couplings, because that is too much overhead
		ss << string("number of raw couplings: ") << num_raw << " number of strong components: " << num_components;
	}
	return ss.str();
}

CouplingPtr Coupling::copy() const{
	CouplingPtr res(new Coupling(*this));
	return res;
}

} /* namespace metaopt */
