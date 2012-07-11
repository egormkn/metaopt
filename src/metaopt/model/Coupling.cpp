/*
 * Coupling.cpp
 *
 *  Created on: 10.07.2012
 *      Author: arnem
 */

#include "Coupling.h"
#include <queue>
#include <iostream>

namespace metaopt {

using namespace std;
using namespace boost;

Coupling::Coupling() {
	// nothing to do
}

Coupling::Coupling(const Coupling &c) : couplings(c.couplings) {
	// nothing else to do
}

Coupling::~Coupling() {
	// nothing to do
}

typedef unordered_multimap<DirectedReaction,DirectedReaction>::iterator citerator;

void Coupling::setCoupled(DirectedReaction a, DirectedReaction b) {
	// TODO: theoretically we could implement additional means to cleanup the data structure to improve the running time of the retrieval methods.

	// first check, if the element has already been inserted
	pair<citerator, citerator> res = couplings.equal_range(b);
	for(citerator i = res.first; i != res.second; i++) {
		if(i->second == a) {
			return; // we don't have to add it
		}
	}
	// the element was not found, so add it.
	couplings.insert(pair<DirectedReaction, DirectedReaction>(b,a));
}

bool Coupling::isCoupled(DirectedReaction a, DirectedReaction b)  {
	// starting from b we do a graph search.
	// since the poset is not necessarily a tree, we will remember the elements we visited
	unordered_set<DirectedReaction> visited;
	visited.insert(b);
	queue<DirectedReaction> q;
	q.push(b);
	while(!q.empty()) {
		DirectedReaction d = q.front();
		q.pop();
		pair<citerator, citerator> res = couplings.equal_range(d);
		for(citerator i = res.first; i != res.second; i++) {
			DirectedReaction& toAdd = i->second;
			if(visited.find(toAdd) == visited.end()){
				if(toAdd == a) {
					return true; // we have found a transitive chain to a
				}
				// we have not yet added i->second to the queue
				q.push(toAdd);
				visited.insert(toAdd);
			}
		}
	}
	// we have not found a, so a is not directionally coupled to b
	return false;
}

inline void insertUnique(unordered_multimap<DirectedReaction, DirectedReaction>& map, DirectedReaction& key, DirectedReaction& val) {
	pair<citerator, citerator> res = map.equal_range(key);
	bool found = false;
	for(citerator i = res.first; !found && i != res.second; i++) {
		found = i->second == val;
	}
	if(!found) {
		map.insert(pair<DirectedReaction, DirectedReaction>(key, val));
	}
}


shared_ptr<vector<CoverReaction> > Coupling::computeCover(boost::unordered_set<DirectedReaction>& reactions) {
	// for each directed reaction a in the target set we do a graph search.
	// If we find another directed reaction b of the set during this search, a is only maximal if a <-> b.

	// first we aggregate the coupling information on only the reactions that we are interested in.
	// in particular, we resolve all transitivities.
	boost::unordered_multimap<DirectedReaction, DirectedReaction> targetCoupled;
	boost::unordered_set<DirectedReaction> analyzed;
#if 0
	boost::unordered_set<DirectedReaction> notMaximal;


	foreach(DirectedReaction rxn, reactions) {
		analyzed.insert(rxn);

		bool canMaximal = true; // indicates that we have not found evidence yet proofing that rxn is not maximal

		// since the poset is not necessarily a tree, we will remember the elements we visited
		unordered_set<DirectedReaction> visited;
		visited.insert(rxn);
		queue<DirectedReaction> q;
		q.push(rxn);
		while(canMaximal && !q.empty()) {
			DirectedReaction d = q.front();
			q.pop();
			pair<citerator, citerator> res = couplings.equal_range(d);
			for(citerator i = res.first; canMaximal && i != res.second; i++) {
				DirectedReaction& toAdd = i->second;
				if(visited.find(toAdd) == visited.end()){
					if(reactions.find(toAdd) != reactions.end()) {
						// we have found a transitive chain to an element of the set
						if(notMaximal.find(toAdd) != notMaximal.end()) {
							notMaximal.insert(rxn);
							canMaximal = false; // rxn cannot be maximal anymore
							pair<citerator, citerator> toAddCoupled = targetCoupled.equal_range(toAdd);
							for(citerator j = toAddCoupled.first; j != toAddCoupled.second; j++) {
								insertUnique(targetCoupled, rxn, j->second); // insert to know from which maximal reaction rxn is covered.
							}
						}
						else if(analyzed.find(toAdd) != analyzed.end()){
							pair<citerator, citerator> back = targetCoupled.equal_range(toAdd);
							bool found = false; // found back coupling
							for(citerator j = back.first; j != back.second; j++) {
								if(j->second == rxn) found = true;
								else {
									// because of transitivity we can reuse the result and don't have to investigate the subtree.
									if(visited.find(j->second) == visited.end()) {
										insertUnique(targetCoupled, rxn, j->second);
										visited.insert(j->second);
									}
								}
							}
							if(!found) {
								notMaximal.insert(rxn);
								canMaximal = false;
							}
						}
						else {
							// we will have to analyze the tree further along this node
							q.push(toAdd);
							visited.insert(toAdd);
						}
						if(canMaximal) {
							insertUnique(targetCoupled, rxn, toAdd);
							// and continue searching
						}
					} // end test if contained in reactions
					else {
						// we will have to analyze the tree further along this node
						q.push(toAdd);
						visited.insert(toAdd);
					}
				} // end test if visited
			}
		}
	}
#endif

	foreach(DirectedReaction a, reactions) {
		analyzed.insert(a);
		foreach(DirectedReaction b, reactions) {
			if(isCoupled(a,b)) {
				targetCoupled.insert(pair<DirectedReaction, DirectedReaction>(b,a));
			}
		}
	}

	// now we have to analyze the aggregated coupling information
	unordered_set<DirectedReaction> maxima;
	foreach(DirectedReaction rxn, reactions) {
		// canMax indicates if it can be the selected maximum
#if 0
		bool canMax = notMaximal.find(rxn) == notMaximal.end();
#else
		bool canMax = true;
#endif
		if(!canMax) continue;

		pair<citerator, citerator> res = targetCoupled.equal_range(rxn);
		for(citerator i = res.first; canMax && i != res.second; i++) {
			DirectedReaction d = i->second;
#if 0
			canMax = notMaximal.find(d) == notMaximal.end();
			if(!canMax) {
				notMaximal.insert(rxn);
				continue;
			}
#endif
			// Question the alg will answer: is d -> rxn strict?
			// Answer: this is the case if we do not find rxn -> d in the aggregated coupling information
			bool found = false;
			pair<citerator, citerator> back = targetCoupled.equal_range(d);
			for(citerator j = back.first; !found && j != back.second; j++) {
				if(j->second == rxn) found = true;
			}
			if(!found) {
				canMax = false;
#if 0
				notMaximal.insert(rxn);
#endif
				continue;
			}
			else {
				if(maxima.find(d) != maxima.end()) {
					canMax = false; // it may be a max, but we already found an equivalent one
				}
			}
		}
		if(canMax) {
			maxima.insert(rxn);
		}
	}


	// we will misuse the list analyzed now to mark which reactions are still uncovered.
	assert(analyzed.size() == reactions.size()); // analyzed must contain every element of reactions

#if 0
	cout << "maxima: ";
	foreach(DirectedReaction d, maxima) {
		cout << d._rxn->getName() << " ";
	}
	cout << endl;
#endif
#if 0
	foreach(DirectedReaction d, reactions) {
		cout << d._rxn->getName() << ": ";
		pair<citerator, citerator> res = targetCoupled.equal_range(d);
		for(citerator j = res.first; j != res.second; j++) {
			cout << j->second._rxn->getName() << " ";
		}
		cout << " ( ";
		res = couplings.equal_range(d);
		for(citerator j = res.first; j != res.second; j++) {
			cout << j->second._rxn->getName() << " ";
		}
		cout << " ) " ;
		cout << endl;
	}
#endif

	shared_ptr<vector<CoverReaction> > cover(new vector<CoverReaction>());
	foreach(DirectedReaction d, maxima) {
		CoverReaction c(d);
		analyzed.erase(d);
		for(unordered_set<DirectedReaction>::iterator i = analyzed.begin(); i != analyzed.end();) {
			// check if d >= *i
			bool found = false;
			pair<citerator, citerator> res = targetCoupled.equal_range(*i);
			for(citerator j = res.first; !found && j != res.second; j++) {
				found = d == j->second;
			}
			if(found) {
				c.covered->push_back(*i);
				i = analyzed.erase(i);
			}
			else {
				i++;
			}
		}
		cover->push_back(c);
	}
	assert(analyzed.empty());
	return cover;
}

string Coupling::getStat() {
	stringstream ss;
	ss << string("number of Couplings: ") << couplings.size();
	return ss.str();
}

CouplingPtr Coupling::copy() const{
	CouplingPtr res(new Coupling(*this));
	return res;
}

} /* namespace metaopt */
