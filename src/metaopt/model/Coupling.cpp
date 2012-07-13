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

Coupling::Coupling(const Coupling &c) : _fwd_couplings(c._fwd_couplings), _bwd_couplings(c._bwd_couplings) {
	// nothing else to do
}

Coupling::~Coupling() {
	// nothing to do
}

typedef unordered_multimap<DirectedReaction,DirectedReaction>::iterator citerator;

void Coupling::setCoupledFast(DirectedReaction a, DirectedReaction b) {
	// we already perform the transitive closure at insertion.
	// this makes setCoupled slow, but this is acceptable, because we will call isCoupled much more often than setCoupled.

	// create: a -> b
	unordered_set<DirectedReaction>& fwd = _fwd_couplings[a];
	unordered_set<DirectedReaction>& bwd = _bwd_couplings[b];
	fwd.insert(b);
	bwd.insert(a);
	/*foreach(DirectedReaction c, _bwd_couplings[a]) {
		// given: c -> a -> b
		// create: c -> b
		bwd.insert(c);
		_fwd_couplings[c].insert(b);
	}
	foreach(DirectedReaction d, _fwd_couplings[b]) {
		// givem: a -> b -> d
		// create: a -> d
		fwd.insert(d);
		_bwd_couplings[d].insert(a);
		foreach(DirectedReaction c, _bwd_couplings[a]) {
			// given: c -> a -> b -> d
			// create: c -> d
			_fwd_couplings[c].insert(d);
			_bwd_couplings[d].insert(c);
		}
	}*/
}
#if 0
void Coupling::computeClosure() {
	//TODO: this can be improved!
	cout << "starting compute closure... ";
	cout.flush();
	typedef pair<const DirectedReaction, unordered_set<DirectedReaction> > Entry;
	boost::unordered_set<DirectedReaction> keys;
	foreach(Entry & e, _fwd_couplings) {
		keys.insert(e.first);
	}
	foreach(Entry & e, _bwd_couplings) {
		keys.insert(e.first);
	}
	bool changeDetected = true;
	int iteration = 0;
	while(changeDetected) {
		iteration++;
		cout << "iteration "<< iteration << endl;
		changeDetected = false;
		foreach(DirectedReaction a, keys) {
			cout << a._rxn->getName() << endl;
			unordered_set<DirectedReaction>& fwdold = _fwd_couplings[a];
			unordered_set<DirectedReaction> fwdneu = fwdold; // copy set
			foreach(DirectedReaction d, fwdold) {
				unordered_set<DirectedReaction>& dfwd = _fwd_couplings[d];
				// given: a -> d -> _fwd_couplings[d]
				// create: a -> _fwd_couplings[d]
				fwdneu.insert(dfwd.begin(), dfwd.end());
			}
			if(fwdneu.size() != fwdold.size()) {
				fwdold = fwdneu;
				changeDetected = true;
			}
		}
		foreach(DirectedReaction a, keys) {
			cout << a._rxn->getName() << endl;
			unordered_set<DirectedReaction>& bwdold = _bwd_couplings[a];
			unordered_set<DirectedReaction> bwdneu = bwdold; // copy set
			foreach(DirectedReaction d, bwdold) {
				unordered_set<DirectedReaction>& dbwd = _bwd_couplings[d];
				bwdneu.insert(dbwd.begin(), dbwd.end());
			}
			if(bwdneu.size() != bwdold.size()) {
				bwdold = bwdneu;
				changeDetected = true;
			}
		}
	}
	cout << "finished" << endl;
}
#else
void Coupling::computeClosure() {
	cout << "starting compute closure";
	cout.flush();
	typedef pair<const DirectedReaction, unordered_set<DirectedReaction> > Entry;
	boost::unordered_set<DirectedReaction> keys;
	foreach(Entry & e, _fwd_couplings) {
		keys.insert(e.first);
	}
	int count = 0;
	foreach(DirectedReaction d, keys) {
		queue<DirectedReaction> q;
		unordered_set<DirectedReaction> visited;
		q.push(d);
		visited.insert(d);
		while(!q.empty()) {
			DirectedReaction a = q.front();
			q.pop();
			foreach(DirectedReaction b, _fwd_couplings[a]) {
				if(visited.find(b) == visited.end()) {
					q.push(b);
					visited.insert(b);
					// mark found couplings in _bwd_couplings to avoid concurrent modification issues
					_bwd_couplings[b].insert(d);
				}
			}
		}
		cout << ".";
		cout.flush();
		count++;
		if(count % 50 == 0) cout << endl;
	}
	// translate results from _bwd_couplings to _fwd_couplings
	foreach(Entry & e, _bwd_couplings) {
		foreach(DirectedReaction d, e.second) {
			_fwd_couplings[d].insert(e.first);
		}
	}
	cout << "finished" << endl;
}
#endif

bool Coupling::isCoupled(DirectedReaction a, DirectedReaction b)  {
	// since we already store the closure of the relation, a single test suffices.
	unordered_set<DirectedReaction>& fwd = _fwd_couplings[a];
	return(fwd.find(b) != fwd.end());
}

shared_ptr<vector<CoverReaction> > Coupling::computeCover(boost::unordered_set<DirectedReaction>& reactions) {
	boost::unordered_set<DirectedReaction> unanalyzed(reactions);

	// now we have to analyze the aggregated coupling information
	unordered_set<DirectedReaction> maxima;
	foreach(DirectedReaction rxn, reactions) {
		// canMax indicates if it can be the selected maximum
		bool canMax = true;

		boost::unordered_set<DirectedReaction>& bwd = _bwd_couplings[rxn];
		boost::unordered_set<DirectedReaction>& fwd = _fwd_couplings[rxn];
		foreach(DirectedReaction f, bwd) {
			if(canMax && reactions.find(f) != reactions.end()) {
				// f is interesting
				if(fwd.find(f) == reactions.end()) {
					// rxn does not cover f, hence rxn is not a maximal element
					canMax = false;
				}
				else {
					if(maxima.find(f) != maxima.end()) {
						canMax = false; // rxn may be a max, but we already found an equivalent one
					}
				}
			}
		}
		if(canMax) {
			maxima.insert(rxn);
		}
	}


	// we will misuse the list analyzed now to mark which reactions are still uncovered.
	assert(unanalyzed.size() == reactions.size()); // analyzed must contain every element of reactions

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
		boost::unordered_set<DirectedReaction>& fwd = _fwd_couplings[d];
		unanalyzed.erase(d);
		foreach(DirectedReaction f, fwd) {
			if(unanalyzed.find(f) != unanalyzed.end()) {
				c.covered->push_back(f);
				unanalyzed.erase(f);
			}
		}

		cover->push_back(c);
	}
	assert(unanalyzed.empty());
	return cover;
}

string Coupling::getStat() {
	stringstream ss;
	ss << string("number of Couplings: ") << _fwd_couplings.size();
	return ss.str();
}

CouplingPtr Coupling::copy() const{
	CouplingPtr res(new Coupling(*this));
	return res;
}

} /* namespace metaopt */
