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
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>


namespace metaopt {

using namespace std;
using namespace boost;

Coupling::Coupling() {
	// nothing to do
}

Coupling::Coupling(const Coupling &c) : _G(c._G), _vertices(c._vertices) {
	// nothing else to do
	israw = true;
}

Coupling::~Coupling() {
	// nothing to do
}

typedef unordered_multimap<DirectedReaction,DirectedReaction>::iterator citerator;

void Coupling::addCoupled(DirectedReaction a, DirectedReaction b) {
	vertex_t va, vb;
	if(_vertices.find(a) == _vertices.end()) {
		va = add_vertex(_G);
		_vertices[a] = va;
	}
	else {
		va = _vertices[a];
	}
	if(_vertices.find(b) == _vertices.end()) {
		vb = add_vertex(_G);
		_vertices[b] = vb;
	}
	else {
		vb = _vertices[b];
	}
	add_edge(va, vb, _G);
	israw = true;
}

void Coupling::computeClosure() {
#ifndef SILENT
	cout << "starting compute closure... ";
	cout.flush();
#endif
	_TC.clear();
	_to_tc_vec = std::vector<tc_vertex_t>(num_vertices(_G));;
	/*_g_to_tc.clear();

	struct map {
		unordered_map<vertex_t, tc_vertex_t>& _tomap;

		map(unordered_map<vertex_t, tc_vertex_t>& tomap) : _tomap(tomap) {}

		tc_vertex_t& operator[](vertex_t v) {
			return _tomap[v];
		}
	};

	map m(_g_to_tc);*/

	typedef boost::property_map<graph_t, boost::vertex_index_t>::const_type VertexIndexMap;
	VertexIndexMap index_map = get(vertex_index, _G);

	iterator_property_map < tc_vertex_t *, VertexIndexMap, tc_vertex_t, tc_vertex_t&> g_to_tc(&_to_tc_vec[0], index_map);


	transitive_closure(_G, _TC, g_to_tc, get(vertex_index, _G));
	//transitive_closure(_G, _TC);
#ifndef SILENT
	cout << "finished" << endl;
#endif
	israw = false;
}


bool Coupling::isCoupled(DirectedReaction a, DirectedReaction b) {
	assert(!israw);
	if(_vertices.find(a) == _vertices.end()) return false;
	if(_vertices.find(b) == _vertices.end()) return false;

	typedef boost::property_map<graph_t, boost::vertex_index_t>::const_type VertexIndexMap;
	VertexIndexMap index_map = get(vertex_index, _G);

	iterator_property_map < tc_vertex_t *, VertexIndexMap, tc_vertex_t, tc_vertex_t&> g_to_tc(&_to_tc_vec[0], index_map);

	try {
		tc_vertex_t va = g_to_tc[_vertices.at(a)];
		tc_vertex_t vb = g_to_tc[_vertices.at(b)];
		bool res = is_adjacent(_TC, va, vb);
		return res;
	}
	catch(std::exception& ex) {
		cout << diagnostic_information(ex) << endl;
		/*typedef pair<DirectedReaction, vertex_t> elem;
		foreach(elem d, _vertices) {
			assert(g_to_tc.find(d.second) != g_to_tc.end());
		}*/
		assert(false);
	}
	return false;
}

shared_ptr<vector<CoverReaction> > Coupling::computeCover(boost::unordered_set<DirectedReaction>& reactions) {
	assert(!israw);

#define NOTCOMPUTE 0
#if NOTCOMPUTE
	shared_ptr<vector<CoverReaction> > cover(new vector<CoverReaction>());
	foreach(DirectedReaction d, reactions) {
		CoverReaction c(d);
		cover->push_back(c);
	}
#else

	boost::unordered_set<DirectedReaction> unanalyzed(reactions);

	// now we have to analyze the aggregated coupling information
	unordered_set<DirectedReaction> maxima;
	foreach(DirectedReaction rxn, reactions) {
		// canMax indicates if it can be the selected maximum
		bool canMax = true;

#if 1
		foreach(DirectedReaction d, reactions) {
			if(isCoupled(d,rxn)) {
				if(isCoupled(rxn,d)) {
					if(maxima.find(d) != maxima.end()) {
						canMax = false; // rxn may be a max, but we already found an equivalent one
					}
				}
				else {
					canMax = false;
				}
			}
		}

#else
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
#endif
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
		unanalyzed.erase(d);
		for(unordered_set<DirectedReaction>::iterator i = unanalyzed.begin(); i != unanalyzed.end();) {
			if(isCoupled(d,*i)) {
				c.covered->push_back(*i);
				i = unanalyzed.erase(i);
			}
			else {
				i++;
			}
		}
#if 0
		boost::unordered_set<DirectedReaction>& fwd = _fwd_couplings[d];
		unanalyzed.erase(d);
		foreach(DirectedReaction f, fwd) {
			if(unanalyzed.find(f) != unanalyzed.end()) {
				c.covered->push_back(f);
				unanalyzed.erase(f);
			}
		}
#endif
		cover->push_back(c);
	}
	assert(unanalyzed.empty());

#endif
	return cover;
}

string Coupling::getStat() {
	stringstream ss;
	ss << string("number of raw couplings: ") << num_edges(_G) << " number of aggregated couplings: " << num_edges(_TC);
	return ss.str();
}

CouplingPtr Coupling::copy() const{
	CouplingPtr res(new Coupling(*this));
	return res;
}

} /* namespace metaopt */
