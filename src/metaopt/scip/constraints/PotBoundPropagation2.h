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
 * PotBoundPropagation2.h
 *
 *  Created on: 18.05.2012
 *      Author: arnem
 */

#ifndef POTBOUNDPROPAGATION2_H_
#define POTBOUNDPROPAGATION2_H_

#include <boost/unordered_map.hpp>
#include <vector>
#include <queue>
#include <cmath>
#include "metaopt/model/Model.h"
#include "metaopt/model/scip/ScipModel.h"
#include "objscip/objscip.h"
#include "metaopt/Properties.h"

namespace metaopt {

class PotBoundPropagation2 {
public:
	PotBoundPropagation2(ModelPtr model);
	virtual ~PotBoundPropagation2();

	/**
	 * prints updated values
	 */
	void print();

	/**
	 * checks if the computed bounds are really a correct bound.
	 */
	void check();

	/**
	 * updates the computed bound to be a good mu-bound.
	 * Uses the fixed directions of current status of scip.
	 * If we cannot be sure if the current bound is a proper mu-bound, computation is started from scratch
	 */
	void update(ScipModelPtr scip);

	/**
	 * updates the computed bound to be a good mu-bound.
	 * Does not use the fixed directions of current status of scip.
	 * This only runs the flow update and uses the current solution as initial solution.
	 */
	void update();

	// from the computed potential bounds lists all blocked reactions
	// the boolean in each pair is true, if the forward direction is blocked.
	// it is false, if the backward direction is blocked.
	// if a reaction is blocked in both directions, we have two entries.
	boost::shared_ptr<std::vector<std::pair<ReactionPtr,bool> > > getBlockedReactions();

private:
	struct Arc;

	struct MetBound {
		MetabolitePtr _met;
		bool _isMinBound;
		// true, if this MetBound is a boundary metabolite (this has nothing to do with potential bounds!)
		bool _isBoundary;

		// value of the bound in transformed setting
		double _bound;
		boost::weak_ptr<Arc> _last_update; // only logging information, so the results can be verified and/or interpreted

		MetBound(MetabolitePtr ptr, bool isMinBound)
			: _met(ptr), _isMinBound(isMinBound) {
			_bound = -INFINITY;
			_isBoundary = false;
		}

		// gives the value of the bound in the original setting
		double orig_value() const {
			if(_isMinBound) {
				return -_bound;
			}
			else {
				return _bound;
			}
		}

		// returns the maximal value this bound can attain
		double getMax() const {
			if(_isMinBound) {
				return -_met->getPotLb();
			}
			else {
				return _met->getPotUb();
			}
		}
	};

	typedef boost::shared_ptr<MetBound> MetBoundPtr;

	std::size_t hash_value(MetBoundPtr const & mb );

	struct Arc {
		MetBoundPtr _target;
		ReactionPtr _creator;
		bool _fwdcreator; // is set to true, iff the arc was created by the forward part of the _creator
		int _num_used;

		//TODO:
		bool _active; // indicates, if this arc can be active

		Arc();

		std::vector<std::pair<MetBoundPtr,double> > _input;

		int unknown_inputs() const;

		double update_value() const;
	};

	typedef boost::shared_ptr<Arc> ArcPtr;

	struct ArcEvent {
		ArcPtr arc;
		int unknown_inputs;
		double gain_val;
		bool first;

		ArcEvent(ArcPtr a) : arc(a) {
			unknown_inputs = a->unknown_inputs();
			gain_val = a->_target->getMax() - a->update_value(); // Why not use the current value of the bound?
			//gain_val = a->_target->getMax() - ldexp(a->update_value(),a->_num_used);
			//gain_val = ldexp(a->update_value() - a->_target->_bound, a->_num_used);
			first = isinf(a->_target->_bound);
		}

		/**
		 * returns true, if the other event should be executed earlier
		 */
		bool operator<(const ArcEvent a) const {
			int diff_inputs = unknown_inputs - a.unknown_inputs;
			if(diff_inputs > 0) {
				return true;
			}
			else if(diff_inputs < 0) {
				return false;
			}
			else if(a.first) {
				return true;
			}
			else if(first) {
				return false;
			}
			else {
				// biggest update corresponds to smallest gain_val
				if(unknown_inputs == 0) {
					// do the biggest update (we have to do it anyway)
					return gain_val > a.gain_val;
				}
				else {
					// do the smallest update (we may still improve)
					return gain_val < a.gain_val;
				}

			}
		}
	};

	const ModelPtr _model;

	boost::unordered_map<MetabolitePtr, MetBoundPtr> _minBounds;
	boost::unordered_map<MetabolitePtr, MetBoundPtr> _maxBounds;

	boost::unordered_map<MetBoundPtr, std::vector<ArcPtr> > _incidence; // use separate map to avoid cyclic dependencies

	std::vector<ArcPtr> _arcs;

	std::priority_queue<ArcEvent> _queue;

	// for the update-lps, the variable indexing is always the same, so we store this globally
	boost::unordered_map<MetBoundPtr, int> _vars;
	typedef std::pair<MetBoundPtr, int> VarEntry;

	// creates arc for reversed reaction
	ArcPtr getReversed(const ArcPtr a) const;

	/// initially fills queue
	void init_Queue();

	void update(ArcPtr a);

	/**
	 * performs one update-step w.r.t. the inflow/outflow property.
	 * Returns true if the update step changed something
	 */
	bool updateStepFlow();

	/**
	 * performs one update-step w.r.t. hard constraints (like bounds or branching-decisions).
	 * Returns true if the update step changed something
	 */
	void updateStepHard(ScipModelPtr scip);

	/**
	 * builds constraint for fixed direction of reaction
	 */
	void buildHardArcConstraint(ArcPtr a, SCIP_LPI* lpi, bool fwd);
};

} /* namespace metaopt */
#endif /* POTBOUNDPROPAGATION2_H_ */
