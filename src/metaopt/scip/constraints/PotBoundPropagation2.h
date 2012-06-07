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
#include "metaopt/model/Model.h"
#include "metaopt/model/scip/ScipModel.h"

namespace metaopt {

class PotBoundPropagation2 {
public:
	PotBoundPropagation2(ModelPtr model);
	virtual ~PotBoundPropagation2();

	void print();

	/**
	 * updates the computed bound to be a good mu-bound.
	 * Uses the fixed directions of current setup.
	 * If we cannot be sure if the current bound is a proper mu-bound, computation is started from scratch
	 */
	void update(ScipModelPtr scip);

	// from the computed potential bounds lists all blocked reactions
	boost::shared_ptr<std::vector<ReactionPtr> > getBlockedReactions();

private:
	struct Arc;

	struct MetBound {
		MetabolitePtr _met;
		bool _isMinBound;
		// value of the bound in transformed setting
		double _bound;
		boost::weak_ptr<Arc> _last_update; // only logging information, so the results can be verified and/or interpreted

		MetBound(MetabolitePtr ptr, bool isMinBound)
			: _met(ptr), _isMinBound(isMinBound) {
			_bound = -INFINITY;
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

	struct Arc {
		MetBoundPtr _target;
		ReactionPtr _creator;

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
			gain_val = a->_target->getMax() - a->update_value();
			first = isinf(a->_target->_bound);
		}

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
				// there are no better arcs that update infinity values,
				// so we can do the biggest update, since entries that are still infinity, stay infinity
				return gain_val < a.gain_val;
			}
		}
	};

	const ModelPtr _model;

	boost::unordered_map<MetabolitePtr, MetBoundPtr> _minBounds;
	boost::unordered_map<MetabolitePtr, MetBoundPtr> _maxBounds;

	boost::unordered_map<MetBoundPtr, std::vector<ArcPtr> > _incidence; // use separate map to avoid cyclic dependencies

	std::vector<ArcPtr> _arcs;

	std::priority_queue<ArcEvent> _queue;

	/// initially fills queue
	void init_Queue();

	void update(ArcPtr a);

	/**
	 * performs one update-step w.r.t. the inflow/outflow property.
	 * Returns true if the update step changed something
	 */
	bool updateStepFlow(ScipModelPtr scip);

	/**
	 * performs one update-step w.r.t. hard constraints (like bounds or branching-decisions).
	 * Returns true if the update step changed something
	 */
	bool updateStepHard(ScipModelPtr scip);
};

} /* namespace metaopt */
#endif /* POTBOUNDPROPAGATION2_H_ */
