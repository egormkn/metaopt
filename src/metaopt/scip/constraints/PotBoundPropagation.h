/*
 * PotBoundPropagation.h
 *
 *  Created on: 18.05.2012
 *      Author: arnem
 */

#ifndef POTBOUNDPROPAGATION_H_
#define POTBOUNDPROPAGATION_H_

#include <boost/unordered_map.hpp>
#include <vector>
#include <queue>
#include "metaopt/model/Model.h"

namespace metaopt {

class PotBoundPropagation {
public:
	PotBoundPropagation(ModelPtr model);
	virtual ~PotBoundPropagation();

	void print();

private:
	struct MetBound {
		MetabolitePtr _met;
		bool _isMinBound;
		double _bound;

		boost::shared_ptr<boost::unordered_set<MetabolitePtr> > _unused;
		boost::unordered_map<ReactionPtr, boost::shared_ptr<boost::unordered_set<MetabolitePtr> > > _unusedFrom;

		MetBound(MetabolitePtr ptr, bool isMinBound);

		// how much more than the default bound is this MetBound restrictive?
		double restrictivity() const;

		void updateUnused();
	};

	typedef boost::shared_ptr<MetBound> MetBoundPtr;

	struct MetBoundCompare {
		bool operator()(const MetBoundPtr a, const MetBoundPtr b) {
			return a->restrictivity() < b->restrictivity();
		}
	};

	const ModelPtr _model;

	boost::unordered_map<MetabolitePtr, MetBoundPtr> _minBounds;
	boost::unordered_map<MetabolitePtr, MetBoundPtr> _maxBounds;

	std::priority_queue<MetBoundPtr, std::vector<MetBoundPtr>, MetBoundCompare> _queue;
	boost::unordered_set<MetBoundPtr> _inQueue;

	/// initially fills queue
	void init_Queue();

	void update(MetBoundPtr b);
	void updateProducts(ReactionPtr& r, MetabolitePtr& ignore, boost::shared_ptr<boost::unordered_set<MetabolitePtr> >& unused, double& val);
	void updateReactants(ReactionPtr& r, MetabolitePtr& ignore, boost::shared_ptr<boost::unordered_set<MetabolitePtr> >& unused, double& val);

	inline void enqueue(MetBoundPtr m);
	inline MetBoundPtr pop();
};

inline void PotBoundPropagation::enqueue(MetBoundPtr m) {
	if(_inQueue.find(m) == _inQueue.end()) {
		_queue.push(m);
		_inQueue.insert(m);
	}
}

inline PotBoundPropagation::MetBoundPtr PotBoundPropagation::pop() {
	PotBoundPropagation::MetBoundPtr out = _queue.top();
	_queue.pop();
	_inQueue.erase(out);
	return out;
}

} /* namespace metaopt */
#endif /* POTBOUNDPROPAGATION_H_ */
