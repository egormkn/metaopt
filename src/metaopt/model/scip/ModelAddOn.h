/*
 * ModelAddOn.h
 *
 *  Created on: 18.04.2012
 *      Author: arnem
 */

#ifndef MODELADDON_H_
#define MODELADDON_H_

#include <boost/enable_shared_from_this.hpp>
#include "metaopt/model/scip/ScipModel.h"
#include "Solution.h"
#include "metaopt/scip/ScipError.h"

namespace metaopt {

/**
 * ModelAddOnS are used to add variables that are needed for certain formulation to the basic ScipModel.
 *
 * Since implementations of Heuristics have to set values for all variables,
 * each add on has to supply a method computeSolutionVals that tries to compute solution values
 * for the additional variables from the basic variables of the ScipModel or other add ons.
 */
class ModelAddOn {
public:
	ModelAddOn(ScipModelPtr model, std::string name);
	virtual ~ModelAddOn();

	/**
	 * Tries to compute solution values for the given solution.
	 *
	 * @returns true, if successfully computes solution values.
	 */
	virtual bool computeSolutionVals(SolutionPtr sol) const = 0;

	inline const std::string& getName() const;

protected:
	inline ScipModelPtr getModel() const;

private:
	boost::weak_ptr<ScipModel> _model;
	std::string _name;
};

typedef boost::error_info<struct tag_addon_name,std::string> addon_name;

inline const std::string& ModelAddOn::getName() const {
	return _name;
}

inline ScipModelPtr ModelAddOn::getModel() const {
	if(_model.expired()) {
		BOOST_THROW_EXCEPTION(ModelOwnershipError() << addon_name(getName()));
	}
	return _model.lock();
}

} /* namespace metaopt */
#endif /* MODELADDON_H_ */
