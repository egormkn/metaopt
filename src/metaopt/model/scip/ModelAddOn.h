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

	/** returns the name of this Addon. Useful for debugging */
	inline const std::string& getName() const;

	/**
	 * Helper method for ScipModel destruction process.
	 * This method should only be called from the destructor of ScipModel.
	 * It processes the destruction of this ModelAddOn that still requires data of ScipModel.
	 * A reference to ScipModel is supplied, since weak_ptrs are already expired during destruction.
	 * After this method has been called, the indicator destroyed is set.
	 *
	 * Important: If this addon has dependencies on other addons, these dependencies must be destroyed, too.
	 */
	virtual void destroy(ScipModel* model);

protected:
	inline ScipModelPtr getModel() const;
	inline bool isDestroyed() const;

private:
	boost::weak_ptr<ScipModel> _model;
	std::string _name;
	bool _destroyed; // indicates if destroy(ScipModel&) has been called
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

inline bool ModelAddOn::isDestroyed() const {
	return _destroyed;
}

} /* namespace metaopt */
#endif /* MODELADDON_H_ */
