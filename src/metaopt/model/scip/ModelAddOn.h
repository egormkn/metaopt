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
#include "metaopt/Uncopyable.h"
#include "metaopt/Properties.h"

namespace metaopt {

/**
 * ModelAddOnS are used to add variables that are needed for certain formulation to the basic ScipModel.
 *
 * Since implementations of Heuristics have to set values for all variables,
 * each add on has to supply a method computeSolutionVals that tries to compute solution values
 * for the additional variables from the basic variables of the ScipModel or other add ons.
 */
class ModelAddOn : Uncopyable {
public:
	ModelAddOn(ScipModelPtr model, std::string name);
	virtual ~ModelAddOn();

	/**
	 * Tries to compute solution values for the given solution.
	 *
	 * @returns true, if successfully computes solution values.
	 */
	virtual bool computeSolutionVals(SolutionPtr sol) const = 0;

	/**
	 * Allows the addon to perform additional fixations according to branching decisions.
	 *
	 * For a complete explanation, see ScipModel#setDirection(SCIP_NODE*, ReactionPtr, bool)
	 *
	 * If this method is not implemented, the addon does nothing.
	 */
	virtual void setDirection(SCIP_NODE* node, ReactionPtr rxn, bool fwd);

	/**
	 * Fetches the list of fixed reactions the addon knows about.
	 * The default implementation simply returns an empty list.
	 */
	virtual boost::shared_ptr<const boost::unordered_set<ReactionPtr> > getFixedDirections();

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
