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
 * Metabolite.h
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#ifndef METABOLITE_H_
#define METABOLITE_H_

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <string>
#include <boost/exception/error_info.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include "metaopt/model/Precision.h"
#include "metaopt/Uncopyable.h"
#include "metaopt/Properties.h"

namespace metaopt {

class Model;
class Reaction;
typedef boost::shared_ptr<Reaction> ReactionPtr;

class Metabolite : public boost::enable_shared_from_this<Metabolite>, Uncopyable
{
public:
	Metabolite(boost::weak_ptr<Model> model, std::string name);
	virtual ~Metabolite();


	/**
	 * tests if the potential of this metabolite participates in objective functions
	 */
	inline bool isObjective() const;

	/** @brief Gets the lower bound on the potential of this metabolite.
	 *
	 * @return the lower bound.
	 */
	double getPotLb() const;

	/** @brief Gets the upper bound on the potential of this metabolite.
	 *
	 * @return the upper bound.
	 */
	double getPotUb() const;

	/** @brief Sets the lower bound on the potential of this metabolite.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the lower bound.
	 */
	void setPotLb(double lb);

	/** @brief Sets the upper bound on the potential of this metabolite.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the upper bound.
	 */
	void setPotUb(double ub);


	/** @brief Gets the objective on the potential of this metabolite.
	 *
	 * @return the objective.
	 */
	double getPotObj() const;

	/** @brief Sets the objective on the potential of this metabolite.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the lower bound.
	 */
	void setPotObj(double obj);


	/** @brief Fetches the name of this reaction.
	 * Attention: Returns a copy of the name!!!
	 * If you access the data of the string, make sure it still exists!
	 *
	 * @return the name of this reaction.
	 */
	std::string getName() const;

	/** @brief Fetches the name of this reaction.
	 * This does the same as getName(), but it does not copy the string.
	 * It only returns a pointer to the internal data.
	 * The pointer will be valid as long as this reaction lives.
	 *
	 * @return the name of this reaction.
	 */
	const char* getCName() const;

	/** @brief tests if this metabolite is isolated.
	 *
	 * A metabolite is isolated, if there exists no reaction in the model that contains it.
	 */
	bool isIsolated() const;

	/**
	 * Does the steady-state assumption have to hold for the metabolite?
	 * If it has the boundary condition, it needs not.
	 */
	bool hasBoundaryCondition() const;

	/** Sets the boundary condition */
	void setBoundaryCondition(bool boundary);

	/** returns the owning model */
	const boost::shared_ptr<Model> getOwner() const;

	/**
	 * lists all reactions producing this metabolite
	 * The list is created anew from a stored list of weak pointers.
	 */
	boost::shared_ptr<std::vector<ReactionPtr> > getProducers() const;
	/**
	 * lists all reactions producing this metabolite
	 * The list is created anew from a stored list of weak pointers.
	 */
	boost::shared_ptr<std::vector<ReactionPtr> > getConsumers() const;

	/** do not call this directly, but use the methods provided by Reaction */
	void removeProducer(ReactionPtr r);
	/** do not call this directly, but use the methods provided by Reaction */
	void removeConsumer(ReactionPtr r);
	/** do not call this directly, but use the methods provided by Reaction */
	void addProducer(ReactionPtr r);
	/** do not call this directly, but use the methods provided by Reaction */
	void addConsumer(ReactionPtr r);

	/**
	 * updates the potential precision (usually only called from Model).
	 */
	void setPotPrecision(PrecisionPtr precision);

private:
	boost::weak_ptr<Model> _model; /** reference to the Model owning this Metabolite */

	double _lb;
	double _ub;
	double _obj;
	std::string _name;
	bool _boundary;

	PrecisionPtr _potPrecision;

	std::vector<boost::weak_ptr<Reaction> > _producers;
	std::vector<boost::weak_ptr<Reaction> > _consumers;

	void notifyChange(); /** notifies the Model of changes performed on this Metabolite. */

};


/** Please use MetabolitePtr to handle pointers to Metabolite instances */
typedef boost::shared_ptr<Metabolite> MetabolitePtr;

/** If an error is thrown by a method of this Metabolite, the name of this Metabolite is added using this tag */
typedef boost::error_info<struct tag_metabolite_name,std::string> metabolite_name;

std::size_t hash_value(MetabolitePtr const & met );

/////////////////////////////////////
// Inline function defs
/////////////////////////////////////

inline double Metabolite::getPotLb() const {
	return _lb;
}

inline double Metabolite::getPotUb() const {
	return _ub;
}

inline void Metabolite::setPotLb(double lb) {
	_lb = lb;
	notifyChange();
}

inline void Metabolite::setPotUb(double ub) {
	_ub = ub;
	notifyChange();
}

inline double Metabolite::getPotObj() const {
	return _obj;
}

inline void Metabolite::setPotObj(double obj) {
	_obj = obj;
	notifyChange();
}

inline bool Metabolite::hasBoundaryCondition() const {
	return _boundary;
}

inline void Metabolite::setBoundaryCondition(bool boundary) {
	_boundary = boundary;
}

inline std::string Metabolite::getName() const {
  return _name;
}

inline const char* Metabolite::getCName() const {
  return _name.c_str();
}

inline bool Metabolite::isObjective() const {
	return(_obj < -_potPrecision->getCheckTol() || _obj > _potPrecision->getCheckTol());
}

} /* namespace metaopt */
#endif /* METABOLITE_H_ */
