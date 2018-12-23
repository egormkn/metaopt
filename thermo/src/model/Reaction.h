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
 * Reaction.h
 *
 *  Created on: 28.01.2012
 *      Author: arne
 */

#ifndef REACTION_H_
#define REACTION_H_

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <string>
#include <boost/exception/all.hpp>
#include <boost/unordered_map.hpp>

#include "Metabolite.h"
#include "Precision.h"
#include "Uncopyable.h"
#include "Properties.h"

namespace metaopt {

class Model;

/** @brief Abstract class of a reaction in a metabolic network.
 *
 * The transformation speed reactants into products of a reaction is described by the flux through the reaction.
 * The flux can be bounded by an upper and lower bound. When we optimize on fluxes, each reaction has an objective coefficient (which may be 0).
 * For debugging purposes each reaction also has a name. This name must not be unique, though (although uniqueness will help debugging a lot!).
 *
 * Each Reaction belongs to a Model. Hence, reactions should be created using the factory method Model::createReaction().
 * Usually, however, the model instantiation already creates the reaction. So there is no need for the user to create reactions manually.
 *
 * This class abstractely represents reactions with these properties. Implementations may use direct implementations or SBML-models as a data-backend.
 */
class Reaction : public boost::enable_shared_from_this<Reaction>, Uncopyable
{
public:
	/** @brief Creates a new Reaction. Please do not use this, but use Model::createReaction() instead.
	 *
	 * @param model a weak pointer to the owning model. The pointer is weak to prevent loops.
	 */
	Reaction(boost::weak_ptr<Model> model, std::string name);
	virtual ~Reaction();

	/**
	 * true if positive and negative flux are possible
	 */
	inline bool isReversible() const;

	/**
	 * true if positive flux is possible
	 */
	inline bool canFwd() const;

	/**
	 * true if only positive flux is possible
	 */
	inline bool isFwdForcing() const;

	/**
	 * true if only negative flux is possible
	 */
	inline bool isBwdForcing() const;

	/**
	 * true if negative flux is possible
	 */
	inline bool canBwd() const;

	/**
	 * true if only zero flux is possible
	 */
	inline bool isBlocked() const;

	/** @brief Gets the lower bound on flux through this reaction.
	 *
	 * @return the lower bound.
	 */
	inline double getLb() const;

	/** @brief Gets the upper bound on flux through this reaction.
	 *
	 * @return the upper bound.
	 */
	inline double getUb() const;

	/** @brief Gets the objective coefficient of the flux through this reaction.
	 *
	 * @return the objective coefficient.
	 */
	inline double getObj() const;

	/** @brief Sets the lower bound on flux through this reaction.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the lower bound.
	 */
	inline void setLb(double lb);

	/** @brief Sets the upper bound on flux through this reaction.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the upper bound.
	 */
	inline void setUb(double ub);

	/** @brief Sets the objective coefficient of the flux through this reaction.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the objective coefficient.
	 */
	inline void setObj(double obj);



	/** @brief Fetches the name of this reaction.
	 * Attention: Returns a copy of the name!!!
	 * If you access the data of the string, make sure it still exists!
	 *
	 * @return the name of this reaction.
	 */
	inline std::string getName() const;

	/** @brief Fetches the name of this reaction.
	 * This does the same as getName(), but it does not copy the string.
	 * It only returns a pointer to the internal data.
	 * The pointer will be valid as long as this reaction lives.
	 *
	 * @return the name of this reaction.
	 */
	inline const char* getCName() const;

	/** true, if this reaction is an exchange reaction */
	inline bool isExchange() const;

	inline void setExchange(bool exchange);

	/**
	 * checks if this reaction is problematic.
	 *
	 * Problematic reactions are reactions whose flux variable is also used in additional constraint (e.g. kinetics).
	 * It is important to mark such reactions as problematic to improve the branching decisions of the constraint handler.
	 */
	inline bool isProblematic() const;

	/**
	 * checks if this reaction is flux forcing.
	 *
	 * A reaction is flux forcing, if it has a positive lower bound or a negative upper bound
	 */
	inline bool isFluxForcing() const;

	/**
	 * checks if this reaction is objective.
	 *
	 * A reaction is objectvie, if it has a non-zero objective coefficient
	 */
	inline bool isObjective() const;

	/**
	 * marks this reaction as problematic.
	 *
	 * Problematic reactions are reactions whose flux variable is also used in additional constraint (e.g. kinetics).
	 * It is important to mark such reactions as problematic to improve the branching decisions of the constraint handler.
	 */
	inline void setProblematic(bool problematic);

	/** @brief returns a map of the metabolites involved in this reaction with their stoichiometric coefficients.
	 *
	 * Reactants have a negative coefficient. Products have a positive coefficient.
	 * Metabolites that appear both as reactants and products are combined in one coefficient.
	 * The list must contain the union of all reactants and products. Hence, it may contain zero entries if the stoichiometries cancel.
	 */
	const boost::unordered_map<MetabolitePtr, double>& getStoichiometries() const;

	/** @brief returns a list of all reactants of this reaction.
	 *
	 * All associated values are positive.
	 */
	const boost::unordered_map<MetabolitePtr, double>& getReactants() const;

	/** @brief returns a list of all products of this reaction.
	 *
	 * All associated values are positive
	 */
	const boost::unordered_map<MetabolitePtr, double>& getProducts() const;

	/** The stoichiometry of the given metabolite in this reaction.
	 *
	 * If the metabolite does not participate in this reaction, the returned value is 0.
	 */
	double getStoichiometry(MetabolitePtr m) const;
	/** The amount the given metabolite is consumed by this reaction.
	 *
	 * If the metabolite is not a reactant of this reaction, the returned value is 0.
	 */
	double getReactant(MetabolitePtr m) const;
	/** The amount the given metabolite is produced by this reaction.
	 *
	 * If the metabolite is not produced by this reaction, the returned value is 0.
	 */
	double getProduct(MetabolitePtr m) const;

	/** Sets the stoichiometry of this metabolite in this reaction.
	 *
	 * If the metabolite was not part of this reaction, it is added.
	 * If the stoichiometry is positive, it is set as a product.
	 * If the stoichiometry is negative, it is set as a reactant.
	 * If the stoichiometry is zero, it will be removed.
	 *
	 * Note, that the metabolite will only be configured either as reactant or product, even if it was previously registered as reactant or product.
	 */
	void setStoichiometry(MetabolitePtr m, double value);

	/** Sets the metabolite as a reactant.
	 *
	 * If the metabolite was not part of this reaction, it is added.
	 * If it was already a reactant, only the stiochiometric value is updated.
	 * If it was already a product, the product value will be retained. This way it is possible to add metabolites that are both consumed and produced by the reaction.
	 * If the value is zero and the product value is also zero, then the metabolite will be removed from this reaction.
	 *
	 * Only non-negative values are allowed
	 */
	void setReactant(MetabolitePtr m, double value);

	/** Sets the metabolite as a reactant.
	 *
	 * If the metabolite was not part of this reaction, it is added.
	 * If it was already a product, only the stiochiometric value is updated.
	 * If it was already a reactant, the reactant value will be retained. This way it is possible to add metabolites that are both consumed and produced by the reaction.
	 * If the value is zero and the reactant value is also zero, then the metabolite will be removed from this reaction.

	 * Only non-negative values are allowed
	 */
	void setProduct(MetabolitePtr m, double value);

	/**
	 * updates the flux precision (usually only called from Model).
	 */
	void setFluxPrecision(PrecisionPtr precision);

	/** returns the owning model */
	const boost::shared_ptr<Model> getOwner() const;

	/** returns a human readable description of this reaction including stoichiometries */
	std::string toString() const;

protected:
	/** empty constructor for ReversedReaction */
	Reaction();

private:
	boost::weak_ptr<Model> _model; /** reference to the Model owning this Reaction */

	boost::unordered_map<MetabolitePtr, double> _stoichiometries; /** list of all metabolites participating in this reaction */
	boost::unordered_map<MetabolitePtr, double> _reactants; /** list of all reactants */
	boost::unordered_map<MetabolitePtr, double> _products; /** list of all products */

	std::string _name;
	double _lb;
	double _ub;
	double _obj;
	bool _exchange;
	bool _problematic;

	PrecisionPtr _fluxPrecision;

	void notifyChange(); /** notifies the Model of changes performed on this Reaction. */
};

/** Please use ReactionPtr to handle pointers to Reaction instances */
typedef boost::shared_ptr<Reaction> ReactionPtr;

/** Used if an reaction is not found */
struct UnknownReactionError : virtual boost::exception, virtual std::exception {};

/** If an error is thrown by a method of this Reaction, the name of this Reaction is added using this tag */
typedef boost::error_info<struct tag_reaction_name,std::string> reaction_name;

/** Typedef for easy use of foreach iteration through stoichiometry data */
typedef std::pair<const MetabolitePtr, double> Stoichiometry ;

std::size_t hash_value(ReactionPtr const & rxn );


/////////////////////////////////////
// Inline function defs
/////////////////////////////////////

double Reaction::getLb() const {
	return _lb;
}

double Reaction::getUb() const {
	return _ub;
}

double Reaction::getObj() const {
	return _obj;
}

void Reaction::setLb(double lb) {
	_lb = lb;
	notifyChange();
}

void Reaction::setUb(double ub) {
	_ub = ub;
	notifyChange();
}

void Reaction::setObj(double obj) {
	_obj = obj;
	notifyChange();
}

std::string Reaction::getName() const {
	return _name;
}

const char* Reaction::getCName() const {
	return _name.c_str();
}

bool Reaction::isExchange() const {
	return _exchange;
}

void Reaction::setExchange(bool exchange) {
	_exchange = exchange;
	notifyChange();
}

bool Reaction::isProblematic() const {
	return _problematic;
}

void Reaction::setProblematic(bool problematic) {
	_problematic = problematic;
	notifyChange();
}

bool Reaction::isObjective() const {
	return _obj < -_fluxPrecision->getCheckTol() || _obj > _fluxPrecision->getCheckTol();
}

bool Reaction::isFluxForcing() const {
	return _lb > _fluxPrecision->getCheckTol() || _ub < -_fluxPrecision->getCheckTol();
}

bool Reaction::isReversible() const {
	return _lb < -_fluxPrecision->getCheckTol() && _ub > _fluxPrecision->getCheckTol();
}

bool Reaction::canFwd() const {
	return _ub > _fluxPrecision->getCheckTol();
}

bool Reaction::canBwd() const {
	return _lb < -_fluxPrecision->getCheckTol();
}

bool Reaction::isFwdForcing() const {
	return _lb > _fluxPrecision->getCheckTol();
}

bool Reaction::isBwdForcing() const {
	return _ub < -_fluxPrecision->getCheckTol();
}

bool Reaction::isBlocked() const {
	return _lb > -_fluxPrecision->getCheckTol() && _ub < _fluxPrecision->getCheckTol();
}



} /* namespace metaopt */


#endif /* REACTION_H_ */
