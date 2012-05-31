/*
 * ReversedReaction.h
 *
 *  Created on: 18.05.2012
 *      Author: arnem
 */

#ifndef REVERSEDREACTION_H_
#define REVERSEDREACTION_H_

#include "Reaction.h"

namespace metaopt {

class ReversedReaction: public metaopt::Reaction {
public:
	ReversedReaction(weak_ptr<Reaction> original);
	virtual ~ReversedReaction();


	virtual bool isReversed() const;

	/** @brief Gets the lower bound on flux through this reaction.
	 *
	 * @return the lower bound.
	 */
	virtual double getLb() const;

	/** @brief Gets the upper bound on flux through this reaction.
	 *
	 * @return the upper bound.
	 */
	virtual double getUb() const;

	/** @brief Gets the objective coefficient of the flux through this reaction.
	 *
	 * @return the objective coefficient.
	 */
	virtual double getObj() const;

	/** @brief Sets the lower bound on flux through this reaction.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the lower bound.
	 */
	virtual void setLb(double lb);

	/** @brief Sets the upper bound on flux through this reaction.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the upper bound.
	 */
	virtual void setUb(double ub);

	/** @brief Sets the objective coefficient of the flux through this reaction.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the objective coefficient.
	 */
	virtual void setObj(double obj);



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

	/** true, if this reaction is an exchange reaction */
	bool isExchange() const;

	void setExchange(bool exchange);

	/** @brief returns a map of the metabolites involved in this reaction with their stoichiometric coefficients.
	 *
	 * Reactants have a negative coefficient. Products have a positive coefficient.
	 * Metabolites that appear both as reactants and products are combined in one coefficient.
	 * The list must contain the union of all reactants and products. Hence, it may contain zero entries if the stoichiometries cancel.
	 */
	virtual const boost::unordered_map<MetabolitePtr, double>& getStoichiometries() const;

	/** @brief returns a list of all reactants of this reaction.
	 *
	 * All associated values are positive.
	 */
	virtual const boost::unordered_map<MetabolitePtr, double>& getReactants() const;

	/** @brief returns a list of all products of this reaction.
	 *
	 * All associated values are positive
	 */
	virtual const boost::unordered_map<MetabolitePtr, double>& getProducts() const;

	/** The stoichiometry of the given metabolite in this reaction.
	 *
	 * If the metabolite does not participate in this reaction, the returned value is 0.
	 */
	virtual double getStoichiometry(MetabolitePtr m) const;
	/** The amount the given metabolite is consumed by this reaction.
	 *
	 * If the metabolite is not a reactant of this reaction, the returned value is 0.
	 */
	virtual double getReactant(MetabolitePtr m) const;
	/** The amount the given metabolite is produced by this reaction.
	 *
	 * If the metabolite is not produced by this reaction, the returned value is 0.
	 */
	virtual double getProduct(MetabolitePtr m) const;

	/** Sets the stoichiometry of this metabolite in this reaction.
	 *
	 * If the metabolite was not part of this reaction, it is added.
	 * If the stoichiometry is positive, it is set as a product.
	 * If the stoichiometry is negative, it is set as a reactant.
	 * If the stoichiometry is zero, it will be removed.
	 *
	 * Note, that the metabolite will only be configured either as reactant or product, even if it was previously registered as reactant or product.
	 */
	virtual void setStoichiometry(MetabolitePtr m, double value);

	/** Sets the metabolite as a reactant.
	 *
	 * If the metabolite was not part of this reaction, it is added.
	 * If it was already a reactant, only the stiochiometric value is updated.
	 * If it was already a product, the product value will be retained. This way it is possible to add metabolites that are both consumed and produced by the reaction.
	 * If the value is zero and the product value is also zero, then the metabolite will be removed from this reaction.
	 *
	 * Only non-negative values are allowed
	 */
	virtual void setReactant(MetabolitePtr m, double value);

	/** Sets the metabolite as a reactant.
	 *
	 * If the metabolite was not part of this reaction, it is added.
	 * If it was already a product, only the stiochiometric value is updated.
	 * If it was already a reactant, the reactant value will be retained. This way it is possible to add metabolites that are both consumed and produced by the reaction.
	 * If the value is zero and the reactant value is also zero, then the metabolite will be removed from this reaction.

	 * Only non-negative values are allowed
	 */
	virtual void setProduct(MetabolitePtr m, double value);

	/** @brief Fetches the name of this reaction.
	 * Attention: Returns a copy of the name!!!
	 * If you access the data of the string, make sure it still exists!
	 *
	 * @return the name of this reaction.
	 */
	virtual std::string getName() const;

	/** @brief Fetches the name of this reaction.
	 * This does the same as getName(), but it does not copy the string.
	 * It only returns a pointer to the internal data.
	 * The pointer will be valid as long as this reaction lives.
	 *
	 * @return the name of this reaction.
	 */
	virtual const char* getCName() const;

	/** true, if this reaction is an exchange reaction */
	virtual bool isExchange() const;

	virtual void setExchange(bool exchange);

	/** returns the owning model */
	virtual const boost::shared_ptr<Model> getOwner() const;

	/** returns a human readable description of this reaction including stoichiometries */
	virtual std::string toString() const;

private:
	boost::weak_ptr<Reaction> _original;

	inline ReactionPtr getOrig() const {
		return _original.lock();
	}
};

} /* namespace metaopt */
#endif /* REVERSEDREACTION_H_ */
