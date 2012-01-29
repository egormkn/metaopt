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

namespace metaopt {

class Model;

class Metabolite : public boost::enable_shared_from_this<Metabolite>
{
public:
	Metabolite();
	virtual ~Metabolite();


	/** @brief Gets the lower bound on the potential of this metabolite.
	 *
	 * @return the lower bound.
	 */
	virtual double getPotLb() const = 0;

	/** @brief Gets the upper bound on the potential of this metabolite.
	 *
	 * @return the upper bound.
	 */
	virtual double getPotUb() const = 0;

	/** @brief Sets the lower bound on the potential of this metabolite.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the lower bound.
	 */
	virtual void setPotLb(double lb) = 0;

	/** @brief Sets the upper bound on the potential of this metabolite.
	 *
	 * Implementations shall call notifyChange() to update the Model of the change.
	 * @param the upper bound.
	 */
	virtual void setPotUb(double ub) = 0;



	/** @brief Fetches the name of this metabolite.
	 *
	 * @return the name of this metabolite.
	 */
	virtual std::string getName() const = 0;


	/** @brief tests if this metabolite is isolated.
	 *
	 * A metabolite is isolated, if there exists no reaction in the model that contains it.
	 */
	bool isIsolated() const;



	/** returns the owning model */
	const boost::shared_ptr<Model> getOwner() const;


private:
	boost::weak_ptr<Model> _model; /** reference to the Model owning this Metabolite */
	void notifyChange(); /** notifies the Model of changes performed on this Metabolite. */

};

/** Please use MetabolitePtr to handle pointers to Metabolite instances */
typedef boost::shared_ptr<Metabolite> MetabolitePtr;

/** If an error is thrown by a method of this Metabolite, the name of this Metabolite is added using this tag */
typedef boost::error_info<struct tag_metabolite_name,std::string> metabolite_name;


} /* namespace metaopt */
#endif /* METABOLITE_H_ */
