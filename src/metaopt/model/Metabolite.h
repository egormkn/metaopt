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

namespace metaopt {

class Model;

class Metabolite : public boost::enable_shared_from_this<Metabolite>
{
public:
	Metabolite(boost::weak_ptr<Model> model, std::string name);
	virtual ~Metabolite();


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



	/** @brief Fetches the name of this metabolite.
	 *
	 * @return the name of this metabolite.
	 */
	std::string getName() const;


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


private:
	boost::weak_ptr<Model> _model; /** reference to the Model owning this Metabolite */

	double _lb;
	double _ub;
	std::string _name;
	bool _boundary;

	void notifyChange(); /** notifies the Model of changes performed on this Metabolite. */

};


/** Please use MetabolitePtr to handle pointers to Metabolite instances */
typedef boost::shared_ptr<Metabolite> MetabolitePtr;

/** If an error is thrown by a method of this Metabolite, the name of this Metabolite is added using this tag */
typedef boost::error_info<struct tag_metabolite_name,std::string> metabolite_name;

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
}

inline void Metabolite::setPotUb(double ub) {
	_ub = ub;
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

} /* namespace metaopt */
#endif /* METABOLITE_H_ */
