/*
 * ReversedReaction.cpp
 *
 *  Created on: 18.05.2012
 *      Author: arnem
 */

#include "ReversedReaction.h"
#include "metaopt/Properties.h"

using namespace boost;

namespace metaopt {

ReversedReaction::ReversedReaction(weak_ptr<Reaction> original) :
	Reaction::Reaction(), _original(original)
{

}

ReversedReaction::~ReversedReaction() {
	// TODO Auto-generated destructor stub
}

double ReversedReaction::getLb() const {
	return -getOrig()->getUb();
}

double ReversedReaction::getUb() const {
	return -getOrig()->getLb();
}

double ReversedReaction::getObj() const {
	return -getOrig()->getObj();
}

void ReversedReaction::setLb(double lb) {
	getOrig()->setUb(-lb);
}

void ReversedReaction::setUb(double ub) {
	getOrig()->setLb(-ub);
}

void ReversedReaction::setObj(double obj) {
	getOrig()->setObj(-obj);
}


const boost::unordered_map<MetabolitePtr, double>& ReversedReaction::getStoichiometries() const {
	boost::unordered_map<MetabolitePtr, double> out;
	// TODO: this is completely inefficient!!!
	foreach(const Stoichiometry s, getOrig()->getStoichiometries) {
		out[s.first] = -s.second;
	}
	return out;
}


const boost::unordered_map<MetabolitePtr, double>& ReversedReaction::getReactants() const {
	boost::unordered_map<MetabolitePtr, double> out;
	// TODO: this is completely inefficient!!!
	foreach(const Stoichiometry s, getOrig()->getProducts) {
		out[s.first] = -s.second;
	}
	return out;
}


const boost::unordered_map<MetabolitePtr, double>& ReversedReaction::getProducts() const {
	boost::unordered_map<MetabolitePtr, double> out;
	// TODO: this is completely inefficient!!!
	foreach(const Stoichiometry s, getOrig()->getReactants) {
		out[s.first] = -s.second;
	}
	return out;
}


double ReversedReaction::getStoichiometry(MetabolitePtr m) const {
	return -getOrig()->getStoichiometry(m);
}

double ReversedReaction::getReactant(MetabolitePtr m) const {
	return getOrig()->getProduct(m);
}

double ReversedReaction::getProduct(MetabolitePtr m) const {
	return getOrig()->getReactant(m);
}


void ReversedReaction::setStoichiometry(MetabolitePtr m, double value) {
	getOrig()->setStoichiometry(m, -value);
}


void ReversedReaction::setReactant(MetabolitePtr m, double value) {
	getOrig()->setProduct(m, value);
}


void ReversedReaction::setProduct(MetabolitePtr m, double value) {
	getOrig()->setReactant(m, value);
}

bool ReversedReaction::isReversed() const {
	return true;
}

std::string ReversedReaction::getName() const {
	return getOrig()->getName();
}

const char* ReversedReaction::getCName() const {
	return getOrig()->getCName();
}

bool ReversedReaction::isExchange() const {
	return getOrig()->isExchange();
}

void ReversedReaction::setExchange(bool exchange) {
	getOrig()->setExchange(exchange);
}

const boost::shared_ptr<Model> ReversedReaction::getOwner() const {
	return getOrig()->getOwner();
}

std::string ReversedReaction::toString() const {
	return "reversed " + getOrig()->toString();
}

} /* namespace metaopt */
