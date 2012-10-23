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
 * SBMLLoader.cpp
 *
 *  Created on: 03.04.2012
 *      Author: arnem
 */

#include "SBMLLoader.h"
#include "metaopt/model/Model.h"
#include <string>
#include <boost/unordered_map.hpp>
#include <metaopt/model/impl/FullModel.h>

using namespace std;
using namespace boost;

namespace metaopt {

SBMLLoader::SBMLLoader() {
	// TODO Auto-generated constructor stub

}

SBMLLoader::~SBMLLoader() {
	// TODO Auto-generated destructor stub
}

void SBMLLoader::buildReaction(const libsbml::Reaction* r) {
	string name = r->getId();
	ReactionPtr rxn = _model->createReaction(name);
	int num_reactants = r->getNumReactants();
	for(int i = 0; i < num_reactants; i++) {
		const libsbml::SpeciesReference* m = r->getReactant(i);
		MetabolitePtr met = _metabolites[m->getSpecies()];
		rxn->setReactant(met, m->getStoichiometry());
	}
	int num_products = r->getNumProducts();
	for(int i = 0; i < num_products; i++) {
		const libsbml::SpeciesReference* m = r->getProduct(i);
		MetabolitePtr met = _metabolites[m->getSpecies()];
		rxn->setProduct(met, m->getStoichiometry());
	}
	// every reaction that only involves one metabolite is an exchange reaction
	rxn->setExchange(num_products + num_reactants <= 1);

	double lb = r->getKineticLaw()->getParameter("LOWER_BOUND")->getValue();
	double ub = r->getKineticLaw()->getParameter("UPPER_BOUND")->getValue();
	double obj = r->getKineticLaw()->getParameter("OBJECTIVE_COEFFICIENT")->getValue();
	rxn->setLb(lb);
	rxn->setUb(ub);
	rxn->setObj(obj);
}

void SBMLLoader::load(const libsbml::Model* m) {
	_model = shared_ptr<FullModel>(new FullModel());
	int num_met = m->getNumSpecies();
	for(int i = 0; i < num_met; i++) {
		const libsbml::Species* s = m->getSpecies(i);
		string name = s->getId();
		MetabolitePtr met = _model->createMetabolite(name);
		met->setBoundaryCondition(s->getBoundaryCondition());
		_metabolites[name] = met;
	}

	int num_rxn = m->getNumReactions();
	for(int i = 0; i < num_rxn; i++) {
		const libsbml::Reaction* r = m->getReaction(i);
		buildReaction(r);
	}
}

ModelPtr SBMLLoader::getModel() const {
  return _model;
}

void SBMLLoader::clear() {
	// clear list of metabolites
	_metabolites.clear();

	// clear loaded model
	_model.reset();
}

ModelPtr loadSBMLModel(const libsbml::Model* m) {
	SBMLLoader loader;
	loader.load(m);
	ModelPtr model = loader.getModel();
	return model;
}

}
