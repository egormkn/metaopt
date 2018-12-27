#include "TextLoader.h"

#include "model/Model.h"
#include "model/impl/FullModel.h"

#include <string>
#include <fstream>
#include <iostream>
#include <tuple>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

namespace metaopt {

    void TextLoader::load(const std::vector<std::tuple<int, int, double>>& stoichiometry,
                          const std::vector<std::tuple<double, double, double>>& reactions, int species) {
        _model = boost::shared_ptr<FullModel>(new FullModel());

        _metabolites.resize(species);
        _reactions.resize(reactions.size());

        int num_met = species;
        for (int i = 0; i < num_met; i++) {
            const string& name = "Met" + std::to_string(i + 1);
            MetabolitePtr met = _model->createMetabolite(name);
            met->setBoundaryCondition(false);
            _metabolites[i] = met;
            _metabolites_map[name] = met;
        }

        int num_rxn = reactions.size();
        for (int i = 0; i < num_rxn; i++) {
            const string& name = "React" + std::to_string(i + 1);
            ReactionPtr rxn = _model->createReaction(name);

            double lb = std::get<0>(reactions[i]);
            double ub = std::get<1>(reactions[i]);
            double obj = std::get<2>(reactions[i]);
            rxn->setLb(lb);
            rxn->setUb(ub);
            rxn->setObj(obj);
            _reactions[i] = rxn;
        }

        std::vector<int> total(reactions.size(), 0);

        for (int i = 0; i < stoichiometry.size(); i++) {
            int row = std::get<0>(stoichiometry[i]) - 1;
            int col = std::get<1>(stoichiometry[i]) - 1;
            double value = std::get<2>(stoichiometry[i]);
            MetabolitePtr met = _metabolites[row];
            ReactionPtr& rxn = _reactions[col];
            rxn->setStoichiometry(met, value);
            total[col] += 1;
        }

        for (int i = 0; i < total.size(); i++) {
            ReactionPtr& rxn = _reactions[i];
            rxn->setExchange(total[i] <= 1);
        }
    }

    ModelPtr TextLoader::getModel() const {
        return _model;
    }

    void TextLoader::clear() {
        // clear list of metabolites
        _metabolites.clear();

        // clear loaded model
        _model.reset();

        _reactions.clear();

        _metabolites_map.clear();
    }

    ReactionPtr TextLoader::getReaction(int i) const {
        return _reactions[i];
    }

}
