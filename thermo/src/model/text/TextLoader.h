#ifndef TEXTLOADER_H
#define TEXTLOADER_H

#include "Properties.h"
#include "Uncopyable.h"
#include "model/Model.h"
#include "model/Metabolite.h"
#include "model/impl/FullModel.h"

#include <boost/unordered_map.hpp>

namespace metaopt {

/**
 * The TextLoader can be used to load networks in form of plain text files.
 */
class TextLoader : Uncopyable {
public:
	void load(const std::vector<std::tuple<int, int, double>>& stoichiometry,
	        const std::vector<std::tuple<double, double, double>>& reactions, int species);

	/** fetch the loaded model */
	ModelPtr getModel() const;

	ReactionPtr getReaction(int i) const;

	/** clear loader for loading another model */
	void clear();

private:
	boost::shared_ptr<FullModel> _model;
	/** map of metabolite indices to created objects */
	std::vector<MetabolitePtr> _metabolites;

    boost::unordered_map<std::string, MetabolitePtr> _metabolites_map;

	/** map of reaction indices to created objects */
	std::vector<ReactionPtr> _reactions;

};

}

#endif // TEXTLOADER_H
