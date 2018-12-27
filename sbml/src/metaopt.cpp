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

#include <iostream>
#include <cinttypes>

#include <boost/unordered_map.hpp>
#include <boost/program_options.hpp>
#include <sbml/Model.h>

#include "Properties.h"
#include "algorithms/FVA.h"
#include "algorithms/BlockingSet.h"

#include "model/Precision.h"
#include "model/sbml/SBMLLoader.h"
#include "model/scip/ScipModel.h"
#include "model/scip/LPFlux.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/constraints/SteadyStateConstraint.h"
#include "scip/constraints/ThermoConstraintHandler.h"
#include "scip/heur/CycleDeletionHeur.h"

namespace metaopt {

    using namespace boost;
    using namespace std;
    using namespace libsbml;

    int fba(const libsbml::Model* m) {
        SBMLLoader loader;
        loader.load(m);
        ModelPtr model = loader.getModel();

        LPFluxPtr flux(new LPFlux(model, true));
        flux->solve();

        cout << "Objective value: " << flux->getObjVal() << endl;

        for (int i = 0; i < model->getReactions().size(); i++) {
            ReactionPtr rxn = loader.getReaction(i);
            cout << flux->getFlux(rxn) << endl;
        }

        return 0;
    }



/*
    void createExtraConstraint(ScipModelPtr scip, OctaveLoader &loader, mxArray *A, mxArray *a) {
        assert(A != NULL); // A does not exist
        assert(a != NULL); // a does not exist
        assert(!mxIsChar(A) && !mxIsComplex(A) && mxIsNumeric(A) &&
               mxIsSparse(A)); // A must be sparse and contain reals
        assert(!mxIsChar(a) && !mxIsComplex(a) && mxIsNumeric(a) && !mxIsSparse(a)); // a must be full and contain reals

        double *data = mxGetPr(A);
        mwIndex *rows = mxGetIr(A);
        mwIndex *pos = mxGetJc(A);

        int numRows = mxGetM(A);
        int numCols = mxGetN(A);

        assert((int) mxGetNumberOfElements(a) == numRows); // a does not have as many elements as A has rows

        vector<vector<SCIP_VAR *> > vars(numRows);
        vector<vector<double> > coef(numRows);

        ModelPtr model = scip->getModel();
        PrecisionPtr prec = model->getCoefPrecision();

        for (int i = 0; i < numCols; i++) {
            ReactionPtr rxn = loader.getReaction(i);
            SCIP_VAR *var = scip->getFlux(rxn);
            for (mwIndex j = pos[i]; j < pos[i + 1]; j++) {
                if (data[j] < -prec->getCheckTol() || data[j] > prec->getCheckTol()) {
                    vars.at(rows[j]).push_back(var);
                    coef.at(rows[j]).push_back(data[j]);
                    rxn->setProblematic(true);
                }
            }
        }

        double *aData = mxGetPr(a);
        double lhs = -INFINITY;

        for (int i = 0; i < numRows; i++) {
            vector<SCIP_VAR *> &cvars = vars[i];
            vector<double> &ccoef = coef[i];
            SCIP_CONS *cons = NULL;
            BOOST_SCIP_CALL(SCIPcreateConsLinear(scip->getScip(), &cons, "additional cons", cvars.size(), cvars.data(),
                                                 ccoef.data(), lhs, aData[i], true, true, true, true, true, false,
                                                 false, false, false, false));
            BOOST_SCIP_CALL(SCIPaddCons(scip->getScip(), cons));
            BOOST_SCIP_CALL(SCIPreleaseCons(scip->getScip(), &cons));
        }
    }
*/

    /**
     * Octave wrapper to thermodynamically constrained FBA
     */
    int tfba(const libsbml::Model* m, int nargout) {
        SBMLLoader loader;
        loader.load(m);
        ModelPtr model = loader.getModel();
        ScipModelPtr scip(new ScipModel(model));
        createSteadyStateConstraint(scip);
        void *extra_A = nullptr;
        void *extra_a = nullptr;
        if (extra_A != NULL) {
            assert(extra_a != NULL);
            // TODO
            // createExtraConstraint(scip, loader, extra_A, extra_a);
        }
        createThermoConstraint(scip);
        createCycleDeletionHeur(scip);
        //registerExitEventHandler(scip);

        scip->solve();

        bool solFound = scip->isOptimal();

        if (nargout >= 1) {
            if (solFound) {
                cout << "Objective value: " << scip->getObjectiveValue() << endl;
            } else if (scip->isUnbounded()) {
                if (scip->isMaximize()) {
                    cout << "Objective value: " << INFINITY << endl;
                } else {
                    cout << "Objective value: " << -INFINITY << endl;
                }
            } else {
                if (!scip->isInfeasible()) {
                    cout << "Warning: Did not compute a solution, but not infeasible nor unbounded" << endl;
                }
                cout << "Objective value: " << NAN << endl;
            }
        }
        if (nargout >= 2) {
            if (solFound) {
                for (int i = 0; i < model->getReactions().size(); i++) {
                    ReactionPtr rxn = loader.getReaction(i);
                    cout << scip->getCurrentFlux(rxn) << endl;
                }
            } else {
                cout << "No flux solution found" << endl;
            }
        }
        if (nargout >= 3 && (solFound || scip->hasPrimalRay())) {
            LPFluxPtr flux(new LPFlux(model, true));
            //LPFluxPtr m(new LPFlux(model, true));
            flux->set(scip); // copy optimal solution
            boost::shared_ptr<std::vector<DirectedReaction> > block = findBlockingSet(flux);
            int num_fwd = 0;
            int num_bwd = 0;
                    foreach(DirectedReaction &d, *block) {
                            if (d._fwd) num_fwd++;
                            else num_bwd++;
                        }

            cout << "num fwd reactions in blocking set: " << num_fwd << endl;
            cout << "num bwd reactions in blocking set: " << num_bwd << endl;

            // TODO
            /*mxArray *blockFwd = mxCreateNumericMatrix(num_fwd, 1, mxUINT32_CLASS, mxREAL);
            mxArray *blockBwd = mxCreateNumericMatrix(num_bwd, 1, mxUINT32_CLASS, mxREAL);

            // it is important to get the same type, else the indices of the array don't match
            uint32_t *fwd = (uint32_t *) mxGetData(blockFwd);
            uint32_t *bwd = (uint32_t *) mxGetData(blockBwd);

            unordered_map<ReactionPtr, int> rxnIndex;
            int numRxns = model->getReactions().size();
            for (int i = 0; i < numRxns; i++) {
                rxnIndex[loader.getReaction(i)] = i;
            }

            int i_fwd = 0;
            int i_bwd = 0;
                    foreach(DirectedReaction &d, *block) {
                            if (d._fwd) {
                                assert(rxnIndex.find(d._rxn) != rxnIndex.end());
                                fwd[i_fwd] = rxnIndex[d._rxn] + 1; // octave indexing starts at 1
                                cout << "blocked fwd reaction: " << d._rxn->getName() << " (" << rxnIndex[d._rxn] << ")"
                                     << endl;
                                i_fwd++;
                            } else {
                                assert(rxnIndex.find(d._rxn) != rxnIndex.end());
                                bwd[i_bwd] = rxnIndex[d._rxn] + 1; // octave indexing starts at 1
                                cout << "blocked bwd reaction: " << d._rxn->getName() << " (" << rxnIndex[d._rxn] << ")"
                                     << endl;
                                i_bwd++;
                            }
                        }
            assert(i_fwd == num_fwd);
            assert(i_bwd == num_bwd);

            const char *fwd_name = "fwd";
            const char *bwd_name = "bwd";
            const char *fieldnames[2] = {fwd_name, bwd_name};

            mxArray *block_m = mxCreateStructMatrix(1, 1, 2, fieldnames);
            mxSetField(block_m, 0, fwd_name, blockFwd);
            mxSetField(block_m, 0, bwd_name, blockBwd);

            mxSetCell(argout, 2, block_m);*/
        }

        return 0;
    }

/**
 * Stores fva results in the octave structure (as a matrix with two columns)
 * The order of the elements is the same as they were loaded by the octave loader
 *
 * @param m the result matrix that is created
 * @param model the model on which FVA had been run
 * @param loader the loader with which the model had been loaded
 * @param min the minimal flux values
 * @param max the maximal flux values
*/

    void convert_fva_result(metaopt::ModelPtr &model, metaopt::SBMLLoader &loader,
                            boost::unordered_map<metaopt::ReactionPtr, double> &min,
                            boost::unordered_map<metaopt::ReactionPtr, double> &max) {
        for (int i = 0; i < model->getReactions().size(); i++) {
            ReactionPtr rxn = loader.getReaction(i);
            unordered_map<ReactionPtr, double>::iterator iter_min = min.find(rxn);
            unordered_map<ReactionPtr, double>::iterator iter_max = max.find(rxn);

            double min_value, max_value;
            if (iter_min == min.end()) {
                cout << "Error: Did not run FVA for a reaction " << i << " - computation aborted?" << endl;
                min_value = -INFINITY;
            } else {
                min_value = iter_min->second;
            }

            if (iter_max == max.end()) {
                cout << "Error: Did not run FVA for a reaction " << i << " - computation aborted?" << endl;
                max_value = INFINITY;
            } else {
                max_value = iter_max->second;
            }
            cout << min_value << " " << max_value << endl;
        }
    }

/**
 * Octave wrapper to plain old FVA

    int fva(mxArray *argout, mxArray *argin) {
        int nargin = mxGetNumberOfElements(argin);

        if (nargin == 0) {
            cout << "fva: To run FVA you must specify a model as the second parameter!" << endl;
            return 18;
        }

        OctaveLoader loader;
        loader.load(mxGetCell(argin, 0));
        ModelPtr model = loader.getModel();

        unordered_map<ReactionPtr, double> min, max;

        fva(model, min, max);

        mxArray *res;

        convert_fva_result(&res, model, loader, min, max);
        mxSetCell(argout, 0, res);

        return 0;
    }

/*class ThermoModelFactory : public ModelFactory {
public:
	ScipModelPtr build(ModelPtr m) {
		ScipModelPtr scip(new ScipModel(m));
		createSteadyStateConstraint(scip);
		createThermoConstraint(scip);
		createCycleDeletionHeur(scip);
		//registerExitEventHandler(scip);

		return scip;
	}
};*/

/**
 * Octave wrapper to thermodynamically constrained FVA
 */
    int tfva(const libsbml::Model* m, int nargout, double timeout) {
        SBMLLoader loader;
        loader.load(m);
        ModelPtr model = loader.getModel();

        // FVA settings
        FVASettingsPtr settings(new FVASettings());

        settings->reactions = model->getReactions();
        if (timeout > 0) {
            settings->timeout = timeout;
        }

        unordered_map<metaopt::ReactionPtr, double> min, max;

        try {
            metaopt::tfva(model, settings, min, max);
        } catch (std::exception &ex) {
            std::cout << diagnostic_information(ex) << std::endl;
            return 19;
        }

        convert_fva_result(model, loader, min, max);

        return 0;
    }

    int help() {
        cout << "Metaopt Version " << VERSION << endl;
        cout << "Copyright (C) 2012 Arne Müller <arne.mueller@fu-berlin.de> " << endl;
        cout << endl;
        cout << "This program is free software: you can redistribute it and/or modify" << endl;
        cout << "it under the terms of the GNU General Public License as published by" << endl;
        cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
        cout << "(at your option) any later version." << endl;
        cout << endl;
        cout << "This program is distributed in the hope that it will be useful," << endl;
        cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
        cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
        cout << "GNU General Public License for more details." << endl;
        cout << endl;
        cout << "You should have received a copy of the GNU General Public License" << endl;
        cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
        cout << endl;
        cout << endl;
        cout << "Usage:" << endl;
        cout << "  <output parameters> = metaopt(<function name>, <input parameters>)" << endl;
        cout << "  <output parameters> = metaopt_time(<function name>, <timeout>, <input parameters>)" << endl;
        cout << "- function name has to be one of the functions listed below." << endl;
        cout
                << "- timout has has to be a number followed by either s, m, h, or d for seconds, minutes, hours, or days respectively."
                << endl;
        cout << "- input parameters depend on the function called." << endl;
        cout << "- output parameters depend on the function called." << endl;
        cout << endl;
        cout << "Available functions for metaopt:" << endl;
        cout << endl;
        cout << "fba: for ordinary fba (without thermodynamic constraints)." << endl;
        cout << "  Input Parameters:" << endl;
        cout << "   - a metabolic network model. See below for a precise specification." << endl;
        cout << "  Output Parameters:" << endl;
        cout << "   - optimal objective value (scalar)" << endl;
        cout << "   - an optimal flux distribution (column-vector)" << endl;
        cout << endl;
        cout << "tfba: for tfba (with thermodynamic constraints)." << endl;
        cout << "  Input Parameters:" << endl;
        cout << "   - a metabolic network model. See below for a precise specification." << endl;
        cout << "  Output Parameters:" << endl;
        cout << "   - optimal objective value (scalar)" << endl;
        cout << "   - an optimal flux distribution (column-vector)" << endl;
        cout << "   - blocking set. This is a struct with two fields 'fwd', 'bwd'." << endl;
        cout << "      - fwd is a list of reactions to be blocked in forward direction" << endl;
        cout << "      - bwd is a list of reactions to be blocked in backward direction" << endl;
        cout
                << "     If this blocking set is enforced, the tfba objective value will equal the fba objective value (only works without potential bounds)."
                << endl;
        cout << endl;
        cout << "fva: for fva (without thermodynamic constraints). This ignores any specified objective function."
             << endl;
        cout
                << "     If you want to run fva only on the optimal flux space, you have to make a separate call to fba and constrain the flux space accordingly."
                << endl;
        cout << "  Input Parameters:" << endl;
        cout << "   - a metabolic network model. See below for a precise specification." << endl;
        cout << "  Output Parameters:" << endl;
        cout << "   - variabilities of fluxes (matrix). " << endl;
        cout << "     The first column gives the minimal flux for each reaction." << endl;
        cout << "     The second column gives maximal flux for each reaction." << endl;
        cout << endl;
        cout << "tfva: for tfva (with thermodynamic constraints). This ignores any specified objective function."
             << endl;
        cout
                << "      If you want to run fva only on the optimal flux space, you have to make a separate call to tfba and constrain the flux space accordingly."
                << endl;
        cout << "  Input Parameters:" << endl;
        cout << "   - a metabolic network model. See below for a precise specification." << endl;
        cout << "   - a struct settings wich can have the following parameters:" << endl;
        cout << "      - reactions: Specifies a list (indices) of reactions for which fva should be performed." << endl;
        cout << "      - timeout: This feature is deprecated and does not work reliably. Use metaopt_timeout instead."
             << endl;
        cout << "  Output Parameters:" << endl;
        cout << "   - variabilities of fluxes (matrix). " << endl;
        cout << "     The first column gives the minimal flux for each reaction." << endl;
        cout << "     The second column gives maximal flux for each reaction." << endl;
        cout << endl;
        cout << "help: Prints this message." << endl << endl;
        cout << "Specification of metabolic network model:" << endl;
        cout << "  The metabolic network model is struct which is basically a COBRA model with additional fields:"
             << endl;
        cout << "   - S         stoichiometric matrix (must be a sparse matrix)" << endl;
        cout << "   - ub        flux upper bounds" << endl;
        cout << "   - lb        flux lower bounds" << endl;
        cout << "   - c         objective function on fluxes" << endl;
        cout << "   - rxns      cell array of reaction names" << endl;
        cout << "   - mets      cell array of metabolite names" << endl;
        cout
                << "   - boundary  list of indices specifying metabolites with boundary condition, i.e. steady-state assumption does not need to hold for this metabolite (experimental)."
                << endl;
        cout << "   - p_ub      potential upper bounds (experimental)" << endl;
        cout << "   - p_lb      potential lower bounds (experimental)" << endl;
        cout << "   - p_c       potential objective function (not yet used)" << endl;
        cout << "   - int_rxns  cell array of internal reaction names" << endl;
        cout << "   - precision precision with which the model should be solved" << endl;
        cout << "   - A         coefficient matrix (sparse) for extra constraints A v <= a (only for tfba)" << endl;
        cout << "   - a         left hand sides of extra constraints A v <= a (only for tfba)" << endl;
        cout << "  Remark: The right-hand side vector b (of S v = b) is assumed to be 0." << endl;
        cout << endl;
        cout << "Please report bugs, feature requests etc. to me: arne.mueller@fu-berlin.de" << endl;

        return 0;
    }

};

namespace opt = boost::program_options;

int main(int argc, const char *argv[]) {
    opt::command_line_parser parser(argc, argv);
    opt::variables_map args;

    try {
        opt::options_description options("Metaopt options");
        options.add_options()
                ("help,h", "Help screen")
                ("reactions,r", opt::value<string>()->required(), "Reactions file")
                ("metabolites,m", opt::value<string>()->required(), "Metabolites file")
                ("solver,s", opt::value<string>()->required(), "Solver type")
                ("output,o", opt::value<string>()->required()->default_value("output.txt"), "Output file");

        opt::parsed_options parsed_options = parser.options(options).run();
        opt::store(parsed_options, args);

        if (args.count("help")) {
            cout << options << endl;
            return 0;
        }

        opt::notify(args);

        string metabolites = args["metabolites"].as<string>(),
                reactions = args["reactions"].as<string>(),
                solver = args["solver"].as<string>(),
                output = args["output"].as<string>();

        cout << "Using metabolites file: " << metabolites << endl;
        cout << "Using reactions file: " << args["reactions"].as<string>() << endl;
        cout << "Using solver type: " << args["solver"].as<string>() << endl;
        cout << "Using output file: " << args["output"].as<string>() << endl;




        libsbml::SBMLDocument* document = libsbml::readSBML(reactions.c_str());
        if (document->getNumErrors() > 0) {
            cerr << "Encountered the following SBML errors:" << endl;
            document->printErrors(cerr);
            return 1;
        }
        unsigned int level   = document->getLevel  ();
        unsigned int version = document->getVersion();
        cout << endl
             << "File: " << reactions
             << " (Level " << level << ", version " << version << ")" << endl;
        libsbml::Model* model = document->getModel();
        if (model == 0)
        {
            cout << "No model present." << endl;
            return 1;
        }
        cout << "               "
             << (level == 1 ? "name: " : "  id: ")
             << (model->isSetId() ? model->getId() : std::string("(empty)")) << endl;
        if (model->isSetSBOTerm())
            cout << "      model sboTerm: " << model->getSBOTerm() << endl;
        cout << "functionDefinitions: " << model->getNumFunctionDefinitions() << endl;
        cout << "    unitDefinitions: " << model->getNumUnitDefinitions    () << endl;
        cout << "   compartmentTypes: " << model->getNumCompartmentTypes   () << endl;
        cout << "        specieTypes: " << model->getNumSpeciesTypes       () << endl;
        cout << "       compartments: " << model->getNumCompartments       () << endl;
        cout << "            species: " << model->getNumSpecies            () << endl;
        cout << "         parameters: " << model->getNumParameters         () << endl;
        cout << " initialAssignments: " << model->getNumInitialAssignments () << endl;
        cout << "              rules: " << model->getNumRules              () << endl;
        cout << "        constraints: " << model->getNumConstraints        () << endl;
        cout << "          reactions: " << model->getNumReactions          () << endl;
        cout << "             events: " << model->getNumEvents             () << endl;
        cout << endl;


        if (solver == "fba") {
            metaopt::fba(model);
        } else if (solver == "tfba") {
            metaopt::tfba(model, 2);
        } else if (solver == "fva") {
            cout << "Not implemented" << endl; // TODO
        } else if (solver == "tfva") {
            metaopt::tfva(model, 2, -1.0);
        }

        delete document;
    } catch (const std::exception &ex) {
        cerr << ex.what() << endl;
        return 1;
    }

    return 0;
}
