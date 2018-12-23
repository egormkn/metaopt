cmake_minimum_required(VERSION 3.12)

project(thermo)

option(MATLAB "Enable support for loading .mat files" OFF)
option(SBML "Enable support for loading SBML models" ON)

set(SRC_METAOPT_MODEL_SCIP
        src/model/scip/AbstractScipFluxModel.cpp
        src/model/scip/DualPotentials.cpp
        src/model/scip/ISSupply.cpp
        src/model/scip/LPFlux.cpp
        src/model/scip/LPPotentials.cpp
        src/model/scip/ModelAddOn.cpp
        src/model/scip/PotSpaceConstraint.cpp
        src/model/scip/ReducedScipFluxModel.cpp
        src/model/scip/ScipModel.cpp
        src/model/scip/Solution.cpp)

set(SRC_METAOPT_MODEL_SCIP_ADDON
        src/model/scip/addon/PotentialDifferences.cpp
        src/model/scip/addon/ReactionDirections.cpp)

set(SRC_METAOPT_SCIP_CONSTRAINTS
        src/scip/constraints/PotBoundPropagation2.cpp
        src/scip/constraints/RelaxedNaiveThermoConstraint.cpp
        src/scip/constraints/SteadyStateConstraint.cpp
        src/scip/constraints/ThermoConstraintHandler.cpp
        src/scip/constraints/ThermoInfeasibleSetPool.cpp)

set(SRC_METAOPT_SCIP_HEUR
        src/scip/heur/CycleDeletionHeur.cpp)

set(SRC_METAOPT_ALGORITHMS
        src/algorithms/BlockingSet.cpp
        src/algorithms/FCA.cpp
        src/algorithms/FluxForcing.cpp
        src/algorithms/FVA.cpp
        src/algorithms/ModelFactory.cpp)

set(SRC_METAOPT
        src/Uncopyable.cpp)

set(SRC_METAOPT_MODEL
        src/model/Coupling.cpp
        src/model/Metabolite.cpp
        src/model/Model.cpp
        src/model/Precision.cpp
        src/model/Reaction.cpp)

set(SRC_METAOPT_MODEL_IMPL
        src/model/impl/FullModel.cpp)

if(SBML)
    set(SRC_METAOPT_MODEL_SBML
            src/model/sbml/SBMLLoader.cpp)
endif()

if(MATLAB)
    set(SRC_METAOPT_MODEL_MATLAB
            src/model/matlab/MatlabLoader.cpp)
endif()


set(SOURCES
        ${SRC_METAOPT_MODEL_SCIP}
        ${SRC_METAOPT_MODEL_SCIP_ADDON}
        ${SRC_METAOPT_SCIP_CONSTRAINTS}
        ${SRC_METAOPT_SCIP_HEUR}
        ${SRC_METAOPT_ALGORITHMS}
        ${SRC_METAOPT}
        ${SRC_METAOPT_MODEL}
        ${SRC_METAOPT_MODEL_IMPL}
        ${SRC_METAOPT_MODEL_SBML}
        ${SRC_METAOPT_MODEL_MATLAB})

add_library(thermo STATIC ${SOURCES})

target_include_directories(thermo PUBLIC src)

add_definitions(-DLIBSBML_USE_CPP_NAMESPACE=1)
target_include_directories(thermo PUBLIC ${LIBSBML_ROOT_SOURCE_DIR}/src)

target_link_libraries(thermo sbml libscip)