# Makefile from scratch so I understand the build process
#
# This make file can be controlled by the following parameters:
#	* SBML=on/off enables/disables support for loading SBML models (on default) 
#	* MATLAB=on/off enables/disables support for loading .mat files (on default)
#

# default path for scip
SCIP_PATH=$(CURDIR)/../third_party/scip
SCIP_LIB=$(SCIP_PATH)/lib
SCIP_H=$(SCIP_PATH)/src

#default path for sbml
SBML_PATH=$(CURDIR)/../third_party/libsbml/build
SBML_LIB=$(SBML_PATH)/lib
SBML_H=$(SBML_PATH)/include

#default path for matlab
MATLAB_PATH=$(CURDIR)/../third_party/matlab
MATLAB_LIB=$(MATLAB_PATH)/bin
MATLAB_H=$(MATLAB_PATH)/include


OBJ_DIR=obj
BIN_DIR=bin

SRC_DIR=src
SRC_METAOPT=Uncopyable.cpp
SRC_METAOPT_DIR=metaopt
SRC_METAOPT_MODEL=Model.cpp Metabolite.cpp Reaction.cpp Coupling.cpp
SRC_METAOPT_MODEL_DIR=model
SRC_METAOPT_MODEL_IMPL=FullModel.cpp
SRC_METAOPT_MODEL_IMPL_DIR=impl
ifeq ($(SBML),off)
	SRC_METAOPT_MODEL_SBML=
else
	SRC_METAOPT_MODEL_SBML=SBMLLoader.cpp
endif
SRC_METAOPT_MODEL_SBML_DIR=sbml
ifeq ($(MATLAB),off)
	SRC_METAOPT_MODEL_MATLAB=
else
	SRC_METAOPT_MODEL_MATLAB=MatlabLoader.cpp
endif
SRC_METAOPT_MODEL_MATLAB_DIR=matlab
SRC_METAOPT_MODEL_SCIP=ScipModel.cpp LPFlux.cpp ModelAddOn.cpp Solution.cpp LPPotentials.cpp DualPotentials.cpp ReducedScipFluxModel.cpp AbstractScipFluxModel.cpp ISSupply.cpp
SRC_METAOPT_MODEL_SCIP_DIR=scip
SRC_METAOPT_MODEL_SCIP_ADDON=PotentialDifferences.cpp ReactionDirections.cpp
SRC_METAOPT_MODEL_SCIP_ADDON_DIR=addon
SRC_METAOPT_SCIP=
SRC_METAOPT_SCIP_DIR=scip
SRC_METAOPT_SCIP_CONSTRAINTS=SteadyStateConstraint.cpp RelaxedNaiveThermoConstraint.cpp ThermoConstraintHandler.cpp PotBoundPropagation2.cpp
SRC_METAOPT_SCIP_CONSTRAINTS_DIR=constraints
SRC_METAOPT_SCIP_HEUR=CycleDeletionHeur.cpp
SRC_METAOPT_SCIP_HEUR_DIR=heur
SRC_METAOPT_ALGORITHMS=FCA.cpp FVA.cpp ModelFactory.cpp FluxForcing.cpp
SRC_METAOPT_ALGORITHMS_DIR=algorithms

SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_SBML_DIR)/%,$(SRC_METAOPT_MODEL_SBML))
SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_MATLAB_DIR)/%,$(SRC_METAOPT_MODEL_MATLAB))
SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_IMPL_DIR)/%,$(SRC_METAOPT_MODEL_IMPL))
SRC_METAOPT_MODEL_SCIP+=$(patsubst %,$(SRC_METAOPT_MODEL_SCIP_ADDON_DIR)/%,$(SRC_METAOPT_MODEL_SCIP_ADDON))
SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_SCIP_DIR)/%,$(SRC_METAOPT_MODEL_SCIP))
SRC_METAOPT+=$(patsubst %,$(SRC_METAOPT_MODEL_DIR)/%,$(SRC_METAOPT_MODEL))
SRC_METAOPT_SCIP+=$(patsubst %,$(SRC_METAOPT_SCIP_CONSTRAINTS_DIR)/%,$(SRC_METAOPT_SCIP_CONSTRAINTS))
SRC_METAOPT_SCIP+=$(patsubst %,$(SRC_METAOPT_SCIP_HEUR_DIR)/%,$(SRC_METAOPT_SCIP_HEUR))
SRC_METAOPT+=$(patsubst %,$(SRC_METAOPT_SCIP_DIR)/%,$(SRC_METAOPT_SCIP))
SRC_METAOPT+=$(patsubst %,$(SRC_METAOPT_ALGORITHMS_DIR)/%,$(SRC_METAOPT_ALGORITHMS))
SRC+=$(patsubst %,$(SRC_METAOPT_DIR)/%,$(SRC_METAOPT))

SOURCES=$(patsubst %,$(SRC_DIR)/%,$(SRC))
OBJECTS=$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC))
LIBRARY=libthermo.so


#depending if SBML is turned on, we have to setup the include and linker correctily
ifeq ($(SBML),off)
	ISBML=
	LSBML=
	RSBML=
	libSBML=
else
	ISBML=-I$(SBML_H)
	LSBML=-L$(SBML_LIB)
	RSBML=$(SBML_LIB)
	libSBML=-lsbml
endif
#depending if MATLAB is turned on, we have to setup the include and linker correctily
ifeq ($(MATLAB),off)
	IMATLAB=
	LMATLAB=
	RMATLAB=
	libMATLAB=
else
	IMATLAB=-I$(MATLAB_H)
	LMATLAB=-L$(MATLAB_LIB)
	RMATLAB=$(MATLAB_LIB)
	libMATLAB=-lmat -lmx -lut
endif


CC=g++
# use g++ as linker to get C++ stuff correctly.
# This has the additional consequence that the -rpath option has to be preceded by -Xlinker, since it is not a direct g++ command
LD=g++

DEBUGFLAGS=-g3 -fno-inline -O0 -D_GLIBCXX_DEBUG

CFLAGS=-Wall -fPIC $(DEBUGFLAGS)
LDFLAGS=-shared

all : obj_dir $(BIN_DIR)/$(LIBRARY)

obj_dir :
	mkdir -p $(dir $(OBJECTS))

clean :
	rm -Rf $(OBJ_DIR)
	rm -Rf $(BIN_DIR)

$(BIN_DIR)/$(LIBRARY) : $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS) -L$(SCIP_LIB) $(LSBML) $(LMATLAB) -Xlinker -rpath=$(SCIP_LIB):$(RSBML):$(RMATLAB) $(libSBML) $(libMATLAB) -lscip -lobjscip -llpispx -lnlpi -lsoplex -lzimpl -lz -lgmp -lreadline -lncurses -lm

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(SRC_DIR) -I$(SCIP_H) $(ISBML) $(IMATLAB)
