# Makefile from scratch so I understand the build process

SCIP_PATH=$(CURDIR)/../third_party/scip
SCIP_LIB=$(SCIP_PATH)/lib
SCIP_H=$(SCIP_PATH)/src

SBML_PATH=$(CURDIR)/../third_party/libsbml/build
SBML_LIB=$(SBML_PATH)/lib
SBML_H=$(SBML_PATH)/include

OBJ_DIR=obj
BIN_DIR=bin

SRC_DIR=src
SRC_METAOPT=
SRC_METAOPT_DIR=metaopt
SRC_METAOPT_MODEL=Model.cpp Metabolite.cpp Reaction.cpp
SRC_METAOPT_MODEL_DIR=model
SRC_METAOPT_MODEL_IMPL=FullModel.cpp
SRC_METAOPT_MODEL_IMPL_DIR=impl
SRC_METAOPT_MODEL_SBML=SBMLLoader.cpp
SRC_METAOPT_MODEL_SBML_DIR=sbml
SRC_METAOPT_MODEL_SCIP=ScipModel.cpp LPFlux.cpp
SRC_METAOPT_MODEL_SCIP_DIR=scip

SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_SBML_DIR)/%,$(SRC_METAOPT_MODEL_SBML))
SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_IMPL_DIR)/%,$(SRC_METAOPT_MODEL_IMPL))
SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_SCIP_DIR)/%,$(SRC_METAOPT_MODEL_SCIP))
SRC_METAOPT+=$(patsubst %,$(SRC_METAOPT_MODEL_DIR)/%,$(SRC_METAOPT_MODEL))
SRC+=$(patsubst %,$(SRC_METAOPT_DIR)/%,$(SRC_METAOPT))

SOURCES=$(patsubst %,$(SRC_DIR)/%,$(SRC))
OBJECTS=$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC))
LIBRARY=libthermo.so


CC=g++
# use g++ as linker to get C++ stuff correctly.
# This has the additional consequence that the -rpath option has to be preceded by -Xlinker, since it is not a direct g++ command
LD=g++

DEBUGFLAGS=-g3 -fno-inline -O0

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
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS) -L$(SCIP_LIB) -L$(SBML_LIB) -Xlinker -rpath=$(SCIP_LIB):$(SBML_LIB) -lsbml -lscip -lobjscip -llpispx -lnlpi -lsoplex -lzimpl -lz -lgmp -lreadline -lncurses -lm

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(SRC_DIR) -I$(SCIP_H) -I$(SBML_H)
