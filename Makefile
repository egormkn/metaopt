# Makefile from scratch so I understand the build process

SCIP_PATH=/home/bude/arnem/RAID/tfva/third_party/ziboptsuite-2.1.1/scip-2.1.1
SCIP_LIB=$(SCIP_PATH)/lib
SCIP_H=$(SCIP_PATH)/src

SBML_PATH=/home/bude/arnem/MolNet/impl2/third_party/libsbml-5.1.0-b0/build
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

SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_SBML_DIR)/%,$(SRC_METAOPT_MODEL_SBML))
SRC_METAOPT_MODEL+=$(patsubst %,$(SRC_METAOPT_MODEL_IMPL_DIR)/%,$(SRC_METAOPT_MODEL_IMPL))
SRC_METAOPT+=$(patsubst %,$(SRC_METAOPT_MODEL_DIR)/%,$(SRC_METAOPT_MODEL))
SRC+=$(patsubst %,$(SRC_METAOPT_DIR)/%,$(SRC_METAOPT))

SOURCES=$(patsubst %,$(SRC_DIR)/%,$(SRC))
OBJECTS=$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC))
LIBRARY=libthermo.so


CC=g++
LD=ld

CFLAGS=-Wall -fPIC
LDFLAGS=-shared

all : obj_dir $(BIN_DIR)/$(LIBRARY)

obj_dir :
	mkdir -p $(dir $(OBJECTS))

clean :
	rm -Rf $(OBJ_DIR)
	rm -Rf $(BIN_DIR)

$(BIN_DIR)/$(LIBRARY) : $(OBJECTS)
	mkdir -p $(BIN_DIR)
	$(LD) -o $@ $(LDFLAGS) $(OBJECTS) -L$(SCIP_LIB) -L$(SBML_LIB) -rpath=$(SCIP_LIB):$(SBML_LIB) -lsbml -lscip -lobjscip -llpispx -lnlpi -lsoplex -lzimpl -lz -lgmp -lreadline -lncurses -lm

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(SRC_DIR) -I$(SCIP_H) -I$(SBML_H)
