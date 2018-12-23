# Makefile from scratch so I understand the build process
#
#

# default path for scip
SCIP_PATH=../third_party/scip
SCIP_LIB=$(SCIP_PATH)/lib
SCIP_H=$(SCIP_PATH)/src
SOPLEX_H=../third_party/soplex/src

#default path for matlab
MATLAB_PATH=../third_party/matlab
MATLAB_LIB=$(MATLAB_PATH)/bin
MATLAB_H=$(MATLAB_PATH)/include

#default path to libthermo
THERMO_PATH=../thermo
THERMO_LIB=$(THERMO_PATH)/bin
THERMO_H=$(THERMO_PATH)/src

OBJ_DIR=obj
BIN_DIR=bin

SRC_DIR=src
SRC=metaopt_mex_octave.cpp

SOURCES=$(patsubst %,$(SRC_DIR)/%,$(SRC))
OBJECTS=$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC))

MEXTGT=metaopt.mex

EXE=metaopt

libMATLAB=-loctave -lmatio


# don't use gcc for matlab linking, but use mex
# mex also wants to do everything in a single step, so no separation between compiling and linking
 MEX=mkoctfile --mex
# remark actually we could also use gcc directly
# we will have to set LD_PRELOAD to use the correct c++-libs anyway

MEXFLAG= -g --output $(BIN_DIR)/$(MEXTGT) -Wl,-rpath=$(SCIP_LIB):$(THERMO_LIB)

# use ordinary g++
CC=g++

#CFLAG=-Wall -g3 -fno-inline -O0 -D_GLIBCXX_DEBUG
ifeq ($(OPT),opt)
DEBUGFLAGS=-O3
else
DEBUGFLAGS=-g3 -fno-inline -O0
endif

CFLAG=-Wall $(DEBUGFLAGS)
RPATH=-Xlinker -rpath=$(SCIP_LIB):$(THERMO_LIB):$(MATLAB_LIB)

all : $(BIN_DIR)/$(MEXTGT)
# all : $(BIN_DIR)/$(EXE)

clean :
	rm -Rf $(BIN_DIR)

$(BIN_DIR)/$(EXE) : $(SOURCES) $(THERMO_LIB)/libthermo.so
	mkdir -p $(BIN_DIR)
	$(CC) -o $(BIN_DIR)/$(EXE) $(CFLAG) -I$(THERMO_H) -I$(SCIP_H) -I$(SOPLEX_H) -I$(MATLAB_H) -I. $(SOURCES) -L$(THERMO_LIB) -L$(SCIP_LIB) -L$(MATLAB_LIB) $(RPATH) $(libMATLAB) -lthermo -lscip -lobjscip -llpispx -lsoplex -lzimpl -lz -lgmpxx -lreadline -lncurses -lm

$(BIN_DIR)/$(MEXTGT) : $(SOURCES) $(THERMO_LIB)/libthermo.so
	mkdir -p $(BIN_DIR)
	$(MEX) $(MEXFLAG) -I$(THERMO_H) -I$(SCIP_H) -I$(SOPLEX_H) -I$(MATLAB_H) -I. $(SOURCES) -L$(THERMO_LIB) -L$(SCIP_LIB) -L$(MATLAB_LIB) $(libMATLAB) -lthermo -lscip -lobjscip -llpispx -lsoplex -lzimpl -lz -lgmp -lreadline -lncurses -lm