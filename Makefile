# build the whole thing

ifeq ($(MATLAB),off)
	MATLAB_DEP=
	TEST_MATLAB=
else
	MATLAB_DEP=matlab
	TEST_MATLAB=test_matlab
endif

ifeq ($(SBML),off)
	SBML_DEP=
	TEST_SBML=
else
	SBML_DEP=sbml
	TEST_SBML=test_sbml
endif

.PHONY: all
all: thermo $(MATLAB_DEP)

.PHONY: thermo
thermo: scip $(SBML_DEP)
	cd thermo; make SBML=$(SBML) MATLAB=$(MATLAB)

.PHONY: scip
scip:
	cd third_party; make scip

.PHONY: sbml
sbml:
	cd third_party; make libsbml

.PHONY: matlab
matlab: thermo
	cd matlab; make

.PHONY: test
test: $(TEST_MATLAB)

.PHONY: test_matlab
test_matlab:
	./matlab/bin/metaopt examples/ecoli.mat results.mat tfva 1

.PHONY: clean
clean: 
	cd thermo; make clean
	cd matlab; make clean
