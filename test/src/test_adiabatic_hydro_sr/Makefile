# Local makefile for compiling unit tests in this directory

# Files for conditional compilation
export COORDINATES_FILE := cartesian.cpp
export CONVERT_VAR_FILE := adiabatic_hydro_sr.cpp
export PROBLEM_FILE := shock_tube_sr.cpp
export RECONSTRUCT_FILE := plm.cpp
export RSOLVER_FILE := hllc_sr.cpp

# Preliminary definitions
export LOCAL_LOCATION := $(CURDIR)

# Call global file
.PHONY : all
all :
	@cp defs.hpp.in ../../../src/defs.hpp
	$(MAKE) -C $(GLOBAL_LOCATION) local
