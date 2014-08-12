# Local makefile for compiling unit tests in this directory

# Files for conditional compilation
export COORDINATES_FILE := minkowski_cartesian.cpp
export CONVERT_VAR_FILE := adiabatic_hydro_gr.cpp
export PROBLEM_FILE := shock_tube_gr.cpp
export RECONSTRUCT_FILE := plm.cpp
export RSOLVER_FILE := hlle_gr.cpp

# Preliminary definitions
export LOCAL_LOCATION := $(CURDIR)

# Call global file
.PHONY : all
all :
	$(MAKE) -C $(GLOBAL_LOCATION) local
