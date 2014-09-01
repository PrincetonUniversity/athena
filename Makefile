# Rough Draft of Athena++ Makefile

# Files for conditional compilation

COORDINATES_FILE = cartesian.cpp
CONVERT_VAR_FILE = isothermal_hydro.cpp
PROBLEM_FILE = shock_tube.cpp
RSOLVER_FILE = hlle.cpp
RECONSTRUCT_FILE = plm.cpp

# General compiler specifications

CXX := g++
CXXFLAGS := -O3
LIBRARY_FLAGS := -lm

# Preliminary definitions

EXE_DIR := bin/
EXECUTABLE := $(EXE_DIR)athena
SRC_FILES := $(wildcard src/*.cpp) \
	     src/coordinates/$(COORDINATES_FILE) \
	     $(wildcard src/fluid/*.cpp) \
	     $(wildcard src/fluid/bvals/*.cpp) \
	     $(wildcard src/fluid/integrators/*.cpp) \
	     $(wildcard src/outputs/*.cpp) \
	     src/fluid/eos/$(CONVERT_VAR_FILE) \
	     src/pgen/$(PROBLEM_FILE) \
	     src/fluid/integrators/reconstruct/$(RECONSTRUCT_FILE) \
	     src/fluid/integrators/rsolvers/$(RSOLVER_FILE)
SRC_DIR := $(dir $(SRC_FILES))
OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))
VPATH := $(SRC_DIR)

# Generally useful targets

.PHONY : all dirs clean

all : dirs $(EXECUTABLE)

dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Link objects into executable

$(EXECUTABLE) : $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(LIBRARY_FLAGS) -o $@ $^

# Create objects from source files

$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

# Cleanup

clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)
