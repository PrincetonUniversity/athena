# Makefile for running unit tests

# Files for conditional compilation
COORDINATES_FILE = cartesian.cpp
CONVERT_VAR_FILE = adiabatic_hydro_sr.cpp
PROBLEM_FILE = shock_tube_sr.cpp
RSOLVER_FILE = hllc_sr.cpp
RECONSTRUCT_FILE = plm.cpp

# General compiler specifications
CXX := g++
CXXFLAGS := -O0 -g $(INCLUDE_FLAGS)
CPPFLAGS :=
INCLUDE_FLAGS := -isystem /usr/local/include
LIBRARY_FLAGS := -lm -lgtest

# Preliminary definitions
EXECUTABLE := run_tests
TEST_SRC := $(wildcard src/*.cpp)
ATHENA_SRC := $(wildcard ../src/*.cpp) \
	      $(wildcard ../src/bvals/*.cpp) \
	      ../src/coordinates/$(COORDINATES_FILE) \
	      $(wildcard ../src/integrators/*.cpp) \
	      $(wildcard ../src/outputs/*.cpp) \
	      ../src/convert_var/$(CONVERT_VAR_FILE) \
	      ../src/pgen/$(PROBLEM_FILE) \
	      ../src/reconstruct/$(RECONSTRUCT_FILE) \
	      ../src/rsolvers/$(RSOLVER_FILE)
SRC_DIR := $(dir $(TEST_SRC)) $(dir $(ATHENA_SRC))
OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(filter-out main.cpp,$(notdir $(ATHENA_SRC))) $(notdir $(TEST_SRC)))
OBJ_FILES := $(OBJ_FILES:.cpp=.o)
VPATH := $(SRC_DIR)

# Generally useful targets
.PHONY : test dirs clean
test : dirs $(EXECUTABLE)
dirs :
	@mkdir -p $(OBJ_DIR)

# Link objects into executable
$(EXECUTABLE) : $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(LIBRARY_FLAGS) -o $@ $^

# Create objects from source files
$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@ 

# Cleanup
clean :
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)
