# Makefile for compiling all unit tests

# General compiler specifications
CXX := g++
CXXFLAGS := -O0 -g $(INCLUDE_FLAGS)
CPPFLAGS :=
INCLUDE_FLAGS := -isystem /usr/local/include
LIBRARY_FLAGS := -lm -lgtest

# Preliminary definitions used in primary instance
export GLOBAL_LOCATION = $(CURDIR)
SUBDIRS := $(filter %/,$(wildcard src/test_*/))  # filter gets around glob() bug 10278 in glibc

# Preliminary definitions used in local instances
BIN_DIR := bin/
EXECUTABLE_NAME := run_tests
EXECUTABLE := $(BIN_DIR)$(EXECUTABLE_NAME)_$(words $(wildcard $(BIN_DIR)$(EXECUTABLE_NAME)_*))
TEST_SRC := $(wildcard src/*.cpp) $(wildcard $(LOCAL_LOCATION)/*.cpp)
ATHENA_SRC := $(wildcard ../src/*.cpp) \
	      $(wildcard ../src/bvals/*.cpp) \
	      $(wildcard ../src/integrators/*.cpp) \
	      $(wildcard ../src/outputs/*.cpp) \
	      ../src/coordinates/$(COORDINATES_FILE) \
	      ../src/convert_var/$(CONVERT_VAR_FILE) \
	      ../src/pgen/$(PROBLEM_FILE) \
	      ../src/reconstruct/$(RECONSTRUCT_FILE) \
	      ../src/rsolvers/$(RSOLVER_FILE)
SRC_DIR := $(dir $(TEST_SRC)) $(dir $(ATHENA_SRC))
OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(TEST_SRC)) $(filter-out main.cpp,$(notdir $(ATHENA_SRC))))
OBJ_FILES := $(OBJ_FILES:.cpp=.o)
VPATH := $(SRC_DIR)

# Generally useful targets
.PHONY : test dirs local clean
test : dirs
	$(foreach SUBDIR,$(SUBDIRS),$(MAKE) -C $(SUBDIR) all;)
dirs :
	@mkdir -p $(OBJ_DIR)
	@mkdir -p data
	@mkdir -p logs
	@mkdir -p plots
local : $(EXECUTABLE)

# Link objects into executable
$(EXECUTABLE) : $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(LIBRARY_FLAGS) -o $@ $^

# Create objects from source files
$(OBJ_DIR)%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@ 

# Cleanup
clean :
	-rm -rf $(OBJ_DIR)*
	-rm -rf $(BIN_DIR)$(EXECUTABLE_NAME)*
