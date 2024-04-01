//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file show_config.cpp

// C headers

// C++ headers
#include <iostream>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ShowConfig()
//! \brief prints diagnostic messages about the configuration of an Athena++ executable

void ShowConfig() {
  // To match configure.py output: use 2 space indent for option, value output starts on
  // column 30
  std::cout<<"This Athena++ executable is configured with:" << std::endl;
  std::cout<<"  Problem generator:          " << PROBLEM_GENERATOR << std::endl;
  std::cout<<"  Coordinate system:          " << COORDINATE_SYSTEM << std::endl;
  std::cout<<"  Equation of state:          " << EQUATION_OF_STATE << std::endl;
  std::cout<<"  Riemann solver:             " << RIEMANN_SOLVER << std::endl;
  if (MAGNETIC_FIELDS_ENABLED) {
    std::cout<<"  Magnetic fields:            ON" << std::endl;
  } else {
    std::cout<<"  Magnetic fields:            OFF" << std::endl;
  }
  if (RELATIVISTIC_DYNAMICS) { // configure.py output: "Special relativity"
    std::cout<<"  Relativistic dynamics:      ON " << std::endl;
  } else {
    std::cout<<"  Relativistic dynamics:      OFF " << std::endl;
  }
  if (GENERAL_RELATIVITY) {
    std::cout<<"  General relativity:         ON " << std::endl;
  } else {
    std::cout<<"  General relativity:         OFF " << std::endl;
  }
  if (NR_RADIATION_ENABLED) {
    std::cout<<"  Radiative Transfer:         ON" << std::endl;
  } else {
    std::cout<<"  Radiative Transfer:         OFF" << std::endl;
  }
  if (IM_RADIATION_ENABLED) {
    std::cout<<"  Implicit Radiation:         ON" << std::endl;
  } else {
    std::cout<<"  Implicit Radiation:         OFF" << std::endl;
  }
  if (CR_ENABLED) {
    std::cout<<"  Cosmic Ray Transport:       ON" << std::endl;
  } else {
    std::cout<<"  Cosmic Ray Transport:       OFF" << std::endl;
  }
  if (CRDIFFUSION_ENABLED) {
    std::cout<<"  Cosmic Ray Diffusion:       ON" << std::endl;
  } else {
    std::cout<<"  Cosmic Ray Diffusion:       OFF" << std::endl;
  }

  // configure.py output: "Frame transformations"
  if (SELF_GRAVITY_ENABLED == 1) {
    std::cout<<"  Self-Gravity:               FFT" << std::endl;
  } else if (SELF_GRAVITY_ENABLED == 2) {
    std::cout<<"  Self-Gravity:               Multigrid" << std::endl;
  } else {
    std::cout<<"  Self-Gravity:               OFF" << std::endl;
  }
  if (STS_ENABLED) {
    std::cout<<"  Super-Time-Stepping:        ON" << std::endl;
  } else {
    std::cout<<"  Super-Time-Stepping:        OFF" << std::endl;
  }
  // configure.py output: +"Debug flags"
  // configure.py output: +"Code coverage flags"
  // configure.py output: +"Linker flags"
  if (SINGLE_PRECISION_ENABLED) {
    std::cout<<"  Floating-point precision:   single" << std::endl;
  } else {
    std::cout<<"  Floating-point precision:   double" << std::endl;
  }
  std::cout<<"  Number of ghost cells:      " << NGHOST << std::endl;
#ifdef MPI_PARALLEL
  std::cout<<"  MPI parallelism:            ON" << std::endl;
#else
  std::cout<<"  MPI parallelism:            OFF" << std::endl;
#endif
#ifdef OPENMP_PARALLEL
  std::cout<<"  OpenMP parallelism:         ON" << std::endl;
#else
  std::cout<<"  OpenMP parallelism:         OFF" << std::endl;
#endif

#ifdef FFT
  std::cout<<"  FFT:                        ON" << std::endl;
#else
  std::cout<<"  FFT:                        OFF" << std::endl;
#endif

#ifdef HDF5OUTPUT
  std::cout<<"  HDF5 output:                ON" << std::endl;
  if (H5_DOUBLE_PRECISION_ENABLED) {
    std::cout<<"  HDF5 precision:             double" << std::endl;
  } else {
    std::cout<<"  HDF5 precision:             single" << std::endl;
  }
#else
  std::cout<<"  HDF5 output:                OFF" << std::endl;
#endif

  std::cout<<"  Compiler:                   " << COMPILED_WITH << std::endl;
  std::cout<<"  Compilation command:        " << COMPILER_COMMAND
           << COMPILED_WITH_OPTIONS << std::endl;
  // configure.py output: Doesnt append "Linker flags" in prev. output (excessive space!)
  return;
}
