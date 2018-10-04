//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file show_config.cpp

// C++ headers
#include <iostream>
#include <sstream>

// Athena headers
#include "../athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ShowConfig(void)
//  \brief prints diagnostic messages about the configuration of an Athena++ executable

void ShowConfig(void) {
  std::cout<<"This Athena++ executable is configured with:" << std::endl;
  std::cout<<"  Problem generator:          " << PROBLEM_GENERATOR << std::endl;
  std::cout<<"  Coordinate system:          " << COORDINATE_SYSTEM << std::endl;
  if (NON_BAROTROPIC_EOS) {
    std::cout<<"  Equation of state:          adiabatic" << std::endl;
  } else {
    std::cout<<"  Equation of state:          isothermal" << std::endl;
  }
  std::cout<<"  Riemann solver:             " << RIEMANN_SOLVER << std::endl;

  if (SHEARING_BOX) {
    std::cout<<"  ShearingBox BCs:            ON" << std::endl;
  } else {
    std::cout<<"  ShearingBox BCs:            OFF" << std::endl;
  }

  if (SELF_GRAVITY_ENABLED == 1) {
    std::cout<<"  Self Gravity:               FFT" << std::endl;
  } else if (SELF_GRAVITY_ENABLED == 2) {
    std::cout<<"  Self Gravity:               Multigrid" << std::endl;
  } else {
    std::cout<<"  Self Gravity:               OFF" << std::endl;
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    std::cout<<"  Magnetic fields:            ON" << std::endl;
  } else {
    std::cout<<"  Magnetic fields:            OFF" << std::endl;
  }
  if (RELATIVISTIC_DYNAMICS) {
    std::cout<<"  Relativistic dynamics:      ON " << std::endl;
  } else {
    std::cout<<"  Relativistic dynamics:      OFF " << std::endl;
  }
  if (GENERAL_RELATIVITY) {
    std::cout<<"  General Relativity:         ON " << std::endl;
  } else {
    std::cout<<"  General Relativity:         OFF " << std::endl;
  }
  if (SINGLE_PRECISION_ENABLED) {
    std::cout<<"  Floating point precision:   single" << std::endl;
  } else {
    std::cout<<"  Floating point precision:   double" << std::endl;
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
  std::cout<<"  HDF5 Output:                ON" << std::endl;
  if (H5_DOUBLE_PRECISION_ENABLED) {
    std::cout<<"  HDF5 precision:             double" << std::endl;
  } else {
    std::cout<<"  HDF5 precision:             single" << std::endl;
  }
#else
  std::cout<<"  HDF5 Output:                OFF" << std::endl;
#endif

  return;
}
