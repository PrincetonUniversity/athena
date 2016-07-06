//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Athena headers
#include "../athena.hpp"

// C++ headers
#include <iostream>
#include <sstream>

//======================================================================================
//! \file show_config.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn void ShowConfig(void)
//  \brief prints diagnostic messages about the configuration of an Athena++ executable

void ShowConfig(void)
{
  std::cout<<"This Athena++ executable is configured with:" << std::endl;
  std::cout<<"  Problem generator:          " << PROBLEM_GENERATOR << std::endl;
  std::cout<<"  Coordinate system:          " << COORDINATE_SYSTEM << std::endl;
  if (NON_BAROTROPIC_EOS) {
    std::cout<<"  Equation of state:          adiabatic" << std::endl;
  } else {
    std::cout<<"  Equation of state:          isothermal" << std::endl;
  }
  std::cout<<"  Riemann solver:             " << RIEMANN_SOLVER << std::endl;
  std::cout<<"  Reconstruction method:      " << RECONSTRUCTION_METHOD << std::endl;
  std::cout<<"  Hydro integrator:           " << HYDRO_TIME_INTEGRATOR << std::endl;
  std::cout<<"  Compiler and flags:         " << COMPILED_WITH << std::endl;
  if (SHEARING_BOX) {
    std::cout<<"  ShearingBox BCs:            ON" << std::endl;
  } else {
    std::cout<<"  ShearingBox BCs:            OFF" << std::endl;
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
#ifdef HDF5OUTPUT
  std::cout<<"  HDF5 Output:                ON" << std::endl;
#else
  std::cout<<"  HDF5 Output:                OFF" << std::endl;
#endif

  return;
}
