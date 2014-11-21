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

// C++ headers
#include <iostream>    // endl, ostream
#include <sstream>     // stringstream
#include <stdexcept>   // runtime_error
#include <sys/stat.h>  // mkdir()
#include <unistd.h>    // chdir()

// Athena headers
#include "athena.hpp"

//======================================================================================
//! \file utils.cpp 
//  \brief contains a variety of utility functions
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
  std::cout<<"  Fluid integrator:           " << FLUID_TIME_INTEGRATOR << std::endl;
  std::cout<<"  Compiler and flags:         " << COMPILED_WITH << std::endl;
  if (MAGNETIC_FIELDS_ENABLED) {
    std::cout<<"  Magnetic fields:            enabled" << std::endl;
  } else {
    std::cout<<"  Magnetic fields:            disabled" << std::endl;
  }
  if (RELATIVISTIC_DYNAMICS) {
    std::cout<<"  Relativistic dynamics:      enabled " << std::endl;
  } else {
    std::cout<<"  Relativistic dynamics:      disabled " << std::endl;
  }
#ifdef OPENMP_PARALLEL
  std::cout<<"  OpenMP parallelism:         enabled" << std::endl;
#else
  std::cout<<"  OpenMP parallelism:         disabled" << std::endl;
#endif


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ChangeToRunDir(const char *pdir)
//  \brief change to input run directory; create if it does not exist yet

void ChangeToRunDir(const char *pdir)
{
  std::stringstream msg;

  if (pdir == NULL || *pdir == '\0') return;

  mkdir(pdir, 0775);
  if(chdir(pdir)) {
    msg << "### FATAL ERROR in function [ChangeToRunDir]" << std::endl
        << "Cannot cd to directory '" << pdir << "'";
    throw std::runtime_error(msg.str().c_str());
  }

  return;
}
