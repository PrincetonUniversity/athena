//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_gravity.cpp
//  \brief create multigrid solver for gravity

// C headers

// C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "gravity.hpp"
#include "mg_gravity.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;

//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::MGGravityDriver(Mesh *pm, ParameterInput *pin)
//  \brief MGGravityDriver constructor

MGGravityDriver::MGGravityDriver(Mesh *pm, ParameterInput *pin)
    : MultigridDriver(pm, pm->MGGravityBoundaryFunction_, 1) {
}


//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::~MGGravityDriver()
//  \brief MGGravityDriver destructor
MGGravityDriver::~MGGravityDriver() {
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::Solve(int stage)
//  \brief load the data and solve

void MGGravityDriver::Solve(int stage) {
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(int color)
//  \brief Red-Black Gauss-Seidel Smoother
void MGGravity::Smooth(int color) {
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravity::CalculateDefect()
//  \brief calculate the residual

void MGGravity::CalculateDefect() {
  return;
}
