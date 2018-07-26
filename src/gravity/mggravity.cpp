//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mggravity.cpp
//  \brief create multigrid solver for gravity

// C/C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Athena++ headers
#include "mggravity.hpp"
#include "gravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../multigrid/multigrid.hpp"
#include "../globals.hpp"

class MeshBlock;

//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::MGGravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary,
//                                   ParameterInput *pin)
//  \brief MGGravityDriver constructor

MGGravityDriver::MGGravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary,
                                 ParameterInput *pin)
    : MultigridDriver(pm, MGBoundary, 1) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::Solve(int stage)
//  \brief load the data and solve

void MGGravityDriver::Solve(int stage) {
}

//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(int color)
//  \brief Red-Black Gauss-Seidel Smoother
void MGGravity::Smooth(int color) {
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravity::CalculateDefect(void)
//  \brief calculate the residual

void MGGravity::CalculateDefect(void) {
}
