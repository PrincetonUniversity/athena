//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_gravity.cpp
//! \brief create multigrid solver for gravity

// C headers

// C++ headers
#include <algorithm>
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
//! \brief MGGravityDriver constructor

MGGravityDriver::MGGravityDriver(Mesh *pm, ParameterInput *pin)
    : MultigridDriver(pm, pm->MGGravityBoundaryFunction_, 1) {
}


//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::~MGGravityDriver()
//! \brief MGGravityDriver destructor

MGGravityDriver::~MGGravityDriver() {
}


//----------------------------------------------------------------------------------------
//! \fn MGGravity::MGGravity(MultigridDriver *pmd, MeshBlock *pmb)
//! \brief MGGravity constructor

MGGravity::MGGravity(MultigridDriver *pmd, MeshBlock *pmb) : Multigrid(pmd, pmb, 1, 1) {
}


//----------------------------------------------------------------------------------------
//! \fn MGGravity::~MGGravity()
//! \brief MGGravity deconstructor

MGGravity::~MGGravity() {
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::Solve(int stage)
//! \brief load the data and solve

void MGGravityDriver::Solve(int stage) {
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src,
//!           int rlev, int il, int iu, int jl, int ju, int kl, int ku, int color)
//! \brief Implementation of the Red-Black Gauss-Seidel Smoother
//!        rlev = relative level from the finest level of this Multigrid block

void MGGravity::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src, int rlev,
                       int il, int iu, int jl, int ju, int kl, int ku, int color) {
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::CalculateDefect(AthenaArray<Real> &def,
//!                      const AthenaArray<Real> &u, const AthenaArray<Real> &src,
//!                      int rlev, int il, int iu, int jl, int ju, int kl, int ku)
//! \brief Implementation of the Defect calculation
//!        rlev = relative level from the finest level of this Multigrid block

void MGGravity::CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
                                const AthenaArray<Real> &src, int rlev,
                                int il, int iu, int jl, int ju, int kl, int ku) {
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::CalculateFASRHS(AthenaArray<Real> &src,
//!  const AthenaArray<Real> &u, int rlev, int il, int iu, int jl, int ju, int kl, int ku)
//! \brief Implementation of the RHS calculation for FAS
//!        rlev = relative level from the finest level of this Multigrid block

void MGGravity::CalculateFASRHS(AthenaArray<Real> &src, const AthenaArray<Real> &u,
                         int rlev, int il, int iu, int jl, int ju, int kl, int ku) {
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::ProlongateOctetBoundariesFluxCons(AthenaArray<Real> &dst)
//! \brief prolongate octet boundaries using the flux conservation formula

void MGGravityDriver::ProlongateOctetBoundariesFluxCons(AthenaArray<Real> &dst) {
}
