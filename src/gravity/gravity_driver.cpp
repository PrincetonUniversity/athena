//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity_driver.cpp
//  \brief implementation of functions in class GravityDriver

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

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;



//----------------------------------------------------------------------------------------
//! \fn GravityDriver::GravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary,
//                                   ParameterInput *pin)
//  \brief GravityDriver constructor

GravityDriver::GravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary, ParameterInput *pin)
 : MultigridDriver(pm, MGBoundary, 1)
{
  four_pi_G_=pin->GetOrAddReal("problem", "four_pi_G", 1.0); // default: 4piG=1
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid* GravityDriver::AllocateNewMultigrid(RegionSize isize,
//                 MGBoundaryFunc_t *MGBoundary, enum BoundaryFlag *input_bcs, bool root)
//  \brief Allocate a MGGravity object
Multigrid* GravityDriver::AllocateNewMultigrid(RegionSize isize,
           MGBoundaryFunc_t *MGBoundary, enum BoundaryFlag *input_bcs, bool root = false)
{
  return new MGGravity(this, isize, MGBoundary, input_bcs, root);
}


//----------------------------------------------------------------------------------------
//! \fn void GravityDriver::LoadSourceAndData(void)
//  \brief load the sourterm and initial guess (if needed)

void GravityDriver::LoadSourceAndData(void)
{
  AthenaArray<Real> tmp;
  Multigrid *pmggrav=pmg_;
  while(pmggrav!=NULL) {
    pmggrav->LoadSource(tmp, IDN, NGHOST, four_pi_G_);
    if(mode_>=2) // iterative mode - load initial guess
      pmggrav->LoadFinestData(tmp, 0, NGHOST);
    pmggrav=pmggrav->next;
  }
  return;
}

