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
//! \fn GravityDriver::GravityDriver(Mesh *pm, MeshBlock *pblock,
//                                   MGBoundaryFunc_t *MGBoundary, ParameterInput *pin)
//  \brief GravityDriver constructor

GravityDriver::GravityDriver(Mesh *pm, MeshBlock *pblock, MGBoundaryFunc_t *MGBoundary,
                             ParameterInput *pin)
 : MultigridDriver(pm, pblock, MGBoundary, 1, pin)
{
  mgroot_ = new MGGravity(NULL,pm->nrbx1,pm->nrbx2,pm->nrbx3,pm->mesh_size,MGBoundary);
  four_pi_G_=pin->GetOrAddReal("problem", "four_pi_G", 1.0); // default: 4piG=1
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid* GravityDriver::GetMultigridBlock(MeshBlock *pmb)
//  \brief returns a pointer to the multigrid gravity object

Multigrid* GravityDriver::GetMultigridBlock(MeshBlock *pmb)
{
  return pmb->pmggrav;
}


//----------------------------------------------------------------------------------------
//! \fn void GravityDriver::LoadSourceAndData(void)
//  \brief load the sourterm and initial guess (if needed)

void GravityDriver::LoadSourceAndData(void)
{
  MeshBlock *pb=pblock_;
  while(pb!=NULL) {
    pb->pmggrav->LoadSource(pb->phydro->u, IDN, NGHOST, four_pi_G_);
    if(mode_>=2) // iterative mode - load initial guess
      pb->pmggrav->LoadFinestData(pb->pgrav->phi, 0, NGHOST);
    pb=pb->next;
  }
  return;
}

