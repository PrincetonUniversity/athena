//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity_driver.cpp
//  \brief implementation of functions in class GravityDriver

// Athena++ headers
#include "mggravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../multigrid/multigrid.hpp"


//! \fn GravityDriver::GravityDriver(Mesh *pm, MeshBlock *pblock,
//                                   MGBoundaryFunc_t *MGBoundary, ParameterInput *pin)
//  \brief GravityDriver constructor

GravityDriver::GravityDriver(Mesh *pm, MeshBlock *pblock, MGBoundaryFunc_t *MGBoundary,
                             ParameterInput *pin)
 : MultigridDriver(pm, pblock, MGBoundary, pin)
{
  mgroot_ = new MGGravity(pm->nrbx1, pm->nrbx2, pm->nrbx3, pm->mesh_size, MGBoundary);
}


//! \fn Multigrid* GravityDriver::GetMultigridBlock (MeshBlock *pmb)
//  \brief returns a pointer to the multigrid gravity object

Multigrid* GravityDriver::GetMultigridBlock (MeshBlock *pmb)
{
  return pmb->pmggrav;
}


