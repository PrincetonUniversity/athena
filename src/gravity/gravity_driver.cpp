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


//! \fn GravityDriver::GravityDriver(Mesh *pm, ParameterInput *pin)
//  \brief GravityDriver constructor

GravityDriver::GravityDriver(Mesh *pm, ParameterInput *pin)
 : MultigridDriver(Mesh *pm, ParameterInput *pin)
{
  Real dx=(pm->mesh_size.x1max-pm->mesh_size.x1min)/pm->mesh_size.nx1;
  mgroot_ = new MGGravity(pm->nrbx1, pm->nrbx2, pm->nrbx3, 1, dx);
}


//! \fn Multigrid* GravityDriver::GetMultigridBlock (MeshBlock *pmb)
//  \brief returns a pointer to the multigrid gravity object

Multigrid* GravityDriver::GetMultigridBlock (MeshBlock *pmb)
{
  return pmb->pmggrav;
}


