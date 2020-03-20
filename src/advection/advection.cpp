//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file advection.hpp
//  \brief implementation of functions in the Advection class

// C++ headers
// #include <algorithm>  // min()
// #include <cmath>      // fabs(), sqrt()
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/cc/bvals_cc.hpp"

#include "advection.hpp"

// constructor, initializes data structures and parameters

Advection::Advection(MeshBlock *pmb, ParameterInput *pin) :
  pmy_block(pmb),
  u(pmb->ncells3, pmb->ncells2, pmb->ncells1),
  empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
  coarse_u_(pmb->ncc3, pmb->ncc2, pmb->ncc1,
            (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
             AthenaArray<Real>::DataStatus::empty)),
  ubvar(pmb, &u, &coarse_u_, empty_flux)
{
  pmy_block = pmb;
  Coordinates * pco = pmb->pcoord;

  int nc1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int nc2 = 1, nc3 = 1;
  if(pmy_block->block_size.nx2 > 1)
    nc2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if(pmy_block->block_size.nx3 > 1)
    nc3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  Mesh *pm = pmy_block->pmy_mesh;

  // inform MeshBlock that this array is the "primary" representation
  // Used for:
  // (1) load-balancing
  // (2) (future) dumping to restart file
  pmb->RegisterMeshBlockData(u);

  // Allocate memory for the solution and its time derivative
  u.NewAthenaArray(nc3, nc2, nc1);
  u1.NewAthenaArray(nc3, nc2, nc1);

  rhs.NewAthenaArray(nc3, nc2, nc1);

  exact.NewAthenaArray(nc3, nc2, nc1);
  error.NewAthenaArray(nc3, nc2, nc1);

  // If user-requested time integrator is type 3S* allocate additional memory
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4")
    u2.NewAthenaArray(nc3, nc2, nc1);

  cx1 = pin->GetOrAddReal("problem", "cx1", 0);
  cx2 = pin->GetOrAddReal("problem", "cx2", 0);
  cx3 = pin->GetOrAddReal("problem", "cx3", 0);
  use_Sommerfeld = pin->GetOrAddInteger("advection", "use_Sommerfeld", 0);

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmy_block->pmr->AddToRefinement(&u, &coarse_u_);
  }

  // enroll CellCenteredBoundaryVariable object
  ubvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&ubvar);
  pmb->pbval->bvars_main_int.push_back(&ubvar);

  // Allocate memory for scratch arrays
  dt1_.NewAthenaArray(nc1);
  dt2_.NewAthenaArray(nc1);
  dt3_.NewAthenaArray(nc1);

  // Set up finite difference operators
  FD.stride[0] = 1;
  FD.stride[1] = 0;
  FD.stride[2] = 0;
  FD.idx[0] = 1.0 / pco->dx1v(0);
  FD.idx[1] = 0.0;
  FD.idx[2] = 0.0;
  if(nc2 > 1) {
    FD.stride[1] = nc1;
    FD.idx[1] = 1.0 / pco->dx2v(0);
  }
  if(nc3 > 1) {
    FD.stride[2] = nc2 * nc1;
    FD.idx[2] = 1.0 / pco->dx3v(0);
  }

}


// destructor
Advection::~Advection()
{
  u.DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  u2.DeleteAthenaArray(); // only allocated in case of 3S*-type of integrator

  rhs.DeleteAthenaArray();

  exact.DeleteAthenaArray();
  error.DeleteAthenaArray();

}


//----------------------------------------------------------------------------------------
//! \fn  void Advection::AddAdvectionRHS
//  \brief Adds RHS to weighted average of variables from
//  previous step(s) of time integrator algorithm

void Advection::AddAdvectionRHS(const Real wght, AthenaArray<Real> &u_out) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // update variables
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        u_out(k, j, i) += wght*(pmb->pmy_mesh->dt)*rhs(k, j, i);
      }
    }
  }

  return;
}
