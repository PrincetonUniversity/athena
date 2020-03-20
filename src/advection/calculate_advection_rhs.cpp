//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_advection_rhs.cpp
//  \brief Calculate advection equation RHS

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "advection.hpp"

//! \fn void Advection::AdvectionRHS
//  \brief Calculate RHS for the advection equation using finite-differencing
void Advection::AdvectionRHS(AthenaArray<Real> & u){
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // internal dimension inferred
  AthenaArray<Real> wu;
  wu.InitWithShallowSlice(u, 0, 1);

  Real cxn[3];
  cxn[0] = pmb->padv->cx1;
  cxn[1] = pmb->padv->cx2;
  cxn[2] = pmb->padv->cx3;

  for(int k = ks; k <= ke; ++k) {
    for(int j = js; j <= je; ++j) {
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        rhs(k,j,i) = 0.0;
      }
      for(int a = 0; a < 3; ++a) {
#pragma omp simd
        for(int i = is; i <= ie; ++i) {
          rhs(k,j,i) += -FD.Da_x(a, cxn[a], wu(k,j,i));
          // rhs(k,j,i) += cxn[a] * FD.Dx(a, wu(k,j,i));
        }
      }
    }
  }
}

//! \fn void Advection:AdvectionBoundaryRHS
//   \brief Calculate the boundary RHS
void Advection::AdvectionBoundaryRHS(AthenaArray<Real> &u){
  MeshBlock * pmb = pmy_block;
  // add Sommerfeld condition here if required
  return;
}
