//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_rhs.cpp
//  \brief Calculate wave equation RHS

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "wave.hpp"

//! \fn void Wave::WaveRHS
//  \brief Calculate RHS for the wave equation using finite-differencing
void Wave::WaveRHS(AthenaArray<Real> & u){
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> wu, wpi;
  wu.InitWithShallowSlice(u, NWAVE_CPT, 0, 1);
  wpi.InitWithShallowSlice(u, NWAVE_CPT, 1, 1);

  for(int k = ks; k <= ke; ++k) {
    for(int j = js; j <= je; ++j) {
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        rhs(0,k,j,i) = wpi(k,j,i);
        rhs(1,k,j,i) = 0.0;
      }
      for(int a = 0; a < 3; ++a) {
#pragma omp simd
        for(int i = is; i <= ie; ++i) {
          rhs(1,k,j,i) += FD.Dxx(a, wu(k,j,i));
        }
      }
    }
  }
}

//! \fn void Wave:WaveBoundaryRHS
//   \brief Calculate the boundary RHS
void Wave::WaveBoundaryRHS(AthenaArray<Real> & u){
  MeshBlock * pmb = pmy_block;

  if(pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow) {
    printf("inner_x1 Sommerfeld activated\n");
    WaveSommerfeld_(u, pmb->is, pmb->is, pmb->js, pmb->je, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow) {
    WaveSommerfeld_(u, pmb->ie, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::outflow) {
    WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->js, pmb->js, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::outflow) {
    WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->je, pmb->je, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::outflow) {
    WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ks);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::outflow) {
    WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ke, pmb->ke);
  }
}

void Wave::WaveSommerfeld_(AthenaArray<Real> & u,
                           int const is, int const ie,
                           int const js, int const je,
                           int const ks, int const ke){
  printf("Sommerfeld_\n");
  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, NWAVE_CPT, 1, 1);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  for(int k = ks; k <= ke; ++k) {
    for(int j = js; j <= je; ++j) {
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        // Derivatives of pi
        Real const wpi_x = FD.Ds(0, wpi(k,j,i));
        Real const wpi_y = FD.Ds(1, wpi(k,j,i));
        Real const wpi_z = FD.Ds(2, wpi(k,j,i));

        Real const rr = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(j)) + SQR(pco->x3v(k)));
        Real const sx = pco->x1v(i)/rr;
        Real const sy = pco->x2v(j)/rr;
        Real const sz = pco->x3v(k)/rr;

        rhs(1,k,j,i) = -(wpi(k,j,i)/rr + sx*wpi_x + sy*wpi_y + sz*wpi_z);
      }
    }
  }
}
