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
  // internal dimension inferred
  wu.InitWithShallowSlice(u, 0, 1);
  wpi.InitWithShallowSlice(u, 1, 1);

  Real c_2 = SQR(pmb->pwave->c);

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
          rhs(1,k,j,i) += c_2 * FD.Dxx(a, wu(k,j,i));
        }
      }
    }
  }
}

//! \fn void Wave:WaveBoundaryRHS
//   \brief Calculate the boundary RHS
void Wave::WaveBoundaryRHS(AthenaArray<Real> & u){
  MeshBlock * pmb = pmy_block;

  if (use_Sommerfeld == 0)
    return;

  // the appropriate condition to apply depends on dimension
  // an assumption is made that any present Sommerfeld conditions are order as
  // x1, x2, ...

  if(pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::extrapolate_outflow ||
     pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow) {
    switch(use_Sommerfeld) {
    case 1:
      WaveSommerfeld_1d_L_(u, pmb->is, pmb->is, pmb->js, pmb->je, pmb->ks, pmb->ke);
      break;
    case 2:
      break;
    case 3:
      WaveSommerfeld_3d_(u, pmb->is, pmb->is, pmb->js, pmb->je, pmb->ks, pmb->ke);
      break;
    }
    // WaveSommerfeld_(u, pmb->is, pmb->is, pmb->js, pmb->je, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::extrapolate_outflow ||
     pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow) {
    switch(use_Sommerfeld) {
    case 1:
      WaveSommerfeld_1d_R_(u, pmb->ie, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ke);
      break;
    case 2:
      break;
    case 3:
      WaveSommerfeld_3d_(u, pmb->ie, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ke);
      break;
    }
    // WaveSommerfeld_(u, pmb->ie, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::extrapolate_outflow ||
     pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::outflow) {
    switch(use_Sommerfeld) {
    case 2:
      break;
    case 3:
      WaveSommerfeld_3d_(u, pmb->is, pmb->ie, pmb->js, pmb->js, pmb->ks, pmb->ke);
      break;
    }
    // WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->js, pmb->js, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::extrapolate_outflow ||
     pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::outflow) {
    switch(use_Sommerfeld) {
    case 2:
      break;
    case 3:
      WaveSommerfeld_3d_(u, pmb->is, pmb->ie, pmb->je, pmb->je, pmb->ks, pmb->ke);
      break;
    }
    // WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->je, pmb->je, pmb->ks, pmb->ke);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::extrapolate_outflow ||
     pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::outflow) {
    if(use_Sommerfeld == 3)
      WaveSommerfeld_3d_(u, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ks);
    // WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ks);
  }
  if(pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::extrapolate_outflow ||
     pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::outflow) {
    if(use_Sommerfeld == 3)
      WaveSommerfeld_3d_(u, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ke, pmb->ke);
    // WaveSommerfeld_(u, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ke, pmb->ke);
  }

}


void Wave::WaveSommerfeld_1d_L_(AthenaArray<Real> & u,
                                int const is, int const ie,
                                int const js, int const je,
                                int const ks, int const ke){

  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)
  //
  // In particular, in 1d, these are exact.

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;
  Real c = pmb->pwave->c;


  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        Real const wpi_x = FD.Ds(0, wpi(k,j,i));
        rhs(1,k,j,i) = pmb->pwave->c * wpi_x;
      }
}


void Wave::WaveSommerfeld_1d_R_(AthenaArray<Real> & u,
                                int const is, int const ie,
                                int const js, int const je,
                                int const ks, int const ke){

  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)
  //
  // In particular, in 1d, these are exact.

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;
  Real c = pmb->pwave->c;


  for(int k = ks; k <= ke; ++k)
    for(int j = js; j <= je; ++j)
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        Real const wpi_x = FD.Ds(0, wpi(k,j,i));
        rhs(1,k,j,i) = -pmb->pwave->c * wpi_x;
      }

}

void Wave::WaveSommerfeld_2d_(AthenaArray<Real> & u,
                              int const is, int const ie,
                              int const js, int const je,
                              int const ks, int const ke){

  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;
  Real c = pmb->pwave->c;
  for(int k = ks; k <= ke; ++k) {
    for(int j = js; j <= je; ++j) {
#pragma omp simd
      for(int i = is; i <= ie; ++i) {
        // Derivatives of pi
        Real const wpi_x = FD.Ds(0, wpi(k,j,i));
        Real const wpi_y = FD.Ds(1, wpi(k,j,i));

        Real const rr = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(j)));
        Real const sx = pco->x1v(i);
        Real const sy = pco->x2v(j);

        rhs(1,k,j,i) = -pmb->pwave->c * (wpi(k,j,i)/rr + 2 * (sx*wpi_x + sy*wpi_y));
      }
    }
  }
}


void Wave::WaveSommerfeld_3d_(AthenaArray<Real> & u,
                              int const is, int const ie,
                              int const js, int const je,
                              int const ks, int const ke){

  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;
  Real c = pmb->pwave->c;
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

        rhs(1,k,j,i) = -pmb->pwave->c * (wpi(k,j,i)/rr + sx*wpi_x + sy*wpi_y + sz*wpi_z);
      }
    }
  }
}
