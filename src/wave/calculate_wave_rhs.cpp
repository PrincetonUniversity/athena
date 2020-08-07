//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_wave_rhs.cpp
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

  AthenaArray<Real> wu, wpi;
  // internal dimension inferred
  wu.InitWithShallowSlice(u, 0, 1);
  wpi.InitWithShallowSlice(u, 1, 1);

  Real c_2 = SQR(c);

  for(int k = mbi.kl; k <= mbi.ku; ++k) {
    for(int j = mbi.jl; j <= mbi.ju; ++j) {
#pragma omp simd
      for(int i = mbi.il; i <= mbi.iu; ++i) {
        rhs(0,k,j,i) = wpi(k,j,i);
        rhs(1,k,j,i) = 0.0;
      }
      for(int a = 0; a < 3; ++a) {
#pragma omp simd
        for(int i = mbi.il; i <= mbi.iu; ++i) {
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
  if (use_Dirichlet) {

    if(pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::extrapolate_outflow ||
       pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow)
        WaveBoundaryDirichlet_(u, mbi.il, mbi.il, mbi.jl, mbi.ju, mbi.kl, mbi.ku);
    if(pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::extrapolate_outflow ||
       pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow)
        WaveBoundaryDirichlet_(u, mbi.iu, mbi.iu, mbi.jl, mbi.ju, mbi.kl, mbi.ku);

    if(pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::extrapolate_outflow ||
       pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::outflow)
        WaveBoundaryDirichlet_(u, mbi.il, mbi.iu, mbi.jl, mbi.jl, mbi.kl, mbi.ku);
    if(pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::extrapolate_outflow ||
       pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::outflow)
        WaveBoundaryDirichlet_(u, mbi.il, mbi.iu, mbi.ju, mbi.ju, mbi.kl, mbi.ku);

    if(pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::extrapolate_outflow ||
       pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::outflow)
        WaveBoundaryDirichlet_(u, mbi.il, mbi.iu, mbi.jl, mbi.ju, mbi.kl, mbi.kl);
    if(pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::extrapolate_outflow ||
       pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::outflow)
        WaveBoundaryDirichlet_(u, mbi.il, mbi.iu, mbi.jl, mbi.ju, mbi.ku, mbi.ku);

  } else if (use_Sommerfeld) {

    if (pmb->pmy_mesh->ndim == 3) {
      if(pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow)
          WaveSommerfeld_3d_(u, mbi.il, mbi.il, mbi.jl, mbi.ju, mbi.kl, mbi.ku);
      if(pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow)
          WaveSommerfeld_3d_(u, mbi.iu, mbi.iu, mbi.jl, mbi.ju, mbi.kl, mbi.ku);

      if(pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::outflow)
          WaveSommerfeld_3d_(u, mbi.il, mbi.iu, mbi.jl, mbi.jl, mbi.kl, mbi.ku);
      if(pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::outflow)
          WaveSommerfeld_3d_(u, mbi.il, mbi.iu, mbi.ju, mbi.ju, mbi.kl, mbi.ku);

      if(pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::outflow)
          WaveSommerfeld_3d_(u, mbi.il, mbi.iu, mbi.jl, mbi.ju, mbi.kl, mbi.kl);
      if(pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::outflow)
          WaveSommerfeld_3d_(u, mbi.il, mbi.iu, mbi.jl, mbi.ju, mbi.ku, mbi.ku);
    } else if (pmb->pmy_mesh->ndim == 2) {
      if(pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow)
          WaveSommerfeld_2d_(u, mbi.il, mbi.il, mbi.jl, mbi.ju, mbi.kl, mbi.ku);
      if(pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow)
          WaveSommerfeld_2d_(u, mbi.iu, mbi.iu, mbi.jl, mbi.ju, mbi.kl, mbi.ku);

      if(pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::outflow)
          WaveSommerfeld_2d_(u, mbi.il, mbi.iu, mbi.jl, mbi.jl, mbi.kl, mbi.ku);
      if(pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::outflow)
          WaveSommerfeld_2d_(u, mbi.il, mbi.iu, mbi.ju, mbi.ju, mbi.kl, mbi.ku);
    } else {
      if(pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow)
          WaveSommerfeld_1d_L_(u, mbi.il, mbi.il, mbi.jl, mbi.ju, mbi.kl, mbi.ku);
      if(pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::extrapolate_outflow ||
        pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow)
          WaveSommerfeld_1d_R_(u, mbi.iu, mbi.iu, mbi.jl, mbi.ju, mbi.kl, mbi.ku);
    }

  } else {
    return;
  }

  return;


}

void Wave::WaveBoundaryDirichlet_(AthenaArray<Real> & u, int il, int iu,
                                  int jl, int ju, int kl, int ku) {
  // u(x; t)|_boundary = 0 for all t
  // => set u_x^n(x; t)|_boundary = 0 for all t

  // BD: for debug
  Real const dconst = (DBG_VC_CONSISTENCY) ? 1.0 : 0.;

  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
#pragma omp simd
      for(int i = il; i <= iu; ++i) {
        rhs(0,k,j,i) = dconst;
        rhs(1,k,j,i) = dconst;
      }
}

void Wave::WaveSommerfeld_1d_L_(AthenaArray<Real> & u,
  int il, int iu, int jl, int ju, int kl, int ku){
  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)
  //
  // In particular, in 1d, these are exact.

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  for(int k = mbi.kl; k <= mbi.ku; ++k)
    for(int j = mbi.jl; j <= mbi.ju; ++j)
#pragma omp simd
      for(int i = mbi.il; i <= mbi.il; ++i) {
        rhs(1,k,j,i) = c * FD.Ds(0, wpi(k,j,i));
      }
}


void Wave::WaveSommerfeld_1d_R_(AthenaArray<Real> & u,
  int il, int iu, int jl, int ju, int kl, int ku){
  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)
  //
  // In particular, in 1d, these are exact.

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  for(int k = kl; k <= ku; ++k)
    for(int j = jl; j <= ju; ++j)
#pragma omp simd
      for(int i = iu; i <= iu; ++i) {
        rhs(1,k,j,i) = -c * FD.Ds(0, wpi(k,j,i));
      }

}

void Wave::WaveSommerfeld_2d_(AthenaArray<Real> & u,
  int il, int iu, int jl, int ju, int kl, int ku){
  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  for(int k = kl; k <= ku; ++k) {
    for(int j = jl; j <= ju; ++j) {
#pragma omp simd
      for(int i = il; i <= iu; ++i) {
        // Derivatives of pi
        Real const wpi_x = FD.Ds(0, wpi(k,j,i));
        Real const wpi_y = FD.Ds(1, wpi(k,j,i));

        Real const rr = std::sqrt(SQR(mbi.x1(i)) + SQR(mbi.x2(j)));
        Real const sx = mbi.x1(i);
        Real const sy = mbi.x2(j);

        rhs(1,k,j,i) = -c * (wpi(k,j,i)/rr + 2. * (sx*wpi_x + sy*wpi_y)) / 2.;
      }
    }
  }
}


void Wave::WaveSommerfeld_3d_(AthenaArray<Real> & u,
  int il, int iu, int jl, int ju, int kl, int ku){
  // For u_tt = c^2 (u_xx1 + .. + u_xxd) with r = sqrt(xx1 ^ 2 + ... + xxd ^2)
  // Sommerfeld conditions should go like:
  // lim r^((d-1)/2)(u_r + 1/c u_t)

  AthenaArray<Real> wpi;
  wpi.InitWithShallowSlice(u, 1, 1);

  for(int k = kl; k <= ku; ++k) {
    for(int j = jl; j <= ju; ++j) {
#pragma omp simd
      for(int i = il; i <= iu; ++i) {
        // Derivatives of pi
        Real const wpi_x = FD.Ds(0, wpi(k,j,i));
        Real const wpi_y = FD.Ds(1, wpi(k,j,i));
        Real const wpi_z = FD.Ds(2, wpi(k,j,i));

        Real const rr = std::sqrt(SQR(mbi.x1(i)) +
                                  SQR(mbi.x2(j)) + SQR(mbi.x3(k)));
        Real const sx = mbi.x1(i)/rr;
        Real const sy = mbi.x2(j)/rr;
        Real const sz = mbi.x3(k)/rr;

        rhs(1,k,j,i) = -c * (wpi(k,j,i)/rr + sx*wpi_x + sy*wpi_y + sz*wpi_z);
      }
    }
  }
}