//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "../fluid_integrator.hpp"

// Athena headers
#include "../../../athena.hpp"         // macros, Real
#include "../../../athena_arrays.hpp"  // AthenaArray
#include "../../fluid.hpp"             // Fluid
#include "../../../mesh.hpp"           // MeshBlock

//======================================================================================
//! \file dc.cpp
//  \brief piecewise constant (donor cell) reconstruction
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::DonorCellX1()
//  \brief 

void FluidIntegrator::DonorCellX1(const int k, const int j,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;

  for (int n=0; n<NFLUID; ++n){
#pragma simd
    for (int i=is; i<=(ie+1); ++i){
      const Real& w_im1 = w(n,k,j,i-1);
      const Real& w_i   = w(n,k,j,i  );

      wl(n,i) = w_im1;
      wr(n,i) = w_i;
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=is; i<=(ie+1); ++i){
      const Real& by_im1 = bcc(IB2,k,j,i-1);
      const Real& bz_im1 = bcc(IB3,k,j,i-1);
      const Real& by_i   = bcc(IB2,k,j,i  );
      const Real& bz_i   = bcc(IB3,k,j,i  );

      wl(IBY,i) = by_im1;
      wl(IBZ,i) = bz_im1;
      wr(IBY,i) = by_i;
      wr(IBZ,i) = bz_i;
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::DonorCellX2()
//  \brief 

void FluidIntegrator::DonorCellX2(const int k, const int j,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;

  for (int n=0; n<NFLUID; ++n){
#pragma simd
    for (int i=is; i<=ie; ++i){
      const Real& w_jm1 = w(n,k,j-1,i);
      const Real& w_j   = w(n,k,j  ,i);

      wl(n,i) = w_jm1;
      wr(n,i) = w_j;
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=is; i<=ie; ++i){
      const Real& bz_jm1 = bcc(IB3,k,j-1,i);
      const Real& bx_jm1 = bcc(IB1,k,j-1,i);
      const Real& bz_j   = bcc(IB3,k,j  ,i);
      const Real& bx_j   = bcc(IB1,k,j  ,i);

      wl(IBY,i) = bz_jm1;
      wl(IBZ,i) = bx_jm1;
      wr(IBY,i) = bz_j;
      wr(IBZ,i) = bx_j;
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::DonorCellX3()
//  \brief 

void FluidIntegrator::DonorCellX3(const int k, const int j,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;

  for (int n=0; n<NFLUID; ++n){
#pragma simd
    for (int i=is; i<=ie; ++i){
      const Real& w_km1 = w(n,k-1,j,i);
      const Real& w_k   = w(n,k  ,j,i);

      wl(n,i) = w_km1;
      wr(n,i) = w_k;
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=is; i<=ie; ++i){
      const Real& bx_km1 = bcc(IB1,k-1,j,i);
      const Real& by_km1 = bcc(IB2,k-1,j,i);
      const Real& bx_k   = bcc(IB1,k  ,j,i);
      const Real& by_k   = bcc(IB2,k  ,j,i);

      wl(IBY,i) = bx_km1;
      wl(IBZ,i) = by_km1;
      wr(IBY,i) = bx_k;
      wr(IBZ,i) = by_k;
    }
  }

  return;
}
