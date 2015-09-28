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
//! \file dc.cpp
//  \brief piecewise constant (donor cell) reconstruction
//======================================================================================

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../hydro.hpp"
#include "../../../mesh.hpp"

// this class header
#include "../hydro_integrator.hpp"

//--------------------------------------------------------------------------------------
//! \fn HydroIntegrator::DonorCellX1()
//  \brief 

void HydroIntegrator::DonorCellX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(n,i) = w(n,k,j,i-1);
      wr(n,i) = w(n,k,j,i  );
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(IBY,i) = bcc(IB2,k,j,i-1);
      wl(IBZ,i) = bcc(IB3,k,j,i-1);
      wr(IBY,i) = bcc(IB2,k,j,i  );
      wr(IBZ,i) = bcc(IB3,k,j,i  );
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn HydroIntegrator::DonorCellX2()
//  \brief 

void HydroIntegrator::DonorCellX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(n,i) = w(n,k,j-1,i);
      wr(n,i) = w(n,k,j  ,i);
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(IBY,i) = bcc(IB3,k,j-1,i);
      wl(IBZ,i) = bcc(IB1,k,j-1,i);
      wr(IBY,i) = bcc(IB3,k,j  ,i);
      wr(IBZ,i) = bcc(IB1,k,j  ,i);
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn HydroIntegrator::DonorCellX3()
//  \brief 

void HydroIntegrator::DonorCellX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(n,i) = w(n,k-1,j,i);
      wr(n,i) = w(n,k  ,j,i);
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(IBY,i) = bcc(IB1,k-1,j,i);
      wl(IBZ,i) = bcc(IB2,k-1,j,i);
      wr(IBY,i) = bcc(IB1,k  ,j,i);
      wr(IBZ,i) = bcc(IB2,k  ,j,i);
    }
  }

  return;
}
