//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "../integrators.hpp"

// Athena headers
#include "../../../athena.hpp"         // macros, Real
#include "../../../athena_arrays.hpp"  // AthenaArray
#include "../../fluid.hpp"          // Fluid
#include "../../../mesh.hpp"           // Block

//======================================================================================
/*! \file plm.cpp
 *  \brief  piecewise linear reconstruction
 *====================================================================================*/

//! \fn FluidIntegrator::ReconstructionFuncX1()
//  \brief 

void FluidIntegrator::ReconstructionFuncX1(
  const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, AthenaArray<Real> *pwl, AthenaArray<Real> *pwr)
{
  for (int n=0; n<NFLUID; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){

      Real dwl = (w(n,k,j,i-1) - w(n,k,j,i-2))/pmy_fluid->pmy_block->dx1v(i-2);
      Real dwc = (w(n,k,j,i)   - w(n,k,j,i-1))/pmy_fluid->pmy_block->dx1v(i-1);
      Real dwr = (w(n,k,j,i+1) - w(n,k,j,i)  )/pmy_fluid->pmy_block->dx1v(i  );

// Apply monotonicity constraints to differences in primitive vars, compute wl_(i-1/2)

      Real dw2 = dwl*dwc;
      Real dwm = dw2/(dwl + dwc + TINY_NUMBER);
      if (dw2 <= 0.0) dwm  = 0.0;

      Real& wli = (*pwl)(n,i);
      wli = w(n,k,j,i-1) + (pmy_fluid->pmy_block->dx1f(i-1))*dwm;
    
// Apply monotonicity constraints to differences in primitive vars, compute wr_(i-1/2)

      dw2 = dwc*dwr;
      dwm = dw2/(dwc + dwr + TINY_NUMBER);
      if (dw2 <= 0.0) dwm  = 0.0;

      Real& wri = (*pwr)(n,i);
      wri = w(n,k,j,i) - (pmy_fluid->pmy_block->dx1f(i))*dwm;
    }
  }

  return;
}

//! \fn FluidIntegrator::ReconstructionFuncX2()
//  \brief 

void FluidIntegrator::ReconstructionFuncX2(
  const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, AthenaArray<Real> *pwl, AthenaArray<Real> *pwr)
{
  for (int n=0; n<NFLUID; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){

      Real dwl = (w(n,k,j-1,i) - w(n,k,j-2,i))/pmy_fluid->pmy_block->dx2v(j-2);
      Real dwc = (w(n,k,j,i)   - w(n,k,j-1,i))/pmy_fluid->pmy_block->dx2v(j-1);
      Real dwr = (w(n,k,j+1,i) - w(n,k,j,i)  )/pmy_fluid->pmy_block->dx2v(j  );

// Apply monotonicity constraints to differences in primitive vars, compute wl_(i-1/2)

      Real dw2 = dwl*dwc;
      Real dwm = dw2/(dwl + dwc + TINY_NUMBER);
      if (dw2 <= 0.0) dwm  = 0.0;

      Real& wli = (*pwl)(n,i);
      wli = w(n,k,j-1,i) + (pmy_fluid->pmy_block->dx2f(j-1))*dwm;
    
// Apply monotonicity constraints to differences in primitive vars, compute wr_(i-1/2)

      dw2 = dwc*dwr;
      dwm = dw2/(dwc + dwr + TINY_NUMBER);
      if (dw2 <= 0.0) dwm  = 0.0;

      Real& wri = (*pwr)(n,i);
      wri = w(n,k,j,i) - (pmy_fluid->pmy_block->dx2f(j))*dwm;
    }
  }

  return;
}

//! \fn FluidIntegrator::ReconstructionFuncX3()
//  \brief 

void FluidIntegrator::ReconstructionFuncX3(
  const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, AthenaArray<Real> *pwl, AthenaArray<Real> *pwr)
{
  for (int n=0; n<NFLUID; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){

      Real dwl = (w(n,k-1,j,i) - w(n,k-2,j,i))/pmy_fluid->pmy_block->dx3v(k-2);
      Real dwc = (w(n,k,j,i)   - w(n,k-1,j,i))/pmy_fluid->pmy_block->dx3v(k-1);
      Real dwr = (w(n,k+1,j,i) - w(n,k,j,i)  )/pmy_fluid->pmy_block->dx3v(k  );

// Apply monotonicity constraints to differences in primitive vars, compute wl_(i-1/2)

      Real dw2 = dwl*dwc;
      Real dwm = dw2/(dwl + dwc + TINY_NUMBER);
      if (dw2 <= 0.0) dwm  = 0.0;

      Real& wli = (*pwl)(n,i);
      wli = w(n,k-1,j,i) + (pmy_fluid->pmy_block->dx3f(k-1))*dwm;
    
// Apply monotonicity constraints to differences in primitive vars, compute wr_(i-1/2)

      dw2 = dwc*dwr;
      dwm = dw2/(dwc + dwr + TINY_NUMBER);
      if (dw2 <= 0.0) dwm  = 0.0;

      Real& wri = (*pwr)(n,i);
      wri = w(n,k,j,i) - (pmy_fluid->pmy_block->dx3f(k))*dwm;
    }
  }

  return;
}
