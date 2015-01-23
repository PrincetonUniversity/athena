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
#include "../../fluid.hpp"          // Fluid
#include "../../../mesh.hpp"           // Block

//======================================================================================
//! \file plm.cpp
//  \brief  piecewise linear reconstruction
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::ReconstructionFuncX1()
//  \brief 

void FluidIntegrator::ReconstructionFuncX1(const int n, const int m, const int k,
  const int j, const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;

#pragma simd
  for (int i=is; i<=(ie+1); ++i){
    const Real& q_im2 = q(n,k,j,i-2);
    const Real& q_im1 = q(n,k,j,i-1);
    const Real& q_i   = q(n,k,j,i  );
    const Real& q_ip1 = q(n,k,j,i+1);
    Real& dx_im2 = pmy_fluid->pmy_block->dx1v(i-2);
    Real& dx_im1 = pmy_fluid->pmy_block->dx1v(i-1);
    Real& dx_i   = pmy_fluid->pmy_block->dx1v(i);

    Real dql = (q_im1 - q_im2)/dx_im2;
    Real dqc = (q_i   - q_im1)/dx_im1;
    Real dqr = (q_ip1 - q_i  )/dx_i;

    // Apply monotonicity constraints, compute ql_(i-1/2)
    Real dq2 = dql*dqc;
    Real dqm = (dq2 > 0.0) ? dqm = dq2/(dql + dqc) : 0.0;

    Real& dxf_im1 = pmy_fluid->pmy_block->dx1f(i-1);
    ql(m,i) = q_im1 + dxf_im1*dqm;
    
    // Apply monotonicity constraints, compute qr_(i-1/2)
    dq2 = dqc*dqr;
    dqm = (dq2 > 0.0) ? dqm = dq2/(dqc + dqr) : 0.0;

    Real& dxf_i = pmy_fluid->pmy_block->dx1f(i);
    qr(m,i) = q_i   - dxf_i*dqm;
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::ReconstructionFuncX2()
//  \brief 

void FluidIntegrator::ReconstructionFuncX2(const int n, const int m, const int k,
  const int j, const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;
  Real dx2jm2i = 1.0/pmy_fluid->pmy_block->dx2v(j-2);
  Real dx2jm1i = 1.0/pmy_fluid->pmy_block->dx2v(j-1);
  Real dx2ji   = 1.0/pmy_fluid->pmy_block->dx2v(j);

#pragma simd
  for (int i=is; i<=ie; ++i){
    const Real& q_jm2 = q(n,k,j-2,i);
    const Real& q_jm1 = q(n,k,j-1,i);
    const Real& q_j   = q(n,k,j  ,i);
    const Real& q_jp1 = q(n,k,j+1,i);

    Real dql = (q_jm1 - q_jm2)*dx2jm2i;
    Real dqc = (q_j   - q_jm1)*dx2jm1i;
    Real dqr = (q_jp1 - q_j  )*dx2ji;

    // Apply monotonicity constraints, compute ql_(i-1/2)
    Real dq2 = dql*dqc;
    Real dqm = (dq2 > 0.0) ? dqm = dq2/(dql + dqc) : 0.0;

    ql(m,i) = q_jm1 + (pmy_fluid->pmy_block->dx2f(j-1))*dqm;
    
    // Apply monotonicity constraints, compute qr_(i-1/2)
    dq2 = dqc*dqr;
    dqm = (dq2 > 0.0) ? dqm = dq2/(dqc + dqr) : 0.0;

    qr(m,i) = q_j   - (pmy_fluid->pmy_block->dx2f(j))*dqm;
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::ReconstructionFuncX3()
//  \brief 

void FluidIntegrator::ReconstructionFuncX3(const int n, const int m, const int k,
  const int j, const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;
  Real dx3km2i = 1.0/pmy_fluid->pmy_block->dx3v(k-2);
  Real dx3km1i = 1.0/pmy_fluid->pmy_block->dx3v(k-1);
  Real dx3ki   = 1.0/pmy_fluid->pmy_block->dx3v(k);

#pragma simd
  for (int i=is; i<=ie; ++i){
    const Real& q_km2 = q(n,k-2,j,i);
    const Real& q_km1 = q(n,k-1,j,i);
    const Real& q_k   = q(n,k  ,j,i);
    const Real& q_kp1 = q(n,k+1,j,i);

    Real dql = (q_km1 - q_km2)*dx3km2i;
    Real dqc = (q_k   - q_km1)*dx3km1i;
    Real dqr = (q_kp1 - q_k  )*dx3ki;

    // Apply monotonicity constraints, compute ql_(i-1/2)
    Real dq2 = dql*dqc;
    Real dqm = (dq2 > 0.0) ? dqm = dq2/(dql + dqc) : 0.0;

    ql(m,i) = q_km1 + (pmy_fluid->pmy_block->dx3f(k-1))*dqm;
    
    // Apply monotonicity constraints, compute qr_(i-1/2)
    dq2 = dqc*dqr;
    dqm = (dq2 > 0.0) ? dqm = dq2/(dqc + dqr) : 0.0;

    qr(m,i) = q_k   - (pmy_fluid->pmy_block->dx3f(k))*dqm;
  }

  return;
}
