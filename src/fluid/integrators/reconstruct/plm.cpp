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

void FluidIntegrator::PiecewiseLinearX1(const int k, const int j,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;

  for (int n=0; n<NWAVE; ++n) {
#pragma simd
<<<<<<< HEAD
  for (int i=is; i<=(ie+1); ++i){

    Real dql = (q(n,k,j,i-1) - q(n,k,j,i-2))/pmy_fluid->pmy_block->dx1v(i-2);
    Real dqc = (q(n,k,j,i)   - q(n,k,j,i-1))/pmy_fluid->pmy_block->dx1v(i-1);
    Real dqr = (q(n,k,j,i+1) - q(n,k,j,i)  )/pmy_fluid->pmy_block->dx1v(i);

// Apply monotonicity constraints to differences in primitive vars, compute ql_(i-1/2)

    Real dq2 = dql*dqc;
    Real dqm = 0.0;
    if (dq2 > 0.0) dqm = dq2/(dql + dqc);

    ql(m,i) = q(n,k,j,i-1) + (pmy_fluid->pmy_block->dx1f(i-1))*dqm;
=======
    for (int i=is; i<=(ie+1); ++i){
      Real& dx_im2 = pmy_fluid->pmy_block->dx1v(i-2);
      Real& dx_im1 = pmy_fluid->pmy_block->dx1v(i-1);
      Real& dx_i   = pmy_fluid->pmy_block->dx1v(i);

      Real dql,dqr,dqc,q_im1,q_i;
      if (n==NFLUID){
        q_im1 = bcc(IB2,k,j,i-1);
        q_i   = bcc(IB2,k,j,i  );
        dql = (bcc(IB2,k,j,i-1) - bcc(IB2,k,j,i-2))/dx_im2;
        dqc = (bcc(IB2,k,j,i  ) - bcc(IB2,k,j,i-1))/dx_im1;
        dqr = (bcc(IB2,k,j,i+1) - bcc(IB2,k,j,i  ))/dx_i;
      } else if (n==(NFLUID+1)) {
        q_im1 = bcc(IB3,k,j,i-1);
        q_i   = bcc(IB3,k,j,i  );
        dql = (bcc(IB3,k,j,i-1) - bcc(IB3,k,j,i-2))/dx_im2;
        dqc = (bcc(IB3,k,j,i  ) - bcc(IB3,k,j,i-1))/dx_im1;
        dqr = (bcc(IB3,k,j,i+1) - bcc(IB3,k,j,i  ))/dx_i;
      } else {
        q_im1 = q(n,k,j,i-1);
        q_i   = q(n,k,j,i  );
        dql = (q(n,k,j,i-1) - q(n,k,j,i-2))/dx_im2;
        dqc = (q(n,k,j,i  ) - q(n,k,j,i-1))/dx_im1;
        dqr = (q(n,k,j,i+1) - q(n,k,j,i  ))/dx_i;
      }

      // Apply monotonicity constraints, compute ql_(i-1/2)
      Real dq2 = dql*dqc;
      Real dqm = (dq2 > 0.0) ? dqm = dq2/(dql + dqc) : 0.0;

      Real& dxf_im1 = pmy_fluid->pmy_block->dx1f(i-1);
      ql(n,i) = q_im1 + dxf_im1*dqm;
>>>>>>> remotes/origin/master
    
      // Apply monotonicity constraints, compute qr_(i-1/2)
      dq2 = dqc*dqr;
      dqm = (dq2 > 0.0) ? dqm = dq2/(dqc + dqr) : 0.0;

<<<<<<< HEAD
    dq2 = dqc*dqr;
    dqm = 0.0;
    if (dq2 > 0.0) dqm = dq2/(dqc + dqr);

    qr(m,i) = q(n,k,j,i) - (pmy_fluid->pmy_block->dx1f(i))*dqm;
=======
      Real& dxf_i = pmy_fluid->pmy_block->dx1f(i);
      qr(n,i) = q_i   - dxf_i*dqm;
    }
>>>>>>> remotes/origin/master
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::ReconstructionFuncX2()
//  \brief 

void FluidIntegrator::PiecewiseLinearX2(const int k, const int j,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;
  Real dx2jm2i = 1.0/pmy_fluid->pmy_block->dx2v(j-2);
  Real dx2jm1i = 1.0/pmy_fluid->pmy_block->dx2v(j-1);
  Real dx2ji   = 1.0/pmy_fluid->pmy_block->dx2v(j);

  for (int n=0; n<NWAVE; ++n) {
#pragma simd
<<<<<<< HEAD
  for (int i=is; i<=ie; ++i){

    Real dql = (q(n,k,j-1,i) - q(n,k,j-2,i))*dx2jm2i;
    Real dqc = (q(n,k,j,i)   - q(n,k,j-1,i))*dx2jm1i;
    Real dqr = (q(n,k,j+1,i) - q(n,k,j,i)  )*dx2ji;

// Apply monotonicity constraints to differences in primitive vars, compute ql_(i-1/2)

    Real dq2 = dql*dqc;
    Real dqm = 0.0;
    if (dq2 > 0.0) dqm = dq2/(dql + dqc);

    ql(m,i) = q(n,k,j-1,i) + (pmy_fluid->pmy_block->dx2f(j-1))*dqm;
    
// Apply monotonicity constraints to differences in primitive vars, compute qr_(i-1/2)

    dq2 = dqc*dqr;
    dqm = 0.0;
    if (dq2 > 0.0) dqm = dq2/(dqc + dqr);

    qr(m,i) = q(n,k,j,i) - (pmy_fluid->pmy_block->dx2f(j))*dqm;
=======
    for (int i=is; i<=ie; ++i){
      Real dql,dqr,dqc,q_jm1,q_j;
      if (n==NFLUID){
        q_jm1 = bcc(IB3,k,j-1,i);
        q_j   = bcc(IB3,k,j  ,i);
        dql = (bcc(IB3,k,j-1,i) - bcc(IB3,k,j-2,i))*dx2jm2i;
        dqc = (bcc(IB3,k,j  ,i) - bcc(IB3,k,j-1,i))*dx2jm1i;
        dqr = (bcc(IB3,k,j+1,i) - bcc(IB3,k,j  ,i))*dx2ji;
      } else if (n==(NFLUID+1)) {
        q_jm1 = bcc(IB1,k,j-1,i);
        q_j   = bcc(IB1,k,j  ,i);
        dql = (bcc(IB1,k,j-1,i) - bcc(IB1,k,j-2,i))*dx2jm2i;
        dqc = (bcc(IB1,k,j  ,i) - bcc(IB1,k,j-1,i))*dx2jm1i;
        dqr = (bcc(IB1,k,j+1,i) - bcc(IB1,k,j  ,i))*dx2ji;
      } else {
        q_jm1 = q(n,k,j-1,i);
        q_j   = q(n,k,j  ,i);
        dql = (q(n,k,j-1,i) - q(n,k,j-2,i))*dx2jm2i;
        dqc = (q(n,k,j  ,i) - q(n,k,j-1,i))*dx2jm1i;
        dqr = (q(n,k,j+1,i) - q(n,k,j  ,i))*dx2ji;
      }

      // Apply monotonicity constraints, compute ql_(i-1/2)
      Real dq2 = dql*dqc;
      Real dqm = (dq2 > 0.0) ? dqm = dq2/(dql + dqc) : 0.0;

      ql(n,i) = q_jm1 + (pmy_fluid->pmy_block->dx2f(j-1))*dqm;
      
      // Apply monotonicity constraints, compute qr_(i-1/2)
      dq2 = dqc*dqr;
      dqm = (dq2 > 0.0) ? dqm = dq2/(dqc + dqr) : 0.0;

      qr(n,i) = q_j   - (pmy_fluid->pmy_block->dx2f(j))*dqm;
    }
>>>>>>> remotes/origin/master
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn FluidIntegrator::ReconstructionFuncX3()
//  \brief 

void FluidIntegrator::PiecewiseLinearX3(const int k, const int j,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  int is = pmy_fluid->pmy_block->is, ie = pmy_fluid->pmy_block->ie;
  Real dx3km2i = 1.0/pmy_fluid->pmy_block->dx3v(k-2);
  Real dx3km1i = 1.0/pmy_fluid->pmy_block->dx3v(k-1);
  Real dx3ki   = 1.0/pmy_fluid->pmy_block->dx3v(k);

  for (int n=0; n<NWAVE; ++n) {
#pragma simd
<<<<<<< HEAD
  for (int i=is; i<=ie; ++i){

    Real dql = (q(n,k-1,j,i) - q(n,k-2,j,i))*dx3km2i;
    Real dqc = (q(n,k,j,i)   - q(n,k-1,j,i))*dx3km1i;
    Real dqr = (q(n,k+1,j,i) - q(n,k,j,i)  )*dx3ki;

// Apply monotonicity constraints to differences in primitive vars, compute ql_(i-1/2)

    Real dq2 = dql*dqc;
    Real dqm = 0.0;
    if (dq2 > 0.0) dqm = dq2/(dql + dqc);

    ql(m,i) = q(n,k-1,j,i) + (pmy_fluid->pmy_block->dx3f(k-1))*dqm;
=======
    for (int i=is; i<=ie; ++i){
      Real dql,dqr,dqc,q_km1,q_k;
      if (n==NFLUID){
        q_km1 = bcc(IB1,k-1,j,i);
        q_k   = bcc(IB1,k  ,j,i);
        dql = (bcc(IB1,k-1,j,i) - bcc(IB1,k-2,j,i))*dx3km2i;
        dqc = (bcc(IB1,k  ,j,i) - bcc(IB1,k-1,j,i))*dx3km1i;
        dqr = (bcc(IB1,k+1,j,i) - bcc(IB1,k  ,j,i))*dx3ki;
      } else if (n==(NFLUID+1)) {
        q_km1 = bcc(IB2,k-1,j,i);
        q_k   = bcc(IB2,k  ,j,i);
        dql = (bcc(IB2,k-1,j,i) - bcc(IB2,k-2,j,i))*dx3km2i;
        dqc = (bcc(IB2,k  ,j,i) - bcc(IB2,k-1,j,i))*dx3km1i;
        dqr = (bcc(IB2,k+1,j,i) - bcc(IB2,k  ,j,i))*dx3ki;
      } else {
        q_km1 = q(n,k-1,j,i);
        q_k   = q(n,k  ,j,i);
        dql = (q(n,k-1,j,i) - q(n,k-2,j,i))*dx3km2i;
        dqc = (q(n,k  ,j,i) - q(n,k-1,j,i))*dx3km1i;
        dqr = (q(n,k+1,j,i) - q(n,k  ,j,i))*dx3ki;
      }

      // Apply monotonicity constraints, compute ql_(i-1/2)
      Real dq2 = dql*dqc;
      Real dqm = (dq2 > 0.0) ? dqm = dq2/(dql + dqc) : 0.0;

      ql(n,i) = q_km1 + (pmy_fluid->pmy_block->dx3f(k-1))*dqm;
>>>>>>> remotes/origin/master
    
      // Apply monotonicity constraints, compute qr_(i-1/2)
      dq2 = dqc*dqr;
      dqm = (dq2 > 0.0) ? dqm = dq2/(dqc + dqr) : 0.0;

<<<<<<< HEAD
    dq2 = dqc*dqr;
    dqm = 0.0;
    if (dq2 > 0.0) dqm = dq2/(dqc + dqr);

    qr(m,i) = q(n,k,j,i) - (pmy_fluid->pmy_block->dx3f(k))*dqm;
=======
      qr(n,i) = q_k   - (pmy_fluid->pmy_block->dx3f(k))*dqm;
    }
>>>>>>> remotes/origin/master
  }

  return;
}
