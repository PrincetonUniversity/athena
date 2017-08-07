//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm-uniform.cpp
//  \brief  piecewise linear reconstruction for a uniform mesh

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX1()
//  \brief 

void Reconstruction::PiecewiseLinearUniformX1(Coordinates *pco,const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  Real dql,dqr,dqc,q_im1,q_i;

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      q_im1 = q(nin,k,j,i-1);
      q_i   = q(nin,k,j,i  );
      dql = (q(nin,k,j,i-1) - q(nin,k,j,i-2));
      dqc = (q(nin,k,j,i  ) - q(nin,k,j,i-1));
      dqr = (q(nin,k,j,i+1) - q(nin,k,j,i  ));

      // compute ql_(i-1/2) using van-Leer limiter
      Real dq2 = dql*dqc;
      ql(nout,k,j,i) = q_im1;
      if(dq2>0.0) {
        ql(nout,k,j,i) += dq2/(dql+dqc);
      }

      // compute qr_(i-1/2) using van-Leer limiter
      dq2 = dqc*dqr;
      qr(nout,k,j,i) = q_i;
      if(dq2>0.0) {
        qr(nout,k,j,i) -= dq2/(dqr+dqc); 
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief 

void Reconstruction::PiecewiseLinearUniformX2(Coordinates *pco,const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  Real dql,dqr,dqc,q_jm1,q_j;

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      q_jm1 = q(nin,k,j-1,i);
      q_j   = q(nin,k,j  ,i);
      dql = (q(nin,k,j-1,i) - q(nin,k,j-2,i));
      dqc = (q(nin,k,j  ,i) - q(nin,k,j-1,i));
      dqr = (q(nin,k,j+1,i) - q(nin,k,j  ,i));

      // compute ql_(j-1/2) using van-Leer limiter
      Real dq2 = dql*dqc;
      ql(nout,k,j,i) = q_jm1;
      if(dq2>0.0) {
        ql(nout,k,j,i) += dq2/(dql+dqc);
      }
       
      // compute qr_(j-1/2) using van-Leer limiter
      dq2 = dqc*dqr;
      qr(nout,k,j,i) = q_j;
      if(dq2>0.0) {
        qr(nout,k,j,i) -= dq2/(dqc+dqr); 
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief 

void Reconstruction::PiecewiseLinearUniformX3(Coordinates *pco,const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  Real dql,dqr,dqc,q_km1,q_k;

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      q_km1 = q(nin,k-1,j,i);
      q_k   = q(nin,k  ,j,i);
      dql = (q(nin,k-1,j,i) - q(nin,k-2,j,i));
      dqc = (q(nin,k  ,j,i) - q(nin,k-1,j,i));
      dqr = (q(nin,k+1,j,i) - q(nin,k  ,j,i));

      // compute ql_(k-1/2) using van-Leer limiter
      Real dq2 = dql*dqc;
      ql(nout,k,j,i) = q_km1;
      if(dq2>0.0) {
        ql(nout,k,j,i) += dq2/(dql+dqc);
      }
      
      // compute qr_(k-1/2) using van-Leer limiter
      dq2 = dqc*dqr;
      qr(nout,k,j,i) = q_k;
      if(dq2>0.0) {
        qr(nout,k,j,i) -= dq2/(dqc+dqr);
      }
    }
  }}

  return;
}
