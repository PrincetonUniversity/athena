//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm.cpp
//  \brief  piecewise linear reconstruction for a non-uniform mesh

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

void Reconstruction::PiecewiseLinearX1(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  Real dql,dqr,dqc,q_im1,q_i;

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      Real& dx_im2 = pco->dx1v(i-2);
      Real& dx_im1 = pco->dx1v(i-1);
      Real& dx_i   = pco->dx1v(i);

      q_im1 = q(nin,k,j,i-1);
      q_i   = q(nin,k,j,i  );
      dql = (q(nin,k,j,i-1) - q(nin,k,j,i-2))/dx_im2;
      dqc = (q(nin,k,j,i  ) - q(nin,k,j,i-1))/dx_im1;
      dqr = (q(nin,k,j,i+1) - q(nin,k,j,i  ))/dx_i;

      // compute ql_(i-1/2) using Mignone 2014's modified van-Leer limiter
      Real dq2 = dql*dqc;
      ql(nout,k,j,i) = q_im1;
      if(dq2>0.0) {
        Real dxfr=pco->x1f(i)-pco->x1v(i-1);
        Real cf=dx_im1/dxfr;
        Real cb=dx_im2/(pco->x1v(i-1)-pco->x1f(i-1));
        ql(nout,k,j,i) += dxfr*dq2*(cf*dql+cb*dqc)/(dql*dql+(cf+cb-2.0)*dq2+dqc*dqc);
      }

      // compute qr_(i-1/2) using Mignone 2014's modified van-Leer limiter
      dq2 = dqc*dqr;
      qr(nout,k,j,i) = q_i;
      if(dq2>0.0) {
        Real dxfl=pco->x1v(i)-pco->x1f(i);
        Real cf=dx_i/(pco->x1f(i+1)-pco->x1v(i));
        Real cb=dx_im1/dxfl;
        qr(nout,k,j,i) -= dxfl*dq2*(cf*dqc+cb*dqr)/(dqc*dqc+(cf+cb-2.0)*dq2+dqr*dqr);
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief

void Reconstruction::PiecewiseLinearX2(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  Real dql,dqr,dqc,q_jm1,q_j;

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
    Real dx2_jm2 = pco->dx2v(j-2);
    Real dx2_jm1 = pco->dx2v(j-1);
    Real dx2_j   = pco->dx2v(j);
    Real dxfr=pco->x2f(j)-pco->x2v(j-1);
    Real dxfl=pco->x2v(j)-pco->x2f(j);
    Real dxfrp=pco->x2f(j+1)-pco->x2v(j);
    Real dxflm=pco->x2v(j-1)-pco->x2f(j-1);
    Real cfm=dx2_jm1/dxfr;
    Real cbm=dx2_jm2/dxflm;
    Real cfp=dx2_j/dxfrp;
    Real cbp=dx2_jm1/dxfl;
#pragma simd
    for (int i=il; i<=iu; ++i){
      q_jm1 = q(nin,k,j-1,i);
      q_j   = q(nin,k,j  ,i);
      dql = (q(nin,k,j-1,i) - q(nin,k,j-2,i))/dx2_jm2;
      dqc = (q(nin,k,j  ,i) - q(nin,k,j-1,i))/dx2_jm1;
      dqr = (q(nin,k,j+1,i) - q(nin,k,j  ,i))/dx2_j;

      // compute ql_(j-1/2) using Mignone 2014's modified van-Leer limiter
      Real dq2 = dql*dqc;
      ql(nout,k,j,i) = q_jm1;
      if(dq2>0.0) {
        ql(nout,k,j,i) += dxfr*dq2*(cfm*dql+cbm*dqc)/(dql*dql+(cfm+cbm-2.0)*dq2+dqc*dqc);
      }

      // compute qr_(j-1/2) using Mignone 2014's modified van-Leer limiter
      dq2 = dqc*dqr;
      qr(nout,k,j,i) = q_j;
      if(dq2>0.0) {
        qr(nout,k,j,i) -= dxfl*dq2*(cfp*dqc+cbp*dqr)/(dqc*dqc+(cfp+cbp-2.0)*dq2+dqr*dqr);
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief

void Reconstruction::PiecewiseLinearX3(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
//  Coordinates *pco = pmy_block_->pcoord;
  Real dql,dqr,dqc,q_km1,q_k;

  for (int k=kl; k<=ku; ++k){
    Real dx3km2i = 1.0/pco->dx3v(k-2);
    Real dx3km1i = 1.0/pco->dx3v(k-1);
    Real dx3ki   = 1.0/pco->dx3v(k);
    Real dxfr=pco->x3f(k)-pco->x3v(k-1);
    Real dxfl=pco->x3v(k)-pco->x3f(k);
    Real cfm=pco->dx3v(k-1)/dxfr;
    Real cbm=pco->dx3v(k-2)/(pco->x3v(k-1)-pco->x3f(k-1));
    Real cfp=pco->dx3v(k)/(pco->x3f(k+1)-pco->x3v(k));
    Real cbp=pco->dx3v(k-1)/dxfl;
    for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      q_km1 = q(nin,k-1,j,i);
      q_k   = q(nin,k  ,j,i);
      dql = (q(nin,k-1,j,i) - q(nin,k-2,j,i))*dx3km2i;
      dqc = (q(nin,k  ,j,i) - q(nin,k-1,j,i))*dx3km1i;
      dqr = (q(nin,k+1,j,i) - q(nin,k  ,j,i))*dx3ki;

      // compute ql_(k-1/2) using Mignone 2014's modified van-Leer limiter
      Real dq2 = dql*dqc;
      ql(nout,k,j,i) = q_km1;
      if(dq2>0.0) {
        ql(nout,k,j,i) += dxfr*dq2*(cfm*dql+cbm*dqc)/(dql*dql+(cfm+cbm-2.0)*dq2+dqc*dqc);
      }

      // compute qr_(k-1/2) using Mignone 2014's modified van-Leer limiter
      dq2 = dqc*dqr;
      qr(nout,k,j,i) = q_k;
      if(dq2>0.0) {
        qr(nout,k,j,i) -= dxfl*dq2*(cfp*dqc+cbp*dqr)/(dqc*dqc+(cfp+cbp-2.0)*dq2+dqr*dqr);
      }
    }}
  }

  return;
}
