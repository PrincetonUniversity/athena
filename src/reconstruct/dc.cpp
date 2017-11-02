//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc.cpp
//  \brief piecewise constant (donor cell) reconstruction

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX1()
//  \brief 

void Reconstruction::DonorCellX1(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(nout,k,j,i) = q(nin,k,j,i-1);
      qr(nout,k,j,i) = q(nin,k,j,i  );
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX2()
//  \brief 

void Reconstruction::DonorCellX2(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(nout,k,j,i) = q(nin,k,j-1,i);
      qr(nout,k,j,i) = q(nin,k,j  ,i);
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX3()
//  \brief 

void Reconstruction::DonorCellX3(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(nout,k,j,i) = q(nin,k-1,j,i);
      qr(nout,k,j,i) = q(nin,k  ,j,i);
    }
  }}

  return;
}
