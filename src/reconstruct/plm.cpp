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

void Reconstruction::PiecewiseLinearX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc, 
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> dwl,dwr,dw2,dwm,wc,bx;
  int ncells1 = (iu-il+1) + 2*(NGHOST);
  dwl.NewAthenaArray(NWAVE,ncells1);
  dwr.NewAthenaArray(NWAVE,ncells1);
  dw2.NewAthenaArray(NWAVE,ncells1);
  dwm.NewAthenaArray(NWAVE,ncells1);
  wc.NewAthenaArray(NWAVE,ncells1);
  bx.NewAthenaArray(ncells1);

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
    // compute L/R slopes for each variable
    for (int n=0; n<(NHYDRO); ++n) {
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        dwl(n,i) = (w(n,k,j,i  ) - w(n,k,j,i-1));
        dwr(n,i) = (w(n,k,j,i+1) - w(n,k,j,i  ));
        wc(n,i) = w(n,k,j,i);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        dwl(IBY,i) = (bcc(IB2,k,j,i  ) - bcc(IB2,k,j,i-1));
        dwr(IBY,i) = (bcc(IB2,k,j,i+1) - bcc(IB2,k,j,i  ));
        wc(IBY,i) = bcc(IB2,k,j,i);
        bx(i) = bcc(IB1,k,j,i);
      }
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        dwl(IBZ,i) = (bcc(IB3,k,j,i  ) - bcc(IB3,k,j,i-1));
        dwr(IBZ,i) = (bcc(IB3,k,j,i+1) - bcc(IB3,k,j,i  ));
        wc(IBZ,i) = bcc(IB3,k,j,i);
      }
    }

    // Project slopes to characteristic variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,dwl);
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,dwr);
    }

    // Apply van Leer limiter for uniform grid
    if (pmb->block_size.x1rat == 1.0) {
      for (int n=0; n<(NWAVE); ++n) {
#pragma simd
        for (int i=il-1; i<=iu; ++i){
          dw2(n,i) = dwl(n,i)*dwr(n,i);
          dwm(n,i) = 2.0*dw2(n,i)/(dwl(n,i) + dwr(n,i));
        }
        for (int i=il-1; i<=iu; ++i){
          if(dw2(n,i) <= 0.0) dwm(n,i) = 0.0;
        }
      }

    // Apply Mignone limiter for non-uniform grid
    } else {
      for (int n=0; n<(NWAVE); ++n) {
#pragma simd
        for (int i=il-1; i<=iu; ++i){
          dw2(n,i) = dwl(n,i)*dwr(n,i);
          Real cf = pco->dx1v(i  )/(pco->x1f(i+1) - pco->x1v(i));
          Real cb = pco->dx1v(i-1)/(pco->x1v(i  ) - pco->x1f(i));
          dwm(n,i) = (dw2(n,i)*(cf*dwl(n,i) + cb*dwr(n,i))/
            (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2(n,i)*(cf + cb - 2.0)));
        }
        for (int i=il-1; i<=iu; ++i){
          if(dw2(n,i) <= 0.0) dwm(n,i) = 0.0;
        }
      }
    }

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      VectorDotRightEigenmatrix(pmb,IVX,il-1,iu,bx,wc,dwm);
    }

    // compute ql_(i+1/2) and qr_(i-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        wl(n,k,j,i+1) = wc(n,i) + ((pco->x1f(i+1)-pco->x1v(i))/pco->dx1f(i))*dwm(n,i);
        wr(n,k,j,i  ) = wc(n,i) - ((pco->x1v(i  )-pco->x1f(i))/pco->dx1f(i))*dwm(n,i);
      }
    }

  }}

  dwl.DeleteAthenaArray();
  dwr.DeleteAthenaArray();
  dw2.DeleteAthenaArray();
  dwm.DeleteAthenaArray();
  wc.DeleteAthenaArray();
  bx.DeleteAthenaArray();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief 

void Reconstruction::PiecewiseLinearX2(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> dwl,dwr,dw2,dwm,wc,bx;
  int ncells1 = (iu-il+1) + 2*(NGHOST);
  dwl.NewAthenaArray(NWAVE,ncells1);
  dwr.NewAthenaArray(NWAVE,ncells1);
  dw2.NewAthenaArray(NWAVE,ncells1);
  dwm.NewAthenaArray(NWAVE,ncells1);
  wc.NewAthenaArray(NWAVE,ncells1);
  bx.NewAthenaArray(ncells1);

  for (int k=kl; k<=ku; ++k){
  for (int j=jl-1; j<=ju; ++j){
    // compute L/R slopes for each variable
    for (int n=0; n<(NHYDRO); ++n) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        dwl(n,i) = (w(n,k,j  ,i) - w(n,k,j-1,i));
        dwr(n,i) = (w(n,k,j+1,i) - w(n,k,j  ,i));
        wc(n,i) = w(n,k,j,i);
      }
    }

    if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        dwl(IBY,i) = (bcc(IB3,k,j  ,i) - bcc(IB3,k,j-1,i));
        dwr(IBY,i) = (bcc(IB3,k,j+1,i) - bcc(IB3,k,j  ,i));
        wc(IBY,i) = bcc(IB3,k,j,i);
        bx(i) = bcc(IB2,k,j,i);
      }  
#pragma simd
      for (int i=il; i<=iu; ++i){
        dwl(IBZ,i) = (bcc(IB1,k,j  ,i) - bcc(IB1,k,j-1,i));
        dwr(IBZ,i) = (bcc(IB1,k,j+1,i) - bcc(IB1,k,j  ,i));
        wc(IBZ,i) = bcc(IB1,k,j,i);
      }
    }

    // Project slopes to characteristic variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,dwl);
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,dwr);
    }

    // Apply van Leer limiter for uniform grid
    if (pmb->block_size.x1rat == 1.0) {
      for (int n=0; n<(NWAVE); ++n) {
#pragma simd
        for (int i=il; i<=iu; ++i){
          dw2(n,i) = dwl(n,i)*dwr(n,i);
          dwm(n,i) = 2.0*dw2(n,i)/(dwl(n,i) + dwr(n,i));
        }
        for (int i=il; i<=iu; ++i){
          if(dw2(n,i) <= 0.0) dwm(n,i) = 0.0;
        }
      }

    // Apply Mignone limiter for non-uniform grid
    } else {
      for (int n=0; n<(NWAVE); ++n) {
#pragma simd
        for (int i=il-1; i<=iu; ++i){
          dw2(n,i) = dwl(n,i)*dwr(n,i);
          Real cf = pco->dx2v(j  )/(pco->x2f(j+1) - pco->x2v(j));
          Real cb = pco->dx2v(j-1)/(pco->x2v(j  ) - pco->x2f(j));
          dwm(n,i) = (dw2(n,i)*(cf*dwl(n,i) + cb*dwr(n,i))/
            (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2(n,i)*(cf + cb - 2.0)));
        }
        for (int i=il-1; i<=iu; ++i){
          if(dw2(n,i) <= 0.0) dwm(n,i) = 0.0;
        }
      }
    }

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      VectorDotRightEigenmatrix(pmb,IVY,il,iu,bx,wc,dwm);
    }

    // compute ql_(j+1/2) and qr_(j-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        wl(n,k,j+1,i) = wc(n,i) + ((pco->x2f(j+1)-pco->x2v(j))/pco->dx2f(j))*dwm(n,i);
        wr(n,k,j  ,i) = wc(n,i) - ((pco->x2v(j  )-pco->x2f(j))/pco->dx2f(j))*dwm(n,i);
      }
    }
  }}

  dwl.DeleteAthenaArray();
  dwr.DeleteAthenaArray();
  dw2.DeleteAthenaArray();
  dwm.DeleteAthenaArray();
  wc.DeleteAthenaArray();
  bx.DeleteAthenaArray();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief 

void Reconstruction::PiecewiseLinearX3(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> dwl,dwr,dw2,dwm,wc,bx;
  int ncells1 = (iu-il+1) + 2*(NGHOST);
  dwl.NewAthenaArray(NWAVE,ncells1);
  dwr.NewAthenaArray(NWAVE,ncells1);
  dw2.NewAthenaArray(NWAVE,ncells1);
  dwm.NewAthenaArray(NWAVE,ncells1);
  wc.NewAthenaArray(NWAVE,ncells1);
  bx.NewAthenaArray(ncells1);

  for (int k=kl-1; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
    // compute L/R slopes for each variable
    for (int n=0; n<(NHYDRO); ++n) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        dwl(n,i) = (w(n,k  ,j,i) - w(n,k-1,j,i));
        dwr(n,i) = (w(n,k+1,j,i) - w(n,k  ,j,i));
        wc(n,i) = w(n,k,j,i);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        dwl(IBY,i) = (bcc(IB1,k  ,j,i) - bcc(IB1,k-1,j,i));
        dwr(IBY,i) = (bcc(IB1,k+1,j,i) - bcc(IB1,k  ,j,i));
        wc(IBY,i) = bcc(IB1,k,j,i);
        bx(i) = bcc(IB3,k,j,i);
      }  
#pragma simd
      for (int i=il; i<=iu; ++i){
        dwl(IBZ,i) = (bcc(IB2,k  ,j,i) - bcc(IB2,k-1,j,i));
        dwr(IBZ,i) = (bcc(IB2,k+1,j,i) - bcc(IB2,k  ,j,i));
        wc(IBZ,i) = bcc(IB2,k,j,i);
      }
    }

    // Project slopes to characteristic variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,dwl);
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,dwr);
    }


    // Apply van Leer limiter for uniform grid
    if (pmb->block_size.x1rat == 1.0) {
      for (int n=0; n<(NWAVE); ++n) {
#pragma simd
        for (int i=il; i<=iu; ++i){
          dw2(n,i) = dwl(n,i)*dwr(n,i);
          dwm(n,i) = 2.0*dw2(n,i)/(dwl(n,i) + dwr(n,i));
        }
        for (int i=il; i<=iu; ++i){
          if(dw2(n,i) <= 0.0) dwm(n,i) = 0.0;
        }
      }

    // Apply Mignone limiter for non-uniform grid
    } else {
      for (int n=0; n<(NWAVE); ++n) {
#pragma simd
        for (int i=il-1; i<=iu; ++i){
          dw2(n,i) = dwl(n,i)*dwr(n,i);
          Real cf = pco->dx3v(k  )/(pco->x3f(k+1) - pco->x3v(k));
          Real cb = pco->dx3v(k-1)/(pco->x3v(k  ) - pco->x3f(k));
          dwm(n,i) = (dw2(n,i)*(cf*dwl(n,i) + cb*dwr(n,i))/
            (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2(n,i)*(cf + cb - 2.0)));
        }
        for (int i=il-1; i<=iu; ++i){
          if(dw2(n,i) <= 0.0) dwm(n,i) = 0.0;
        }
      }
    }

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      VectorDotRightEigenmatrix(pmb,IVZ,il,iu,bx,wc,dwm);
    }

    // compute ql_(k+1/2) and qr_(k-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        wl(n,k+1,j,i) = wc(n,i) + ((pco->x3f(k+1)-pco->x3v(k))/pco->dx3f(k))*dwm(n,i);
        wr(n,k  ,j,i) = wc(n,i) - ((pco->x3v(k  )-pco->x3f(k))/pco->dx3f(k))*dwm(n,i);
      }
    }
  }}

  dwl.DeleteAthenaArray();
  dwr.DeleteAthenaArray();
  dw2.DeleteAthenaArray();
  dwm.DeleteAthenaArray();
  wc.DeleteAthenaArray();
  bx.DeleteAthenaArray();

  return;
}
