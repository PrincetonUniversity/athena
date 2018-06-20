//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm-uniform.cpp
//  \brief  piecewise linear reconstruction for both uniform and non-uniform meshes

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"


//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX1()
//  \brief

void Reconstruction::PiecewiseLinearX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Coordinates *pco = pmb->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> bx,dw2,wc,dwl,dwr,dwm;
  bx.InitWithShallowCopy(pmb->precon->scr01_i_);
  dw2.InitWithShallowCopy(pmb->precon->scr02_i_);
  wc.InitWithShallowCopy(pmb->precon->scr1_ni_);
  dwl.InitWithShallowCopy(pmb->precon->scr2_ni_);
  dwr.InitWithShallowCopy(pmb->precon->scr3_ni_);
  dwm.InitWithShallowCopy(pmb->precon->scr4_ni_);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    // compute L/R slopes for each variable
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        dwl(n,i) = (w(n,k,j,i  ) - w(n,k,j,i-1));
        dwr(n,i) = (w(n,k,j,i+1) - w(n,k,j,i  ));
        wc(n,i) = w(n,k,j,i);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        bx(i) = bcc(IB1,k,j,i);

        dwl(IBY,i) = (bcc(IB2,k,j,i  ) - bcc(IB2,k,j,i-1));
        dwr(IBY,i) = (bcc(IB2,k,j,i+1) - bcc(IB2,k,j,i  ));
        wc(IBY,i) = bcc(IB2,k,j,i);

        dwl(IBZ,i) = (bcc(IB3,k,j,i  ) - bcc(IB3,k,j,i-1));
        dwr(IBZ,i) = (bcc(IB3,k,j,i+1) - bcc(IB3,k,j,i  ));
        wc(IBZ,i) = bcc(IB3,k,j,i);
      }
    }

    // Project slopes to characteristic variables, if necessary
    // Note order of characteristic fields in output vect corresponds to (IVX,IVY,IVZ)
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,dwl);
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,dwr);
    }

    // Apply van Leer limiter for uniform grid
    if (pmb->precon->uniform_limiter[X1DIR]) {
      for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il-1; i<=iu; ++i) {
          dw2(i) = dwl(n,i)*dwr(n,i);
          dwm(n,i) = 2.0*dw2(i)/(dwl(n,i) + dwr(n,i));
          if (dw2(i) <= 0.0) dwm(n,i) = 0.0;
        }
      }

    // Apply Mignone limiter for non-uniform grid
    } else {
      for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il-1; i<=iu; ++i) {
          dw2(i) = dwl(n,i)*dwr(n,i);
          Real cf = pco->dx1v(i  )/(pco->x1f(i+1) - pco->x1v(i));
          Real cb = pco->dx1v(i-1)/(pco->x1v(i  ) - pco->x1f(i));
          dwm(n,i) = (dw2(i)*(cf*dwl(n,i) + cb*dwr(n,i))/
            (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2(i)*(cf + cb - 2.0)));
          if (dw2(i) <= 0.0) dwm(n,i) = 0.0;
        }
      }
    }

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      RightEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,dwm);
    }

    // compute ql_(i+1/2) and qr_(i-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il-1; i<=iu; ++i) {
        wl(n,k,j,i+1) = wc(n,i) + ((pco->x1f(i+1)-pco->x1v(i))/pco->dx1f(i))*dwm(n,i);
        wr(n,k,j,i  ) = wc(n,i) - ((pco->x1v(i  )-pco->x1f(i))/pco->dx1f(i))*dwm(n,i);
        if (pmb->precon->characteristic_reconstruction) {
          // Reapply EOS floors to both L/R reconstructed primitive states
          pmb->peos->ApplyPrimitiveFloors(wl, k, j, i+1);
          pmb->peos->ApplyPrimitiveFloors(wr, k, j, i);
        }
      }
    }

  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX2()
//  \brief

void Reconstruction::PiecewiseLinearX2(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Coordinates *pco = pmb->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> bx,dw2,wc,dwl,dwr,dwm;
  bx.InitWithShallowCopy(pmb->precon->scr01_i_);
  dw2.InitWithShallowCopy(pmb->precon->scr02_i_);
  wc.InitWithShallowCopy(pmb->precon->scr1_ni_);
  dwl.InitWithShallowCopy(pmb->precon->scr2_ni_);
  dwr.InitWithShallowCopy(pmb->precon->scr3_ni_);
  dwm.InitWithShallowCopy(pmb->precon->scr4_ni_);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl-1; j<=ju; ++j) {
    // compute L/R slopes for each variable
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dwl(n,i) = (w(n,k,j  ,i) - w(n,k,j-1,i));
        dwr(n,i) = (w(n,k,j+1,i) - w(n,k,j  ,i));
        wc(n,i) = w(n,k,j,i);
      }
    }

    if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        bx(i) = bcc(IB2,k,j,i);

        dwl(IBY,i) = (bcc(IB3,k,j  ,i) - bcc(IB3,k,j-1,i));
        dwr(IBY,i) = (bcc(IB3,k,j+1,i) - bcc(IB3,k,j  ,i));
        wc(IBY,i) = bcc(IB3,k,j,i);

        dwl(IBZ,i) = (bcc(IB1,k,j  ,i) - bcc(IB1,k,j-1,i));
        dwr(IBZ,i) = (bcc(IB1,k,j+1,i) - bcc(IB1,k,j  ,i));
        wc(IBZ,i) = bcc(IB1,k,j,i);
      }
    }

    // Project slopes to characteristic variables, if necessary
    // Note order of characteristic fields in output vect corresponds to (IVY,IVZ,IVX)
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,dwl);
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,dwr);
    }

    // Apply van Leer limiter for uniform grid
    if (pmb->precon->uniform_limiter[X2DIR]) {
      for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          dw2(i) = dwl(n,i)*dwr(n,i);
          dwm(n,i) = 2.0*dw2(i)/(dwl(n,i) + dwr(n,i));
          if (dw2(i) <= 0.0) dwm(n,i) = 0.0;
        }
      }

    // Apply Mignone limiter for non-uniform grid
    } else {
      for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          dw2(i) = dwl(n,i)*dwr(n,i);
          Real cf = pco->dx2v(j  )/(pco->x2f(j+1) - pco->x2v(j));
          Real cb = pco->dx2v(j-1)/(pco->x2v(j  ) - pco->x2f(j));
          dwm(n,i) = (dw2(i)*(cf*dwl(n,i) + cb*dwr(n,i))/
            (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2(i)*(cf + cb - 2.0)));
          if (dw2(i) <= 0.0) dwm(n,i) = 0.0;
        }
      }
    }

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      RightEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,dwm);
    }

    // compute ql_(j+1/2) and qr_(j-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        wl(n,k,j+1,i) = wc(n,i) + ((pco->x2f(j+1)-pco->x2v(j))/pco->dx2f(j))*dwm(n,i);
        wr(n,k,j  ,i) = wc(n,i) - ((pco->x2v(j  )-pco->x2f(j))/pco->dx2f(j))*dwm(n,i);
        if (pmb->precon->characteristic_reconstruction) {
          // Reapply EOS floors to both L/R reconstructed primitive states
          pmb->peos->ApplyPrimitiveFloors(wl, k, j+1, i);
          pmb->peos->ApplyPrimitiveFloors(wr, k, j, i);
        }
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX3()
//  \brief

void Reconstruction::PiecewiseLinearX3(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Coordinates *pco = pmb->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> bx,dw2,wc,dwl,dwr,dwm;
  bx.InitWithShallowCopy(pmb->precon->scr01_i_);
  dw2.InitWithShallowCopy(pmb->precon->scr02_i_);
  wc.InitWithShallowCopy(pmb->precon->scr1_ni_);
  dwl.InitWithShallowCopy(pmb->precon->scr2_ni_);
  dwr.InitWithShallowCopy(pmb->precon->scr3_ni_);
  dwm.InitWithShallowCopy(pmb->precon->scr4_ni_);

  for (int k=kl-1; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    // compute L/R slopes for each variable
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dwl(n,i) = (w(n,k  ,j,i) - w(n,k-1,j,i));
        dwr(n,i) = (w(n,k+1,j,i) - w(n,k  ,j,i));
        wc(n,i) = w(n,k,j,i);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        bx(i) = bcc(IB3,k,j,i);

        dwl(IBY,i) = (bcc(IB1,k  ,j,i) - bcc(IB1,k-1,j,i));
        dwr(IBY,i) = (bcc(IB1,k+1,j,i) - bcc(IB1,k  ,j,i));
        wc(IBY,i) = bcc(IB1,k,j,i);

        dwl(IBZ,i) = (bcc(IB2,k  ,j,i) - bcc(IB2,k-1,j,i));
        dwr(IBZ,i) = (bcc(IB2,k+1,j,i) - bcc(IB2,k  ,j,i));
        wc(IBZ,i) = bcc(IB2,k,j,i);
      }
    }

    // Project slopes to characteristic variables, if necessary
    // Note order of characteristic fields in output vect corresponds to (IVZ,IVX,IVY)
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,dwl);
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,dwr);
    }


    // Apply van Leer limiter for uniform grid
    if (pmb->precon->uniform_limiter[X3DIR]) {
      for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          dw2(i) = dwl(n,i)*dwr(n,i);
          dwm(n,i) = 2.0*dw2(i)/(dwl(n,i) + dwr(n,i));
          if (dw2(i) <= 0.0) dwm(n,i) = 0.0;
        }
      }

    // Apply Mignone limiter for non-uniform grid
    } else {
      for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          dw2(i) = dwl(n,i)*dwr(n,i);
          Real cf = pco->dx3v(k  )/(pco->x3f(k+1) - pco->x3v(k));
          Real cb = pco->dx3v(k-1)/(pco->x3v(k  ) - pco->x3f(k));
          dwm(n,i) = (dw2(i)*(cf*dwl(n,i) + cb*dwr(n,i))/
            (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2(i)*(cf + cb - 2.0)));
          if (dw2(i) <= 0.0) dwm(n,i) = 0.0;
        }
      }
    }

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      RightEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,dwm);
    }

    // compute ql_(k+1/2) and qr_(k-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        wl(n,k+1,j,i) = wc(n,i) + ((pco->x3f(k+1)-pco->x3v(k))/pco->dx3f(k))*dwm(n,i);
        wr(n,k  ,j,i) = wc(n,i) - ((pco->x3v(k  )-pco->x3f(k))/pco->dx3f(k))*dwm(n,i);
        if (pmb->precon->characteristic_reconstruction) {
          // Reapply EOS floors to both L/R reconstructed primitive states
          pmb->peos->ApplyPrimitiveFloors(wl, k+1, j, i);
          pmb->peos->ApplyPrimitiveFloors(wr, k, j, i);
        }
      }
    }
  }}

  return;
}
