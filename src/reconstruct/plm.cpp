//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm-uniform.cpp
//  \brief  piecewise linear reconstruction for both uniform and non-uniform meshes

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX1()
//  \brief

void Reconstruction::PiecewiseLinearX1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Coordinates *pco = pmy_block_->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &bx = scr01_i_, &wc = scr1_ni_, &dwl = scr2_ni_, &dwr = scr3_ni_,
                   &dwm = scr4_ni_;

  // compute L/R slopes for each variable
  for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      dwl(n,i) = (w(n,k,j,i  ) - w(n,k,j,i-1));
      dwr(n,i) = (w(n,k,j,i+1) - w(n,k,j,i  ));
      wc(n,i) = w(n,k,j,i);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
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
  if (characteristic_reconstruction) {
    LeftEigenmatrixDotVector(IVX,il,iu,bx,wc,dwl);
    LeftEigenmatrixDotVector(IVX,il,iu,bx,wc,dwr);
  }

  // Apply van Leer limiter for uniform grid
  if (uniform_limiter[X1DIR]) {
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dw2 = dwl(n,i)*dwr(n,i);
        dwm(n,i) = 2.0*dw2/(dwl(n,i) + dwr(n,i));
        if (dw2 <= 0.0) dwm(n,i) = 0.0;
      }
    }

    // Apply Mignone limiter for non-uniform grid
  } else {
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dw2 = dwl(n,i)*dwr(n,i);
        Real cf = pco->dx1v(i  )/(pco->x1f(i+1) - pco->x1v(i));
        Real cb = pco->dx1v(i-1)/(pco->x1v(i  ) - pco->x1f(i));
        dwm(n,i) = (dw2*(cf*dwl(n,i) + cb*dwr(n,i))/
                    (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2*(cf + cb - 2.0)));
        if (dw2 <= 0.0) dwm(n,i) = 0.0;
      }
    }
  }

  // Project limited slope back to primitive variables, if necessary
  if (characteristic_reconstruction) {
    RightEigenmatrixDotVector(IVX,il,iu,bx,wc,dwm);
  }

  // compute ql_(i+1/2) and qr_(i-1/2) using monotonized slopes
  for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      wl(n,i+1) = wc(n,i) + ((pco->x1f(i+1)-pco->x1v(i))/pco->dx1f(i))*dwm(n,i);
      wr(n,i  ) = wc(n,i) - ((pco->x1v(i  )-pco->x1f(i))/pco->dx1f(i))*dwm(n,i);
    }
  }

  if (characteristic_reconstruction) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      // Reapply EOS floors to both L/R reconstructed primitive states
      // TODO(felker): check if fused loop with NWAVE redundant application is slower
      pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i+1);
      pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX2()
//  \brief

void Reconstruction::PiecewiseLinearX2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Coordinates *pco = pmy_block_->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &bx = scr01_i_, &wc = scr1_ni_, &dwl = scr2_ni_,
                   &dwr = scr3_ni_, &dwm = scr4_ni_;

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
  if (characteristic_reconstruction) {
    LeftEigenmatrixDotVector(IVY,il,iu,bx,wc,dwl);
    LeftEigenmatrixDotVector(IVY,il,iu,bx,wc,dwr);
  }

  // Apply van Leer limiter for uniform grid
  if (uniform_limiter[X2DIR]) {
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dw2 = dwl(n,i)*dwr(n,i);
        dwm(n,i) = 2.0*dw2/(dwl(n,i) + dwr(n,i));
        if (dw2 <= 0.0) dwm(n,i) = 0.0;
      }
    }

    // Apply Mignone limiter for non-uniform grid
  } else {
    Real cf = pco->dx2v(j  )/(pco->x2f(j+1) - pco->x2v(j));
    Real cb = pco->dx2v(j-1)/(pco->x2v(j  ) - pco->x2f(j));
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dw2 = dwl(n,i)*dwr(n,i);
        dwm(n,i) = (dw2*(cf*dwl(n,i) + cb*dwr(n,i))/
                    (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2*(cf + cb - 2.0)));
        if (dw2 <= 0.0) dwm(n,i) = 0.0;
      }
    }
  }

  // Project limited slope back to primitive variables, if necessary
  if (characteristic_reconstruction) {
    RightEigenmatrixDotVector(IVY,il,iu,bx,wc,dwm);
  }

  // compute ql_(j+1/2) and qr_(j-1/2) using monotonized slopes
  Real dxp = (pco->x2f(j+1)-pco->x2v(j))/pco->dx2f(j);
  Real dxm = (pco->x2v(j  )-pco->x2f(j))/pco->dx2f(j);
  for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = wc(n,i) + dxp*dwm(n,i);
      wr(n,i) = wc(n,i) - dxm*dwm(n,i);
    }
  }

  if (characteristic_reconstruction) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      // Reapply EOS floors to both L/R reconstructed primitive states
      pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i);
      pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseLinearX3()
//  \brief

void Reconstruction::PiecewiseLinearX3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Coordinates *pco = pmy_block_->pcoord;
  // set work arrays to shallow copies of scratch arrays
  AthenaArray<Real> &bx = scr01_i_, &wc = scr1_ni_, &dwl = scr2_ni_, &dwr = scr3_ni_,
                   &dwm = scr4_ni_;

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
  if (characteristic_reconstruction) {
    LeftEigenmatrixDotVector(IVZ,il,iu,bx,wc,dwl);
    LeftEigenmatrixDotVector(IVZ,il,iu,bx,wc,dwr);
  }


  // Apply van Leer limiter for uniform grid
  if (uniform_limiter[X3DIR]) {
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dw2 = dwl(n,i)*dwr(n,i);
        dwm(n,i) = 2.0*dw2/(dwl(n,i) + dwr(n,i));
        if (dw2 <= 0.0) dwm(n,i) = 0.0;
      }
    }

    // Apply Mignone limiter for non-uniform grid
  } else {
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real dw2 = dwl(n,i)*dwr(n,i);
        Real cf = pco->dx3v(k  )/(pco->x3f(k+1) - pco->x3v(k));
        Real cb = pco->dx3v(k-1)/(pco->x3v(k  ) - pco->x3f(k));
        dwm(n,i) = (dw2*(cf*dwl(n,i) + cb*dwr(n,i))/
                    (SQR(dwl(n,i)) + SQR(dwr(n,i)) + dw2*(cf + cb - 2.0)));
        if (dw2 <= 0.0) dwm(n,i) = 0.0;
      }
    }
  }

  // Project limited slope back to primitive variables, if necessary
  if (characteristic_reconstruction) {
    RightEigenmatrixDotVector(IVZ,il,iu,bx,wc,dwm);
  }

  // compute ql_(k+1/2) and qr_(k-1/2) using monotonized slopes
  Real dxp = (pco->x3f(k+1)-pco->x3v(k))/pco->dx3f(k);
  Real dxm = (pco->x3v(k  )-pco->x3f(k))/pco->dx3f(k);
  for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = wc(n,i) + dxp*dwm(n,i);
      wr(n,i) = wc(n,i) - dxm*dwm(n,i);
    }
  }

  if (characteristic_reconstruction) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      // Reapply EOS floors to both L/R reconstructed primitive states
      pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i);
      pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
    }
  }
  return;
}
