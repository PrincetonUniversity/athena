//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno.cpp
//! \brief WENO-Z type reconstriction functions

// C headers
#include <algorithm>      // max(), min()

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::WENOZ()
//! \brief WENO-Z: Borges et al. JCoPh, 227, 3191 (2008)
//! Reconstructs 5th-order polynomial in cell i to compute ql(i+1) and qr(i).
//! Works for any dimension by passing in the appropriate q_im2,...,q _ip2.
inline void Reconstruction::WENOZ(const Real &q_im2, const Real &q_im1, const Real &q_i,
            const Real &q_ip1, const Real &q_ip2, Real &ql_ip1, Real &qr_i) noexcept  {
  constexpr Real beta_coeff[2]{13.0/12.0, 0.25};

  Real beta[3];
  beta[0] = beta_coeff[0] * SQR(q_im2 +     q_i - 2.0*q_im1) +
            beta_coeff[1] * SQR(q_im2 + 3.0*q_i - 4.0*q_im1);

  beta[1] = beta_coeff[0] * SQR(q_im1 + q_ip1 - 2.0*q_i) +
            beta_coeff[1] * SQR(q_im1 - q_ip1);

  beta[2] = beta_coeff[0] * SQR(q_ip2 +      q_i - 2.0*q_ip1) +
            beta_coeff[1] * SQR(q_ip2 + 3.0* q_i - 4.0*q_ip1);

  // Rescale epsilon
  constexpr Real epsL = 1.0e-42;

  const Real tau_5 = std::abs(beta[0] - beta[2]);

  Real indicator[3];
  indicator[0] = tau_5 / (beta[0] + epsL);
  indicator[1] = tau_5 / (beta[1] + epsL);
  indicator[2] = tau_5 / (beta[2] + epsL);

  // compute qL_ip1
  // Factor of 1/6 in coefficients of f[] array applied to alpha_sum to reduce divisions
  Real f[3];
  f[0] = ( 2.0*q_im2 - 7.0*q_im1 + 11.0*q_i  );
  f[1] = (-1.0*q_im1 + 5.0*q_i   + 2.0 *q_ip1);
  f[2] = ( 2.0*q_i   + 5.0*q_ip1 -      q_ip2);

  Real alpha[3];
  alpha[0] = 0.1*(1.0 + SQR(indicator[0]));
  alpha[1] = 0.6*(1.0 + SQR(indicator[1]));
  alpha[2] = 0.3*(1.0 + SQR(indicator[2]));
  Real alpha_sum = 6.0*(alpha[0] + alpha[1] + alpha[2]);

  ql_ip1 = (f[0]*alpha[0] + f[1]*alpha[1] + f[2]*alpha[2])/alpha_sum;

  // compute qR_i
  // Factor of 1/6 in coefficients of f[] array applied to alpha_sum to reduce divisions
  f[0] = ( 2.0*q_ip2 - 7.0*q_ip1 + 11.0*q_i  );
  f[1] = (-1.0*q_ip1 + 5.0*q_i   + 2.0 *q_im1);
  f[2] = ( 2.0*q_i   + 5.0*q_im1 -      q_im2);

  alpha[0] = 0.1*(1.0 + SQR(indicator[2]));
  alpha[1] = 0.6*(1.0 + SQR(indicator[1]));
  alpha[2] = 0.3*(1.0 + SQR(indicator[0]));
  alpha_sum = 6.0*(alpha[0] + alpha[1] + alpha[2]);

  qr_i = (f[0]*alpha[0] + f[1]*alpha[1] + f[2]*alpha[2])/alpha_sum;

  return;
}


//----------------------------------------------------------------------------------------
//! \fn Reconstruction::WENOMZ()
//! \brief WENO-MZ: Wang et al., Int. J. for Numerical Methods in Fluids, 96, 1631 (2024)
//! Reconstructs 5th-order polynomial in cell i to compute ql(i+1) and qr(i).
//! Works for any dimension by passing in the appropriate q_im2,...,q _ip2.
inline void Reconstruction::WENOMZ(const Real &q_im2, const Real &q_im1, const Real &q_i,
            const Real &q_ip1, const Real &q_ip2, Real &ql_ip1, Real &qr_i) noexcept  {
  constexpr Real beta_coeff0 = 13.0/12.0, beta_coeff1 = 0.25, beta_coeff4 = 1.0/12.0;
  constexpr Real epsL = 1.0e-40;

  Real b4 = beta_coeff4 * SQR(q_im1 + q_ip1 - 2.0*q_i);
  Real b0 = beta_coeff0 * SQR(q_im2 +     q_i - 2.0*q_im1) +
            beta_coeff1 * SQR(q_im2 + 3.0*q_i - 4.0*q_im1);
  Real b1 = 13.0 * b4 + beta_coeff1 * SQR(q_im1 - q_ip1);
  Real b2 = beta_coeff0 * SQR(q_ip2 +      q_i - 2.0*q_ip1) +
            beta_coeff1 * SQR(q_ip2 + 3.0* q_i - 4.0*q_ip1);

  Real tau_5 = std::abs(b0 - b2);
  Real r = (std::abs(b2 - b1) + epsL) / (std::abs(b1 - b0) + epsL);
  Real t0 = 1.0 + r, t2 = 1.0 + 1.0 / r;
  constexpr Real t1 = 2.0;

  Real tmp = tau_5 / (std::max(b0, b2) + epsL);
  Real eta = SQR(SQR(tmp));
  Real tmpA = eta * tau_5;
  Real tmpB = (1.0 - eta) * tau_5;

  Real indA0 = tmpA / (b0 + epsL);
  Real indA1 = tmpA / (b1 + epsL);
  Real indA2 = tmpA / (b2 + epsL);
  Real indB0 = tmpB / (t0 * b4 + epsL);
  Real indB1 = tmpB / (t1 * b4 + epsL);
  Real indB2 = tmpB / (t2 * b4 + epsL);

  Real a1  = 0.6 * (1.0 + indA1 + indB1);
  Real a0 = 1.0 + indA0 + indB0;
  Real a2 = 1.0 + indA2 + indB2;
  Real al0 = 0.1 * a0;
  Real al2 = 0.3 * a2;
  Real ar0 = 0.1 * a2;
  Real ar2 = 0.3 * a0;

  Real fl0 = ( 2.0*q_im2 - 7.0*q_im1 + 11.0*q_i  );
  Real fl1 = (-1.0*q_im1 + 5.0*q_i   + 2.0 *q_ip1);
  Real fl2 = ( 2.0*q_i   + 5.0*q_ip1 -      q_ip2);
  ql_ip1 = (fl0*al0 + fl1*a1 + fl2*al2) / (6.0 * (al0 + a1 + al2));

  Real fr0 = ( 2.0*q_ip2 - 7.0*q_ip1 + 11.0*q_i  );
  Real fr1 = (-1.0*q_ip1 + 5.0*q_i   + 2.0 *q_im1);
  Real fr2 = ( 2.0*q_i   + 5.0*q_im1 -      q_im2);
  qr_i = (fr0*ar0 + fr1*a1 + fr2*ar2) / (6.0 * (ar0 + a1 + ar2));

  return;
}


//----------------------------------------------------------------------------------------
//! \fn WENOZX1
//! \brief Wrapper function for the WENO-Z reconstruction in x1-direction.
//! This function should be called over [is-1,ie+1] to get BOTH L/R states over [is,ie]
void Reconstruction::WENOZX1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Real dfloor = pmy_block_->peos->GetDensityFloor();
  Real pfloor = pmy_block_->peos->GetPressureFloor();

  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = w(n,k,j,i-2);
      Real qim1 = w(n,k,j,i-1);
      Real qi   = w(n,k,j,i  );
      Real qip1 = w(n,k,j,i+1);
      Real qip2 = w(n,k,j,i+2);
      WENOZ(qim2, qim1, qi, qip1, qip2, wl(n,i+1), wr(n,i));
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = bcc(IB2,k,j,i-2);
      Real qim1 = bcc(IB2,k,j,i-1);
      Real qi   = bcc(IB2,k,j,i  );
      Real qip1 = bcc(IB2,k,j,i+1);
      Real qip2 = bcc(IB2,k,j,i+2);
      WENOZ(qim2, qim1, qi, qip1, qip2, wl(IBY,i+1), wr(IBY,i));
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = bcc(IB3,k,j,i-2);
      Real qim1 = bcc(IB3,k,j,i-1);
      Real qi   = bcc(IB3,k,j,i  );
      Real qip1 = bcc(IB3,k,j,i+1);
      Real qip2 = bcc(IB3,k,j,i+2);
      WENOZ(qim2, qim1, qi, qip1, qip2, wl(IBZ,i+1), wr(IBZ,i));
    }
  }

#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i+1);
    pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOZX2
//! \brief Wrapper function for the WENO-Z reconstruction in x2-direction.
//! This function should be called over [js-1,je+1] to get BOTH L/R states over [js,je]
void Reconstruction::WENOZX2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {

  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = w(n,k,j-2,i);
      Real qjm1 = w(n,k,j-1,i);
      Real qj   = w(n,k,j  ,i);
      Real qjp1 = w(n,k,j+1,i);
      Real qjp2 = w(n,k,j+2,i);
      WENOZ(qjm2, qjm1, qj, qjp1, qjp2, wl(n,i), wr(n,i));
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = bcc(IB3,k,j-2,i);
      Real qjm1 = bcc(IB3,k,j-1,i);
      Real qj   = bcc(IB3,k,j  ,i);
      Real qjp1 = bcc(IB3,k,j+1,i);
      Real qjp2 = bcc(IB3,k,j+2,i);
      WENOZ(qjm2, qjm1, qj, qjp1, qjp2, wl(IBY,i), wr(IBY,i));
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = bcc(IB1,k,j-2,i);
      Real qjm1 = bcc(IB1,k,j-1,i);
      Real qj   = bcc(IB1,k,j  ,i);
      Real qjp1 = bcc(IB1,k,j+1,i);
      Real qjp2 = bcc(IB1,k,j+2,i);
      WENOZ(qjm2, qjm1, qj, qjp1, qjp2, wl(IBZ,i), wr(IBZ,i));
    }
  }

#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i);
    pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOZX3
//! \brief Wrapper function for the WENO-Z reconstruction in x3-direction.
//! This function should be called over [ks-1,ke+1] to get BOTH L/R states over [ks,ke]
void Reconstruction::WENOZX3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {

  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = w(n,k-2,j,i);
      Real qkm1 = w(n,k-1,j,i);
      Real qk   = w(n,k  ,j,i);
      Real qkp1 = w(n,k+1,j,i);
      Real qkp2 = w(n,k+2,j,i);
      WENOZ(qkm2, qkm1, qk, qkp1, qkp2, wl(n,i), wr(n,i));
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = bcc(IB1,k-2,j,i);
      Real qkm1 = bcc(IB1,k-1,j,i);
      Real qk   = bcc(IB1,k  ,j,i);
      Real qkp1 = bcc(IB1,k+1,j,i);
      Real qkp2 = bcc(IB1,k+2,j,i);
      WENOZ(qkm2, qkm1, qk, qkp1, qkp2, wl(IBY,i), wr(IBY,i));
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = bcc(IB2,k-2,j,i);
      Real qkm1 = bcc(IB2,k-1,j,i);
      Real qk   = bcc(IB2,k  ,j,i);
      Real qkp1 = bcc(IB2,k+1,j,i);
      Real qkp2 = bcc(IB2,k+2,j,i);
      WENOZ(qkm2, qkm1, qk, qkp1, qkp2, wl(IBZ,i), wr(IBZ,i));
    }
  }

#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i);
    pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn WENOZX1
//! \brief WENO-Z reconstruction in x1-direction for non-hydro variables
void Reconstruction::WENOZX1(const int k, const int j, const int il, const int iu,
       const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  for (int n=0; n<nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = q(n,k,j,i-2);
      Real qim1 = q(n,k,j,i-1);
      Real qi   = q(n,k,j,i  );
      Real qip1 = q(n,k,j,i+1);
      Real qip2 = q(n,k,j,i+2);
      WENOZ(qim2, qim1, qi, qip1, qip2, ql(n,i+1), qr(n,i));
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOZX2
//! \brief WENO-Z reconstruction in x2-direction for non-hydro variables
void Reconstruction::WENOZX2(const int k, const int j, const int il, const int iu,
       const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  for (int n=0; n<nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = q(n,k,j-2,i);
      Real qjm1 = q(n,k,j-1,i);
      Real qj   = q(n,k,j  ,i);
      Real qjp1 = q(n,k,j+1,i);
      Real qjp2 = q(n,k,j+2,i);
      WENOZ(qjm2, qjm1, qj, qjp1, qjp2, ql(n,i), qr(n,i));
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOZX3
//! \brief WENO-Z reconstruction in x3-direction for non-hydro variables
void Reconstruction::WENOZX3(const int k, const int j, const int il, const int iu,
       const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  for (int n=0; n<nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = q(n,k-2,j,i);
      Real qkm1 = q(n,k-1,j,i);
      Real qk   = q(n,k  ,j,i);
      Real qkp1 = q(n,k+1,j,i);
      Real qkp2 = q(n,k+2,j,i);
      WENOZ(qkm2, qkm1, qk, qkp1, qkp2, ql(n,i), qr(n,i));
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn WENOMZX1
//! \brief Wrapper function for the WENO-MZ reconstruction in x1-direction.
//! This function should be called over [is-1,ie+1] to get BOTH L/R states over [is,ie]
void Reconstruction::WENOMZX1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Real dfloor = pmy_block_->peos->GetDensityFloor();
  Real pfloor = pmy_block_->peos->GetPressureFloor();

  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = w(n,k,j,i-2);
      Real qim1 = w(n,k,j,i-1);
      Real qi   = w(n,k,j,i  );
      Real qip1 = w(n,k,j,i+1);
      Real qip2 = w(n,k,j,i+2);
      WENOMZ(qim2, qim1, qi, qip1, qip2, wl(n,i+1), wr(n,i));
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = bcc(IB2,k,j,i-2);
      Real qim1 = bcc(IB2,k,j,i-1);
      Real qi   = bcc(IB2,k,j,i  );
      Real qip1 = bcc(IB2,k,j,i+1);
      Real qip2 = bcc(IB2,k,j,i+2);
      WENOMZ(qim2, qim1, qi, qip1, qip2, wl(IBY,i+1), wr(IBY,i));
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = bcc(IB3,k,j,i-2);
      Real qim1 = bcc(IB3,k,j,i-1);
      Real qi   = bcc(IB3,k,j,i  );
      Real qip1 = bcc(IB3,k,j,i+1);
      Real qip2 = bcc(IB3,k,j,i+2);
      WENOMZ(qim2, qim1, qi, qip1, qip2, wl(IBZ,i+1), wr(IBZ,i));
    }
  }

#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i+1);
    pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOMZX2
//! \brief Wrapper function for the WENO-MZ reconstruction in x2-direction.
//! This function should be called over [js-1,je+1] to get BOTH L/R states over [js,je]
void Reconstruction::WENOMZX2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {

  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = w(n,k,j-2,i);
      Real qjm1 = w(n,k,j-1,i);
      Real qj   = w(n,k,j  ,i);
      Real qjp1 = w(n,k,j+1,i);
      Real qjp2 = w(n,k,j+2,i);
      WENOMZ(qjm2, qjm1, qj, qjp1, qjp2, wl(n,i), wr(n,i));
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = bcc(IB3,k,j-2,i);
      Real qjm1 = bcc(IB3,k,j-1,i);
      Real qj   = bcc(IB3,k,j  ,i);
      Real qjp1 = bcc(IB3,k,j+1,i);
      Real qjp2 = bcc(IB3,k,j+2,i);
      WENOMZ(qjm2, qjm1, qj, qjp1, qjp2, wl(IBY,i), wr(IBY,i));
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = bcc(IB1,k,j-2,i);
      Real qjm1 = bcc(IB1,k,j-1,i);
      Real qj   = bcc(IB1,k,j  ,i);
      Real qjp1 = bcc(IB1,k,j+1,i);
      Real qjp2 = bcc(IB1,k,j+2,i);
      WENOMZ(qjm2, qjm1, qj, qjp1, qjp2, wl(IBZ,i), wr(IBZ,i));
    }
  }

#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i);
    pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOMZX3
//! \brief Wrapper function for the WENO-MZ reconstruction in x3-direction.
//! This function should be called over [ks-1,ke+1] to get BOTH L/R states over [ks,ke]
void Reconstruction::WENOMZX3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr) {

  for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = w(n,k-2,j,i);
      Real qkm1 = w(n,k-1,j,i);
      Real qk   = w(n,k  ,j,i);
      Real qkp1 = w(n,k+1,j,i);
      Real qkp2 = w(n,k+2,j,i);
      WENOMZ(qkm2, qkm1, qk, qkp1, qkp2, wl(n,i), wr(n,i));
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = bcc(IB1,k-2,j,i);
      Real qkm1 = bcc(IB1,k-1,j,i);
      Real qk   = bcc(IB1,k  ,j,i);
      Real qkp1 = bcc(IB1,k+1,j,i);
      Real qkp2 = bcc(IB1,k+2,j,i);
      WENOMZ(qkm2, qkm1, qk, qkp1, qkp2, wl(IBY,i), wr(IBY,i));
    }
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = bcc(IB2,k-2,j,i);
      Real qkm1 = bcc(IB2,k-1,j,i);
      Real qk   = bcc(IB2,k  ,j,i);
      Real qkp1 = bcc(IB2,k+1,j,i);
      Real qkp2 = bcc(IB2,k+2,j,i);
      WENOMZ(qkm2, qkm1, qk, qkp1, qkp2, wl(IBZ,i), wr(IBZ,i));
    }
  }

#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    pmy_block_->peos->ApplyPrimitiveFloors(wl, k, j, i);
    pmy_block_->peos->ApplyPrimitiveFloors(wr, k, j, i);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn WENOMZX1
//! \brief WENO-MZ reconstruction in x1-direction for non-hydro variables
void Reconstruction::WENOMZX1(const int k, const int j, const int il, const int iu,
       const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  for (int n=0; n<nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qim2 = q(n,k,j,i-2);
      Real qim1 = q(n,k,j,i-1);
      Real qi   = q(n,k,j,i  );
      Real qip1 = q(n,k,j,i+1);
      Real qip2 = q(n,k,j,i+2);
      WENOMZ(qim2, qim1, qi, qip1, qip2, ql(n,i+1), qr(n,i));
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOMZX2
//! \brief WENO-MZ reconstruction in x2-direction for non-hydro variables
void Reconstruction::WENOMZX2(const int k, const int j, const int il, const int iu,
       const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  for (int n=0; n<nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qjm2 = q(n,k,j-2,i);
      Real qjm1 = q(n,k,j-1,i);
      Real qj   = q(n,k,j  ,i);
      Real qjp1 = q(n,k,j+1,i);
      Real qjp2 = q(n,k,j+2,i);
      WENOMZ(qjm2, qjm1, qj, qjp1, qjp2, ql(n,i), qr(n,i));
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn WENOMZX3
//! \brief WENO-MZ reconstruction in x3-direction for non-hydro variables
void Reconstruction::WENOMZX3(const int k, const int j, const int il, const int iu,
       const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;

  for (int n=0; n<nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real qkm2 = q(n,k-2,j,i);
      Real qkm1 = q(n,k-1,j,i);
      Real qk   = q(n,k  ,j,i);
      Real qkp1 = q(n,k+1,j,i);
      Real qkp2 = q(n,k+2,j,i);
      WENOMZ(qkm2, qkm1, qk, qkp1, qkp2, ql(n,i), qr(n,i));
    }
  }
  return;
}
