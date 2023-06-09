//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file characteristic.cpp
//! \brief Functions to transform vectors between primitive and characteristic variables

// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::LeftEigenmatrixDotVector()
//! \brief Computes inner-product of left-eigenmatrix of Roe's matrix A in the primitive
//! variables and an input vector.  This operation converts primitive to characteristic
//! variables.  The result is returned in the input vector, with the components of the
//! characteristic field stored such that vect(1,i) is in the direction of the sweep.
//!
//! The order of the components in the input vector should be:
//!    (IDN,IVX,IVY,IVZ,[IPR],[IBY,IBZ])
//! and these are permuted according to the direction specified by the input flag "ivx".
//!
//! REFERENCES:
//! - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//!   astrophysical MHD", ApJS, (2008), Appendix A. (Equation numbers refer to this paper)
//!
//! - Roe, Philip L., and Dinshaw S. Balsara. "Notes on the eigensystem of
//!   magnetohydrodynamics." SIAM Journal on Applied Mathematics 56.1 (1996): 57-67.
//!
//! - Brio, Moysey, and Cheng Chin Wu. "An upwind differencing scheme for the equations of
//!   ideal magnetohydrodynamics." Journal of computational physics 75.2 (1988): 400-422.


void Reconstruction::LeftEigenmatrixDotVector(
    const int ivx, const int il, const int iu,
    const AthenaArray<Real> &b1, const AthenaArray<Real> &w, AthenaArray<Real> &vect) {
  // permute components of input primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  if (MAGNETIC_FIELDS_ENABLED) {
    // Adiabatic MHD ---------------------------------------------------------------------
    if (NON_BAROTROPIC_EOS) {
      Real gamma = pmy_block_->peos->GetGamma();
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real id = 1.0/w(IDN,i);
        Real sqrtd = std::sqrt(w(IDN,i));
        Real isqrtd = 1.0/sqrtd;

        Real btsq = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real bxsq = b1(i)*b1(i);
        Real gamp = gamma*w(IPR,i);

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = bxsq + btsq - gamp;
        Real cf2_cs2 = std::sqrt(tdif*tdif + 4.0*gamp*btsq);

        Real cfsq = 0.5*(bxsq + btsq + gamp + cf2_cs2);
        Real cssq = gamp*bxsq/cfsq;

        cfsq *= id;
        Real cf = std::sqrt(cfsq);

        cssq *= id;
        Real cs = std::sqrt(cssq);

        Real asq = gamp*id;
        Real a = std::sqrt(asq);

        // Compute beta(s) (eq A17)
        Real bt  = std::sqrt(btsq);
        // edge case when transverse fields disappear: they can be set arbitrarily, except
        // cannot both be 0. this preserves orthonormality of the eigenvectors.
        // See BW88 eq 45, Roe96 pg 60
        Real bet2 = 1.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        // handle special cases in which bet2=be3=0
        if ((cfsq - cssq) <= 0.0) { // Roe96 case V (the triple umbilic)
          // degenerate case mentioned before A18
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (asq - cssq) <= 0.0) { // Roe96 case IV
          // low beta; slow waves are regular acoustic waves
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - asq) <= 0.0) { // Roe96 case III
          // high beta; fast waves are regular acoustic waves
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = std::sqrt((asq - cssq)/(cfsq - cssq));
          alpha_s = std::sqrt((cfsq - asq)/(cfsq - cssq));
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real nf = 0.5/asq;
        Real qf = nf*cf*alpha_f*s;
        Real qs = nf*cs*alpha_s*s;
        Real af_prime = 0.5*alpha_f/(a*sqrtd);
        Real as_prime = 0.5*alpha_s/(a*sqrtd);

        // Multiply row of L-eigenmatrix with vector using matrix elements from eq. A18
        Real v_0 = nf*alpha_f*(vect(IPR,i)*id - cf*vect(ivx,i)) +
                   qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                   as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_1 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd + vect(ivz,i)) -
                        bet3*(vect(IBY,i)*s*isqrtd + vect(ivy,i)));
        Real v_2 = nf*alpha_s*(vect(IPR,i)*id - cs*vect(ivx,i)) -
                   qf*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) -
                   af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_3 = vect(IDN,i) - vect(IPR,i)/asq;
        Real v_4 = nf*alpha_s*(vect(IPR,i)*id + cs*vect(ivx,i)) +
                   qf*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) -
                   af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_5 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd - vect(ivz,i)) -
                        bet3*(vect(IBY,i)*s*isqrtd - vect(ivy,i)));
        Real v_6 = nf*alpha_f*(vect(IPR,i)*id + cf*vect(ivx,i)) -
                   qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                   as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));

        vect(0,i) = v_0;
        vect(1,i) = v_1;
        vect(2,i) = v_2;
        vect(3,i) = v_3;
        vect(4,i) = v_4;
        vect(5,i) = v_5;
        vect(6,i) = v_6;
      }

      // Isothermal MHD ------------------------------------------------------------------
    } else {
      Real iso_cs = pmy_block_->peos->GetIsoSoundSpeed();
      Real iso_cs2 = SQR(iso_cs);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real id = 1.0/w(IDN,i);
        Real sqrtd = std::sqrt(w(IDN,i));
        Real isqrtd = 1.0/sqrtd;

        Real btsq = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real bxsq = b1(i)*b1(i);
        Real gamp = iso_cs2*w(IDN,i);

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = bxsq + btsq - gamp;
        Real cf2_cs2 = std::sqrt(tdif*tdif + 4.0*gamp*btsq);

        Real cfsq = 0.5*(bxsq + btsq + gamp + cf2_cs2);
        Real cssq = gamp*bxsq/cfsq;

        cfsq *= id;
        Real cf = std::sqrt(cfsq);

        cssq *= id;
        Real cs = std::sqrt(cssq);

        // Compute beta(s) (eq A17)
        Real bt  = std::sqrt(btsq);
        Real bet2 = 1.0;  // see above note
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if ((cfsq-cssq) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (iso_cs2 - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - iso_cs2) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = std::sqrt((iso_cs2 - cssq)/(cfsq - cssq));
          alpha_s = std::sqrt((cfsq - iso_cs2)/(cfsq - cssq));
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real qf = 0.5*(cf*alpha_f*s)/iso_cs2;
        Real qs = 0.5*(cs*alpha_s*s)/iso_cs2;
        Real af_prime = 0.5*alpha_f/(iso_cs*sqrtd);
        Real as_prime = 0.5*alpha_s/(iso_cs*sqrtd);

        // Multiply row of L-eigenmatrix with vector using matrix elements from eq. A22
        Real v_0 = 0.5*alpha_f*(vect(IDN,i)*id - cf*vect(ivx,i)/iso_cs2) +
                   qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                   as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_1 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd + vect(ivz,i)) -
                        bet3*(vect(IBY,i)*s*isqrtd + vect(ivy,i)));
        Real v_2 = 0.5*alpha_s*(vect(IDN,i)*id - cs*vect(ivx,i)/iso_cs2) -
                   qf*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) -
                   af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_3 = 0.5*alpha_s*(vect(IDN,i)*id + cs*vect(ivx,i)/iso_cs2) +
                   qf*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) -
                   af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_4 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd - vect(ivz,i)) -
                        bet3*(vect(IBY,i)*s*isqrtd - vect(ivy,i)));
        Real v_5 = 0.5*alpha_f*(vect(IDN,i)*id + cf*vect(ivx,i)/iso_cs2) -
                   qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                   as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));

        vect(0,i) = v_0;
        vect(1,i) = v_1;
        vect(2,i) = v_2;
        vect(3,i) = v_3;
        vect(4,i) = v_4;
        vect(5,i) = v_5;
      }
    }

  } else {
    // Adiabatic hydrodynamics -----------------------------------------------------------
    if (NON_BAROTROPIC_EOS) {
      Real gamma = pmy_block_->peos->GetGamma();
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real asq = gamma*w(IPR,i)/w(IDN,i);
        Real a   = std::sqrt(asq);

        // Multiply row of L-eigenmatrix with vector using matrix elements from eq. A4
        Real v_0 = 0.5*(vect(IPR,i)/asq - w(IDN,i)*vect(ivx,i)/a);
        Real v_1 = vect(IDN,i) - vect(IPR,i)/asq;
        Real v_2 = vect(ivy,i);
        Real v_3 = vect(ivz,i);
        Real v_4 = 0.5*(vect(IPR,i)/asq + w(IDN,i)*vect(ivx,i)/a);

        vect(0,i) = v_0;
        vect(1,i) = v_1;
        vect(2,i) = v_2;
        vect(3,i) = v_3;
        vect(4,i) = v_4;
      }

      // Isothermal hydrodynamics --------------------------------------------------------
    } else {
      Real iso_cs = pmy_block_->peos->GetIsoSoundSpeed();
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        // Multiply row of L-eigenmatrix with vector using matrix elements from eq. A7
        Real v_0 = 0.5*(vect(IDN,i) - w(IDN,i)*vect(ivx,i)/iso_cs);
        Real v_1 = vect(ivy,i);
        Real v_2 = vect(ivz,i);
        Real v_3 = 0.5*(vect(IDN,i) + w(IDN,i)*vect(ivx,i)/iso_cs);

        vect(0,i) = v_0;
        vect(1,i) = v_1;
        vect(2,i) = v_2;
        vect(3,i) = v_3;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::RightEigenmatrixDotVector()
//! \brief Computes inner-product of right-eigenmatrix of Roe's matrix A in the primitive
//! variables and an input vector.  This operation converts characteristic to primitive
//! variables.  The result is returned in the input vector.
//!
//! The order of the components in the input vector (characteristic fields) should be:
//!    (IDN,ivx,ivy,ivz,[IPR],[IBY,IBZ])
//! where the lower-case indices indicate that the characteristic field in the direction
//! of the sweep (designated by the input flag "ivx") is stored first.  On output, the
//! components of velocity are in the standard order used for primitive variables.
//!
//! REFERENCES:
//! - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//!   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.

void Reconstruction::RightEigenmatrixDotVector(
    const int ivx, const int il, const int iu,
    const AthenaArray<Real> &b1, const AthenaArray<Real> &w, AthenaArray<Real> &vect) {
  // permute components of output primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  if (MAGNETIC_FIELDS_ENABLED) {
    // Adiabatic MHD ---------------------------------------------------------------------
    if (NON_BAROTROPIC_EOS) {
      Real gamma = pmy_block_->peos->GetGamma();
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real id = 1.0/w(IDN,i);
        Real sqrtd = std::sqrt(w(IDN,i));

        Real btsq = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real bxsq = b1(i)*b1(i);
        Real gamp = gamma*w(IPR,i);

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = bxsq + btsq - gamp;
        Real cf2_cs2 = std::sqrt(tdif*tdif + 4.0*gamp*btsq);

        Real cfsq = 0.5*(bxsq + btsq + gamp + cf2_cs2);
        Real cssq = gamp*bxsq/cfsq;

        cfsq *= id;
        Real cf = std::sqrt(cfsq);

        cssq *= id;
        Real cs = std::sqrt(cssq);

        Real asq = gamp*id;
        Real a = std::sqrt(asq);

        // Compute beta(s) (eq A17)
        Real bt  = std::sqrt(btsq);
        Real bet2 = 0.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if ((cfsq - cssq) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (asq - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - asq) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = std::sqrt((asq - cssq)/(cfsq - cssq));
          alpha_s = std::sqrt((cfsq - asq)/(cfsq - cssq));
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real qf = cf*alpha_f*s;
        Real qs = cs*alpha_s*s;
        Real af = a*alpha_f*sqrtd;
        Real as = a*alpha_s*sqrtd;

        // Multiply row of R-eigenmatrix with vector using matrix elements from eq. A12
        // Components of vect() are addressed directly as they are input in permuted order
        Real v_0 = w(IDN,i)*(alpha_f*(vect(0,i) + vect(6,i)) +
                             alpha_s*(vect(2,i) + vect(4,i))) + vect(3,i);
        Real v_1 = cf*alpha_f*(vect(6,i)-vect(0,i)) + cs*alpha_s*(vect(4,i)-vect(2,i));
        Real v_2 = bet2*(qs*(vect(0,i) - vect(6,i)) + qf*(vect(4,i) - vect(2,i)))
                   + bet3*(vect(5,i) - vect(1,i));
        Real v_3 = bet3*(qs*(vect(0,i) - vect(6,i)) + qf*(vect(4,i) - vect(2,i)))
                   + bet2*(vect(1,i) - vect(5,i));
        Real v_4 = w(IDN,i)*asq*(alpha_f*(vect(0,i) + vect(6,i)) +
                                 alpha_s*(vect(2,i) + vect(4,i)));
        Real v_5 = bet2*(as*(vect(0,i) + vect(6,i)) - af*(vect(2,i) + vect(4,i)))
                   - bet3*s*sqrtd*(vect(5,i) + vect(1,i));
        Real v_6 = bet3*(as*(vect(0,i) + vect(6,i)) - af*(vect(2,i) + vect(4,i)))
                   + bet2*s*sqrtd*(vect(5,i) + vect(1,i));

        // Permute components back into standard order for primitives on output
        vect(IDN,i) = v_0;
        vect(ivx,i) = v_1;
        vect(ivy,i) = v_2;
        vect(ivz,i) = v_3;
        vect(IPR,i) = v_4;
        vect(IBY,i) = v_5;
        vect(IBZ,i) = v_6;
      }

      // Isothermal MHD ------------------------------------------------------------------
    } else {
      Real iso_cs = pmy_block_->peos->GetIsoSoundSpeed();
      Real iso_cs2 = SQR(iso_cs);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real id = 1.0/w(IDN,i);
        Real sqrtd = std::sqrt(w(IDN,i));

        Real btsq = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real bxsq = b1(i)*b1(i);
        Real gamp = iso_cs2*w(IDN,i);

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = bxsq + btsq - gamp;
        Real cf2_cs2 = std::sqrt(tdif*tdif + 4.0*gamp*btsq);

        Real cfsq = 0.5*(bxsq + btsq + gamp + cf2_cs2);
        Real cssq = gamp*bxsq/cfsq;

        cfsq *= id;
        Real cf = std::sqrt(cfsq);

        cssq *= id;
        Real cs = std::sqrt(cssq);

        // Compute beta(s) (eq A17)
        Real bt  = std::sqrt(btsq);
        Real bet2 = 0.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if ((cfsq - cssq) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (iso_cs2 - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - iso_cs2) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = std::sqrt((iso_cs2 - cssq)/(cfsq - cssq));
          alpha_s = std::sqrt((cfsq - iso_cs2)/(cfsq - cssq));
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real qf = cf*alpha_f*s;
        Real qs = cs*alpha_s*s;
        Real af = iso_cs*alpha_f*sqrtd;
        Real as = iso_cs*alpha_s*sqrtd;

        // Multiply row of R-eigenmatrix with vector using matrix elements from eq. A12
        // Components of vect() are addressed directly as they are input in permuted order
        Real v_0 = w(IDN,i)*(alpha_f*(vect(0,i) + vect(5,i)) +
                             alpha_s*(vect(2,i) + vect(3,i)));
        Real v_1 = cf*alpha_f*(vect(5,i) - vect(0,i)) + cs*alpha_s*(vect(3,i)-vect(2,i));
        Real v_2 = bet2*(qs*(vect(0,i) - vect(5,i)) + qf*(vect(3,i) - vect(2,i)))
                   + bet3*(vect(4,i) - vect(1,i));
        Real v_3 = bet3*(qs*(vect(0,i) - vect(5,i)) + qf*(vect(3,i) - vect(2,i)))
                   + bet2*(vect(1,i) - vect(4,i));
        Real v_4 = bet2*(as*(vect(0,i) + vect(5,i)) - af*(vect(2,i) + vect(3,i)))
                   - bet3*s*sqrtd*(vect(4,i) + vect(1,i));
        Real v_5 = bet3*(as*(vect(0,i) + vect(5,i)) - af*(vect(2,i) + vect(3,i)))
                   + bet2*s*sqrtd*(vect(4,i) + vect(1,i));

        // Permute components back into standard order for primitives on output
        vect(IDN,i) = v_0;
        vect(ivx,i) = v_1;
        vect(ivy,i) = v_2;
        vect(ivz,i) = v_3;
        vect(IBY,i) = v_4;
        vect(IBZ,i) = v_5;
      }
    }
  } else {
    // Adiabatic hydrodynamics -----------------------------------------------------------
    if (NON_BAROTROPIC_EOS) {
      Real gamma = pmy_block_->peos->GetGamma();
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real asq = gamma*w(IPR,i)/w(IDN,i);
        Real a   = std::sqrt(asq);

        // Multiply row of R-eigenmatrix with vector using matrix elements from eq. A3
        // Components of vect() are addressed directly as they are input in permuted order
        Real v_0 = vect(0,i) + vect(1,i) + vect(4,i);
        Real v_1 = a*(vect(4,i) - vect(0,i))/w(IDN,i);
        Real v_2 = vect(2,i);
        Real v_3 = vect(3,i);
        Real v_4 = asq*(vect(0,i) + vect(4,i));

        // Permute components back into standard order for primitives on output
        vect(IDN,i) = v_0;
        vect(ivx,i) = v_1;
        vect(ivy,i) = v_2;
        vect(ivz,i) = v_3;
        vect(IPR,i) = v_4;
      }

      // Isothermal hydrodynamics --------------------------------------------------------
    } else {
      Real iso_cs = pmy_block_->peos->GetIsoSoundSpeed();
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        // Multiply row of R-eigenmatrix with vector using matrix elements from eq. A3
        // Components of vect() are addressed directly as they are input in permuted order
        Real v_0 = vect(0,i) + vect(3,i);
        Real v_1 = iso_cs*(vect(3,i) - vect(0,i))/w(IDN,i);
        Real v_2 = vect(1,i);
        Real v_3 = vect(2,i);

        // Permute components back into standard order for primitives on output
        vect(IDN,i) = v_0;
        vect(ivx,i) = v_1;
        vect(ivy,i) = v_2;
        vect(ivz,i) = v_3;
      }
    }
  }
  return;
}
