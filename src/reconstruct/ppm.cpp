//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppm.cpp
//  \brief piecewise parabolic reconstruction with modified McCorquodale/Colella limiter
//         for a uniform Cartesian mesh, Mignone limiter for nonuniform mesh
//
// REFERENCES:
// (CW) P. Colella & P. Woodward, "The Piecewise Parabolic Method (PPM) for Gas-Dynamical
// Simulations", JCP, 54, 174 (1984)
//
// (CS) P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at smooth
// extrema", JCP, 227, 7069 (2008)
//
// (MC) P. McCorquodale & P. Colella,  "A high-order finite-volume method for conservation
// laws on locally refined grids", CAMCoS, 6, 1 (2011)
//
// (CD) P. Colella, M.R. Dorr, J. Hittinger, D. Martin, "High-order, finite-volume methods
// in mapped coordinates", JCP, 230, 2952 (2011)
//
// (Mignone) A. Mignone, "High-order conservative reconstruction schemes for finite volume
// methods in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//========================================================================================

// C++ headers
#include <algorithm>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "reconstruction.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX1()
//  \brief Returns L/R interface values in X1-dir constructed using fourth-order PPM and
//         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Reconstruction* prec = pmb->precon;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> bx,wc,q_im2,q_im1,q,q_ip1,q_ip2,qr_imh,ql_iph;
  bx.InitWithShallowCopy(pmb->precon->scr01_i_);
  wc.InitWithShallowCopy(pmb->precon->scr1_ni_);
  q_im2.InitWithShallowCopy(pmb->precon->scr2_ni_);
  q_im1.InitWithShallowCopy(pmb->precon->scr3_ni_);
  q.InitWithShallowCopy(pmb->precon->scr4_ni_);
  q_ip1.InitWithShallowCopy(pmb->precon->scr5_ni_);
  q_ip2.InitWithShallowCopy(pmb->precon->scr6_ni_);
  qr_imh.InitWithShallowCopy(pmb->precon->scr7_ni_);
  ql_iph.InitWithShallowCopy(pmb->precon->scr8_ni_);

  // set work PPM work arrays to shallow copies of scratch arrays:
  AthenaArray<Real> dd,dd_im1,dd_ip1,dph,dph_ip1;
  dd.InitWithShallowCopy(pmb->precon->scr02_i_);
  dd_im1.InitWithShallowCopy(pmb->precon->scr03_i_);
  dd_ip1.InitWithShallowCopy(pmb->precon->scr04_i_);
  dph.InitWithShallowCopy(pmb->precon->scr05_i_);
  dph_ip1.InitWithShallowCopy(pmb->precon->scr06_i_);

  AthenaArray<Real> d2qc_im1,d2qc,d2qc_ip1,d2qf;
  d2qc_im1.InitWithShallowCopy(pmb->precon->scr07_i_);
  d2qc.InitWithShallowCopy(pmb->precon->scr08_i_);
  d2qc_ip1.InitWithShallowCopy(pmb->precon->scr09_i_);
  d2qf.InitWithShallowCopy(pmb->precon->scr10_i_);

  AthenaArray<Real> qplus,qminus,dqf_plus,dqf_minus;
  qplus.InitWithShallowCopy(pmb->precon->scr11_i_);
  qminus.InitWithShallowCopy(pmb->precon->scr12_i_);
  dqf_plus.InitWithShallowCopy(pmb->precon->scr13_i_);
  dqf_minus.InitWithShallowCopy(pmb->precon->scr14_i_);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    // cache the x1-sliced primitive states for eigensystem calculation
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        wc(n,i) = w(n,k,j,i);
        q    (n,i) = w(n,k,j,i  );
        q_im2(n,i) = w(n,k,j,i-2);
        q_im1(n,i) = w(n,k,j,i-1);
        q_ip1(n,i) = w(n,k,j,i+1);
        q_ip2(n,i) = w(n,k,j,i+2);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        bx(i) = bcc(IB1,k,j,i);

        wc(IBY,i) = bcc(IB2,k,j,i);
        q    (IBY,i) = bcc(IB2,k,j,i  );
        q_im2(IBY,i) = bcc(IB2,k,j,i-2);
        q_im1(IBY,i) = bcc(IB2,k,j,i-1);
        q_ip1(IBY,i) = bcc(IB2,k,j,i+1);
        q_ip2(IBY,i) = bcc(IB2,k,j,i+2);

        wc(IBZ,i) = bcc(IB3,k,j,i);
        q    (IBZ,i) = bcc(IB3,k,j,i  );
        q_im2(IBZ,i) = bcc(IB3,k,j,i-2);
        q_im1(IBZ,i) = bcc(IB3,k,j,i-1);
        q_ip1(IBZ,i) = bcc(IB3,k,j,i+1);
        q_ip2(IBZ,i) = bcc(IB3,k,j,i+2);
      }
    }

    // Project cell-averages to characteristic variables, if necessary
    // Note order of characteristic fields in output vect corresponds to (IVX,IVY,IVZ)
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,q_im2);
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,q_im1);
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,q);
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,q_ip1);
      LeftEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,q_ip2);
    }

//--- Step 1. ----------------------------------------------------------------------------
// Reconstruct interface averages <a>_{i-1/2} and <a>_{i+1/2}
    for (int n=0; n<(NWAVE); ++n) {

      // Compute average slope in i-1, i, i+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il-1; i<=iu; ++i) {
        Real qa = (q(n,i) - q_im1(n,i));
        Real qb = (q_ip1(n,i) - q(n,i));
        dd_im1(i) = prec->c1i(i-1)*qa + prec->c2i(i-1)*(q_im1(n,i) - q_im2(n,i));
        dd    (i) = prec->c1i(i  )*qb + prec->c2i(i  )*qa;
        dd_ip1(i) = prec->c1i(i+1)*(q_ip2(n,i) - q_ip1(n,i)) + prec->c2i(i+1)*qb;

        // Approximate interface average at i-1/2 and i+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph(i)= (prec->c3i(i)*q_im1(n,i) + prec->c4i(i)*q(n,i)) +
            (prec->c5i(i)*dd_im1(i) + prec->c6i(i)*dd(i));
        dph_ip1(i)= (prec->c3i(i+1)*q(n,i) + prec->c4i(i+1)*q_ip1(n,i)) +
            (prec->c5i(i+1)*dd(i) + prec->c6i(i+1)*dd_ip1(i) );
      }

//--- Step 2a. ---------------------------------------------------------------------------
      // Uniform Cartesian grid: limit interpolated interface states as in CD 4.3.1
      if (pmb->precon->uniform_limiter[X1DIR]) {
        // approximate second derivative at interfaces for smooth extrema preservation
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il-1; i<=iu+1; ++i) {
          // KGF: add the off-centered quantities first to preserve FP symmetry
          d2qc_im1(i) = q_im2(n,i) + q    (n,i) - 2.0*q_im1(n,i);
          d2qc    (i) = q_im1(n,i) + q_ip1(n,i) - 2.0*q    (n,i); //(CD eq 85a) (no 1/2)
          d2qc_ip1(i) = q    (n,i) + q_ip2(n,i) - 2.0*q_ip1(n,i);
        }

        // i-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il-1; i<=(iu+1); ++i) {
          Real qa_tmp = dph(i) - q_im1(n,i); // (CD eq 84a)
          Real qb_tmp = q(n,i) - dph(i);     // (CD eq 84b)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real qa = 3.0*(q_im1(n,i) + q(n,i)  - 2.0*dph(i));  // (CD eq 85b)
          Real qb = d2qc_im1(i);    // (CD eq 85a) (no 1/2)
          Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
          Real qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          Real dph_tmp = 0.5*(q_im1(n,i)+q(n,i)) - qd/6.0;
          if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at i-1/2 face
            dph(i) = dph_tmp;
          }
        }

        // i+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il-1; i<=(iu+1); ++i) {
          Real qa_tmp = dph_ip1(i) - q(n,i);       // (CD eq 84a)
          Real qb_tmp = q_ip1(n,i) - dph_ip1(i);   // (CD eq 84b)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real qa = 3.0*(q(n,i) + q_ip1(n,i) - 2.0*dph_ip1(i));  // (CD eq 85b)
          Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
          Real qc = d2qc_ip1(i);   // (CD eq 85c) (no 1/2)
          Real qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          Real dphip1_tmp = 0.5*(q(n,i)+q_ip1(n,i)) - qd/6.0;
          if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at i+1/2 face
            dph_ip1(i) = dphip1_tmp;
          }
        }

#pragma omp simd
        for (int i=il-1; i<=iu; ++i) {
          // KGF: add the off-centered quantities first to preserve FP symmetry
          d2qf(i) = 6.0*(dph(i) + dph_ip1(i) - 2.0*q(n,i)); // a6 coefficient * -2
        }

//--- Step 2b. ---------------------------------------------------------------------------
      // Non-uniform/curvilinear: apply strict monotonicity constraints (Mignone eq 45)
      } else {
#pragma omp simd
        for (int i=il-1; i<=iu; ++i) {
          dph    (i) = std::min(dph    (i), std::max(q(n,i),q_im1(n,i)));
          dph_ip1(i) = std::min(dph_ip1(i), std::max(q(n,i),q_ip1(n,i)));

          dph    (i) = std::max(dph    (i), std::min(q(n,i),q_im1(n,i)));
          dph_ip1(i) = std::max(dph_ip1(i), std::min(q(n,i),q_ip1(n,i)));
        }
      }

      // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        qminus(i) = dph(i  );
        qplus(i) =  dph_ip1(i );
      }

//--- Step 3. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        dqf_minus(i) = q(n,i) - qminus(i); // (CS eq 25)
        dqf_plus(i)  = qplus(i) - q(n,i);
      }

//--- Step 4a. ---------------------------------------------------------------------------
      // For uniform Cartesian mesh: apply CS limiters to parabolic interpolant
      if (pmb->precon->uniform_limiter[X1DIR]) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il-1; i<=iu; ++i) {
          Real qa_tmp = dqf_minus(i)*dqf_plus(i);
          Real qb_tmp = (q_ip1(n,i) - q(n,i))*(q(n,i) - q_im1(n,i));

          Real qa = d2qc_im1(i);
          Real qb = d2qc(i);
          Real qc = d2qc_ip1(i);
          Real qd = d2qf(i);
          Real qe = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          }

          // Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q_im1(n,i)),fabs(q_im2(n,i)));
          qb = std::max(std::max(fabs(q(n,i)),fabs(q_ip1(n,i))), fabs(q_ip2(n,i)));

          Real rho = 0.0;
          if (fabs(qd) > (1.0e-12)*std::max(qa,qb)) {
            // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
            rho = qe/qd;
          }

          Real tmp_m = q(n,i) - rho*dqf_minus(i);
          Real tmp_p = q(n,i) + rho*dqf_plus(i);
          Real tmp2_m = q(n,i) - 2.0*dqf_plus(i);
          Real tmp2_p = q(n,i) + 2.0*dqf_minus(i);

          // Check for local extrema
          if ((qa_tmp <= 0.0 || qb_tmp <=0.0)) {
            // Check if relative change in limited 2nd deriv is > roundoff
            if (rho <= (1.0 - (1.0e-12))) {
              // Limit smooth extrema
              qminus(i) = tmp_m;// (CS eq 23)
              qplus(i) = tmp_p;
            }
            // No extrema detected
          } else {
            // Overshoot i-1/2,R / i,(-) state
            if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
              qminus(i) = tmp2_m;
            }
            // Overshoot i+1/2,L / i,(+) state
            if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
              qplus(i) = tmp2_p;
            }
          }
        }

//--- Step 4b. ---------------------------------------------------------------------------
      // Non-uniform/curvilinear mesh: apply Mignone limiters to parabolic interpolant
      // Note, the Mignone limiter does not check for cell-averaged extrema:
      } else {
        for (int i=il-1; i<=iu; ++i) {
          Real qa = dqf_minus(i)*dqf_plus(i);
          if (qa <= 0.0) { // Local extrema detected
            qminus(i) = q(n,i);
            qplus(i) = q(n,i);
          } else { // No extrema detected
            // Overshoot i-1/2,R / i,(-) state
            if (fabs(dqf_minus(i)) >= prec->hplus_ratio_i(i)*fabs(dqf_plus(i))) {
              qminus(i) = q(n,i) - prec->hplus_ratio_i(i)*dqf_plus(i);
            }
            // Overshoot i+1/2,L / i,(+) state
            if (fabs(dqf_plus(i)) >= prec->hminus_ratio_i(i)*fabs(dqf_minus(i))) {
              qplus(i) = q(n,i) + prec->hminus_ratio_i(i)*dqf_minus(i);
            }
          }
        }
      }

//--- Step 5. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [il,iu]
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        ql_iph(n,i ) = qplus(i);
        qr_imh(n,i ) = qminus(i);
      }
    } // end char PPM loop over NWAVE

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      RightEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,ql_iph);
      RightEigenmatrixDotVector(pmb,IVX,il-1,iu,bx,wc,qr_imh);
    }

    // compute ql_(i+1/2) and qr_(i-1/2)
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd
      for (int i=il-1; i<=iu; ++i) {
        wl(n,k,j,i+1) = ql_iph(n,i);
        wr(n,k,j,i  ) = qr_imh(n,i);
        // Reapply EOS floors to both L/R reconstructed primitive states
        pmb->peos->ApplyPrimitiveFloors(wl, k, j, i+1);
        pmb->peos->ApplyPrimitiveFloors(wr, k, j, i);
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX2()
//  \brief Returns L/R interface values in X2-dir constructed using fourth-order PPM and
//         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX2(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Reconstruction* prec = pmb->precon;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> bx,wc,q_jm2,q_jm1,q,q_jp1,q_jp2,qr_jmh,ql_jph;
  bx.InitWithShallowCopy(pmb->precon->scr01_i_);
  wc.InitWithShallowCopy(pmb->precon->scr1_ni_);
  q_jm2.InitWithShallowCopy(pmb->precon->scr2_ni_);
  q_jm1.InitWithShallowCopy(pmb->precon->scr3_ni_);
  q.InitWithShallowCopy(pmb->precon->scr4_ni_);
  q_jp1.InitWithShallowCopy(pmb->precon->scr5_ni_);
  q_jp2.InitWithShallowCopy(pmb->precon->scr6_ni_);
  qr_jmh.InitWithShallowCopy(pmb->precon->scr7_ni_);
  ql_jph.InitWithShallowCopy(pmb->precon->scr8_ni_);

  // set work PPM work arrays to shallow copies of scratch arrays:
  AthenaArray<Real> dd,dd_jm1,dd_jp1,dph,dph_jp1;
  dd.InitWithShallowCopy(pmb->precon->scr02_i_);
  dd_jm1.InitWithShallowCopy(pmb->precon->scr03_i_);
  dd_jp1.InitWithShallowCopy(pmb->precon->scr04_i_);
  dph.InitWithShallowCopy(pmb->precon->scr05_i_);
  dph_jp1.InitWithShallowCopy(pmb->precon->scr06_i_);

  AthenaArray<Real> d2qc_jm1,d2qc,d2qc_jp1,d2qf;
  d2qc_jm1.InitWithShallowCopy(pmb->precon->scr07_i_);
  d2qc.InitWithShallowCopy(pmb->precon->scr08_i_);
  d2qc_jp1.InitWithShallowCopy(pmb->precon->scr09_i_);
  d2qf.InitWithShallowCopy(pmb->precon->scr10_i_);

  AthenaArray<Real> qplus,qminus,dqf_plus,dqf_minus;
  qplus.InitWithShallowCopy(pmb->precon->scr11_i_);
  qminus.InitWithShallowCopy(pmb->precon->scr12_i_);
  dqf_plus.InitWithShallowCopy(pmb->precon->scr13_i_);
  dqf_minus.InitWithShallowCopy(pmb->precon->scr14_i_);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl-1; j<=ju; ++j) {
    // cache the x1-sliced primitive states for eigensystem calculation
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wc(n,i) = w(n,k,j,i);
        q    (n,i) = w(n,k,j  ,i);
        q_jm2(n,i) = w(n,k,j-2,i);
        q_jm1(n,i) = w(n,k,j-1,i);
        q_jp1(n,i) = w(n,k,j+1,i);
        q_jp2(n,i) = w(n,k,j+2,i);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        bx(i) = bcc(IB2,k,j,i);

        wc(IBY,i) = bcc(IB3,k,j,i);
        q    (IBY,i) = bcc(IB3,k,j  ,i);
        q_jm2(IBY,i) = bcc(IB3,k,j-2,i);
        q_jm1(IBY,i) = bcc(IB3,k,j-1,i);
        q_jp1(IBY,i) = bcc(IB3,k,j+1,i);
        q_jp2(IBY,i) = bcc(IB3,k,j+2,i);

        wc(IBZ,i) = bcc(IB1,k,j,i);
        q    (IBZ,i) = bcc(IB1,k,j  ,i);
        q_jm2(IBZ,i) = bcc(IB1,k,j-2,i);
        q_jm1(IBZ,i) = bcc(IB1,k,j-1,i);
        q_jp1(IBZ,i) = bcc(IB1,k,j+1,i);
        q_jp2(IBZ,i) = bcc(IB1,k,j+2,i);
      }
    }

    // Project cell-averages to characteristic variables, if necessary
    // Note order of characteristic fields in output vect corresponds to (IVY,IVZ,IVX)
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,q_jm2);
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,q_jm1);
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,q);
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,q_jp1);
      LeftEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,q_jp2);
    }

//--- Step 1. ----------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} and <a>_{j+1/2}
    for (int n=0; n<(NWAVE); ++n) {

      // Compute average slope in j-1, j, j+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa = (q(n,i) - q_jm1(n,i));
        Real qb = (q_jp1(n,i) - q(n,i));
        dd_jm1(i) = prec->c1j(j-1)*qa + prec->c2j(j-1)*(q_jm1(n,i) - q_jm2(n,i));
        dd    (i) = prec->c1j(j  )*qb + prec->c2j(j  )*qa;
        dd_jp1(i) = prec->c1j(j+1)*(q_jp2(n,i) - q_jp1(n,i)) + prec->c2j(j+1)*qb;

        // Approximate interface average at j-1/2 and j+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph(i)= (prec->c3j(j)*q_jm1(n,i) + prec->c4j(j)*q(n,i)) +
            (prec->c5j(j)*dd_jm1(i) + prec->c6j(j)*dd(i));
        dph_jp1(i)= (prec->c3j(j+1)*q(n,i) + prec->c4j(j+1)*q_jp1(n,i)) +
            (prec->c5j(j+1)*dd(i) + prec->c6j(j+1)*dd_jp1(i));
      }

//--- Step 2a. ---------------------------------------------------------------------------
      // Uniform Cartesian grid: limit interpolated interface states as in CD 4.3.1
      if (pmb->precon->uniform_limiter[X2DIR]) {
        // approximate second derivative at interfaces for smooth extrema preservation
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          // KGF: add the off-centered quantities first to preserve FP symmetry
          d2qc_jm1(i) = q_jm2(n,i) + q    (n,i) - 2.0*q_jm1(n,i);
          d2qc    (i) = q_jm1(n,i) + q_jp1(n,i) - 2.0*q    (n,i); //(CD eq 85a) (no 1/2)
          d2qc_jp1(i) = q    (n,i) + q_jp2(n,i) - 2.0*q_jp1(n,i);
        }

        // j-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          Real qa_tmp = dph(i) - q_jm1(n,i); // (CD eq 84a)
          Real qb_tmp = q(n,i) - dph(i);     // (CD eq 84b)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real qa = 3.0*(q_jm1(n,i) + q(n,i) - 2.0*dph(i));  // (CD eq 85b)
          Real qb = d2qc_jm1(i);    // (CD eq 85a) (no 1/2)
          Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
          Real qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          Real dph_tmp = 0.5*(q_jm1(n,i)+q(n,i)) - qd/6.0;
          if (qa_tmp*qb_tmp < 0.0) { //Local extrema detected at j-1/2 face
            dph(i) = dph_tmp;
          }
        }
        // j+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          Real qa_tmp = dph_jp1(i) - q(n,i);       // (CD eq 84a)
          Real qb_tmp = q_jp1(n,i) - dph_jp1(i);   // (CD eq 84b)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real qa = 3.0*(q(n,i) + q_jp1(n,i)  - 2.0*dph_jp1(i));  // (CD eq 85b)
          Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
          Real qc = d2qc_jp1(i);   // (CD eq 85c) (no 1/2)
          Real qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          Real dphjp1_tmp = 0.5*(q(n,i)+q_jp1(n,i)) - qd/6.0;
          if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at j+1/2 face
            dph_jp1(i) = dphjp1_tmp;
          }
        }

#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // KGF: add the off-centered quantities first to preserve FP symmetry
          d2qf(i) = 6.0*(dph(i) + dph_jp1(i) - 2.0*q(n,i)); // a6 coefficient * -2
        }

//--- Step 2b. ---------------------------------------------------------------------------
      // Non-uniform/curvilinear: apply strict monotonicity constraints (Mignone eq 45)
      } else {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          dph    (i) = std::min(dph    (i), std::max(q(n,i),q_jm1(n,i)));
          dph_jp1(i) = std::min(dph_jp1(i), std::max(q(n,i),q_jp1(n,i)));

          dph    (i) = std::max(dph    (i), std::min(q(n,i),q_jm1(n,i)));
          dph_jp1(i) = std::max(dph_jp1(i), std::min(q(n,i),q_jp1(n,i)));
        }
      }

      // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        qminus(i) = dph(i  );
        qplus(i) =  dph_jp1(i );
      }

//--- Step 3. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dqf_minus(i) = q(n,i) - qminus(i); // (CS eq 25)
        dqf_plus(i)  = qplus(i) - q(n,i);
      }

//--- Step 4a. ---------------------------------------------------------------------------
      // For uniform Cartesian mesh: apply CS limiters to parabolic interpolant
      if (pmb->precon->uniform_limiter[X2DIR]) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          Real qa_tmp = dqf_minus(i)*dqf_plus(i);
          Real qb_tmp = (q_jp1(n,i) - q(n,i))*(q(n,i) - q_jm1(n,i));

          Real qa = d2qc_jm1(i);
          Real qb = d2qc(i);
          Real qc = d2qc_jp1(i);
          Real qd = d2qf(i);
          Real qe = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          }

          // Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q_jm1(n,i)),fabs(q_jm2(n,i)));
          qb = std::max(std::max(fabs(q(n,i)),fabs(q_jp1(n,i))), fabs(q_jp2(n,i)));

          Real rho = 0.0;
          if (fabs(qd) > (1.0e-12)*std::max(qa,qb)) {
              // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
            rho = qe/qd;
          }

          Real tmp_m = q(n,i) - rho*dqf_minus(i);
          Real tmp_p = q(n,i) + rho*dqf_plus(i);
          Real tmp2_m = q(n,i) - 2.0*dqf_plus(i);
          Real tmp2_p = q(n,i) + 2.0*dqf_minus(i);

          // Check if relative change in limited 2nd deriv is > roundoff
          // Check for local extrema
          if ((qa_tmp <= 0.0 || qb_tmp <=0.0)) {
            if (rho <= (1.0 - (1.0e-12))) {
              // Limit smooth extrema
              qminus(i) = tmp_m; // (CS eq 23)
              qplus(i) = tmp_p;
            }
            // No extrema detected
          } else {
            // Overshoot j-1/2,R / j,(-) state
            if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
              qminus(i) = tmp2_m;
            }
            // Overshoot j+1/2,L / j,(+) state
            if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
              qplus(i) = tmp2_p;
            }
          }
        }

//--- Step 4b. ---------------------------------------------------------------------------
      // Non-uniform/curvilinear mesh: apply Mignone limiters to parabolic interpolant
      // Note, the Mignone limiter does not check for cell-averaged extrema:
      } else {
        for (int i=il; i<=iu; ++i) {
          Real qa = dqf_minus(i)*dqf_plus(i);
          if (qa <= 0.0) { // Local extrema detected
            qminus(i) = q(n,i);
            qplus(i) = q(n,i);
          } else { // No extrema detected
            // Overshoot j-1/2,R / j,(-) state
            if (fabs(dqf_minus(i)) >= prec->hplus_ratio_j(j)*fabs(dqf_plus(i))) {
              qminus(i) = q(n,i) - prec->hplus_ratio_j(j)*dqf_plus(i);
            }
            // Overshoot j+1/2,L / j,(+) state
            if (fabs(dqf_plus(i)) >= prec->hminus_ratio_j(j)*fabs(dqf_minus(i))) {
              qplus(i) = q(n,i) + prec->hminus_ratio_j(j)*dqf_minus(i);
            }
          }
        }
      }

//--- Step 5. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [il,iu]
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        ql_jph(n,i ) = qplus(i);
        qr_jmh(n,i ) = qminus(i);
      }
    } // end char PPM loop over NWAVE

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      RightEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,ql_jph);
      RightEigenmatrixDotVector(pmb,IVY,il,iu,bx,wc,qr_jmh);
    }

    // compute ql_(j+1/2) and qr_(j-1/2)
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,k,j+1,i) = ql_jph(n,i);
        wr(n,k,j  ,i) = qr_jmh(n,i);
        // Reapply EOS floors to both L/R reconstructed primitive states
        pmb->peos->ApplyPrimitiveFloors(wl, k, j+1, i);
        pmb->peos->ApplyPrimitiveFloors(wr, k, j, i);
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX3()
//  \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX3(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  Reconstruction* prec = pmb->precon;
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;

  // set work arrays used for primitive/characterstic cell-averages to scratch
  AthenaArray<Real> bx,wc,q_km2,q_km1,q,q_kp1,q_kp2,qr_kmh,ql_kph;
  bx.InitWithShallowCopy(pmb->precon->scr01_i_);
  wc.InitWithShallowCopy(pmb->precon->scr1_ni_);
  q_km2.InitWithShallowCopy(pmb->precon->scr2_ni_);
  q_km1.InitWithShallowCopy(pmb->precon->scr3_ni_);
  q.InitWithShallowCopy(pmb->precon->scr4_ni_);
  q_kp1.InitWithShallowCopy(pmb->precon->scr5_ni_);
  q_kp2.InitWithShallowCopy(pmb->precon->scr6_ni_);
  qr_kmh.InitWithShallowCopy(pmb->precon->scr7_ni_);
  ql_kph.InitWithShallowCopy(pmb->precon->scr8_ni_);

  // set work PPM work arrays to shallow copies of scratch arrays:
  AthenaArray<Real> dd,dd_km1,dd_kp1,dph,dph_kp1;
  dd.InitWithShallowCopy(pmb->precon->scr02_i_);
  dd_km1.InitWithShallowCopy(pmb->precon->scr03_i_);
  dd_kp1.InitWithShallowCopy(pmb->precon->scr04_i_);
  dph.InitWithShallowCopy(pmb->precon->scr05_i_);
  dph_kp1.InitWithShallowCopy(pmb->precon->scr06_i_);

  AthenaArray<Real> d2qc_km1,d2qc,d2qc_kp1,d2qf;
  d2qc_km1.InitWithShallowCopy(pmb->precon->scr07_i_);
  d2qc.InitWithShallowCopy(pmb->precon->scr08_i_);
  d2qc_kp1.InitWithShallowCopy(pmb->precon->scr09_i_);
  d2qf.InitWithShallowCopy(pmb->precon->scr10_i_);

  AthenaArray<Real> qplus,qminus,dqf_plus,dqf_minus;
  qplus.InitWithShallowCopy(pmb->precon->scr11_i_);
  qminus.InitWithShallowCopy(pmb->precon->scr12_i_);
  dqf_plus.InitWithShallowCopy(pmb->precon->scr13_i_);
  dqf_minus.InitWithShallowCopy(pmb->precon->scr14_i_);

  for (int k=kl-1; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
    // cache the x1-sliced primitive states for eigensystem calculation
    for (int n=0; n<(NHYDRO); ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wc(n,i) = w(n,k,j,i);
        q    (n,i) = w(n,k  ,j,i);
        q_km2(n,i) = w(n,k-2,j,i);
        q_km1(n,i) = w(n,k-1,j,i);
        q_kp1(n,i) = w(n,k+1,j,i);
        q_kp2(n,i) = w(n,k+2,j,i);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        bx(i) = bcc(IB3,k,j,i);

        wc(IBY,i) = bcc(IB1,k,j,i);
        q    (IBY,i) = bcc(IB1,k  ,j,i);
        q_km2(IBY,i) = bcc(IB1,k-2,j,i);
        q_km1(IBY,i) = bcc(IB1,k-1,j,i);
        q_kp1(IBY,i) = bcc(IB1,k+1,j,i);
        q_kp2(IBY,i) = bcc(IB1,k+2,j,i);

        wc(IBZ,i) = bcc(IB2,k,j,i);
        q    (IBZ,i) = bcc(IB2,k  ,j,i);
        q_km2(IBZ,i) = bcc(IB2,k-2,j,i);
        q_km1(IBZ,i) = bcc(IB2,k-1,j,i);
        q_kp1(IBZ,i) = bcc(IB2,k+1,j,i);
        q_kp2(IBZ,i) = bcc(IB2,k+2,j,i);
      }
    }

    // Project cell-averages to characteristic variables, if necessary
    // Note order of characteristic fields in output vect corresponds to (IVZ,IVX,IVY)
    if (pmb->precon->characteristic_reconstruction) {
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,q_km2);
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,q_km1);
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,q);
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,q_kp1);
      LeftEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,q_kp2);
    }

//--- Step 1. ----------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} and <a>_{k+1/2}
    for (int n=0; n<(NWAVE); ++n) {

      // Compute average slope in k-1, k, k+1 zones
#pragma omp simd simdlen(SIMD_WIDTH)
      for (int i=il; i<=iu; ++i) {
        Real qa = (q(n,i) - q_km1(n,i));
        Real qb = (q_kp1(n,i) - q(n,i));
        dd_km1(i) = prec->c1k(k-1)*qa + prec->c2k(k-1)*(q_km1(n,i) - q_km2(n,i));
        dd    (i) = prec->c1k(k  )*qb + prec->c2k(k  )*qa;
        dd_kp1(i) = prec->c1k(k+1)*(q_kp2(n,i) - q_kp1(n,i)) + prec->c2k(k+1)*qb;

        // Approximate interface average at k-1/2 and k+1/2 using PPM (CW eq 1.6)
        // KGF: group the biased stencil quantities to preserve FP symmetry
        dph(i)= (prec->c3k(k)*q_km1(n,i) + prec->c4k(k)*q(n,i)) +
          (prec->c5k(k)*dd_km1(i) + prec->c6k(k)*dd(i));
        dph_kp1(i)= (prec->c3k(k+1)*q(n,i) + prec->c4k(k+1)*q_kp1(n,i)) +
          (prec->c5k(k+1)*dd(i) + prec->c6k(k+1)*dd_kp1(i));
      }

//--- Step 2a. ---------------------------------------------------------------------------
      // Uniform Cartesian grid: limit interpolated interface states as in CD 4.3.1
      if (pmb->precon->uniform_limiter[X3DIR]) {
        // approximate second derivative at interfaces for smooth extrema preservation
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          // KGF: add the off-centered quantities first to preserve FP symmetry
          d2qc_km1(i) = q_km2(n,i) + q    (n,i) - 2.0*q_km1(n,i) ;
          d2qc    (i) = q_km1(n,i) + q_kp1(n,i) - 2.0*q    (n,i) ; //(CD eq 85a) (no 1/2)
          d2qc_kp1(i) = q    (n,i) + q_kp2(n,i) - 2.0*q_kp1(n,i) ;
        }

        // k-1/2
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          Real qa_tmp = dph(i) - q_km1(n,i); // (CD eq 84a)
          Real qb_tmp = q(n,i) - dph(i);     // (CD eq 84b)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real qa = 3.0*(q_km1(n,i) + q(n,i) - 2.0*dph(i));  // (CD eq 85b)
          Real qb = d2qc_km1(i);    // (CD eq 85a) (no 1/2)
          Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
          Real qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          Real dph_tmp = 0.5*(q_km1(n,i)+q(n,i)) - qd/6.0;
          if (qa_tmp*qb_tmp < 0.0) {  // Local extrema detected at k-1/2 face
            dph(i) = dph_tmp;
          }
        }
        // k+1/2
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          Real qa_tmp = dph_kp1(i) - q(n,i);       // (CD eq 84a)
          Real qb_tmp = q_kp1(n,i) - dph_kp1(i);   // (CD eq 84b)
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real qa = 3.0*(q(n,i) + q_kp1(n,i) - 2.0*dph_kp1(i));  // (CD eq 85b)
          Real qb = d2qc(i);            // (CD eq 85a) (no 1/2)
          Real qc = d2qc_kp1(i);   // (CD eq 85c) (no 1/2)
          Real qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          Real dphkp1_tmp = 0.5*(q(n,i)+q_kp1(n,i)) - qd/6.0;
          if (qa_tmp*qb_tmp < 0.0) { // Local extrema detected at k+1/2 face
            dph_kp1(i) = dphkp1_tmp;
          }
        }

#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // KGF: add the off-centered quantities first to preserve FP symmetry
          d2qf(i) = 6.0*(dph(i) + dph_kp1(i) - 2.0*q(n,i)); // a6 coefficient * -2
        }

//--- Step 2b. ---------------------------------------------------------------------------
      // Non-uniform/curvilinear: apply strict monotonicity constraints (Mignone eq 45)
      } else {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          dph    (i) = std::min(dph    (i), std::max(q(n,i),q_km1(n,i)));
          dph_kp1(i) = std::min(dph_kp1(i), std::max(q(n,i),q_kp1(n,i)));

          dph    (i) = std::max(dph    (i), std::min(q(n,i),q_km1(n,i)));
          dph_kp1(i) = std::max(dph_kp1(i), std::min(q(n,i),q_kp1(n,i)));
        }
      }

      // Cache Riemann states for both non-/uniform limiters
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        qminus(i) = dph(i  );
        qplus(i) =  dph_kp1(i );
      }

//--- Step 3. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        dqf_minus(i) = q(n,i) - qminus(i); // (CS eq 25)
        dqf_plus(i)  = qplus(i) - q(n,i);
      }

//--- Step 4a. ---------------------------------------------------------------------------
      // For uniform Cartesian mesh: apply CS limiters to parabolic interpolant
      if (pmb->precon->uniform_limiter[X3DIR]) {
#pragma omp simd simdlen(SIMD_WIDTH)
        for (int i=il; i<=iu; ++i) {
          Real qa_tmp = dqf_minus(i)*dqf_plus(i);
          Real qb_tmp = (q_kp1(n,i) - q(n,i))*(q(n,i) - q_km1(n,i));

          // Check if extrema is smooth
          Real qa = d2qc_km1(i);
          Real qb = d2qc(i);
          Real qc = d2qc_kp1(i);
          Real qd = d2qf(i);
          Real qe = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
              // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          }

            // Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q_km1(n,i)),fabs(q_km2(n,i)));
          qb = std::max(std::max(fabs(q(n,i)),fabs(q_kp1(n,i))), fabs(q_kp2(n,i)));

          Real rho = 0.0;
          if (fabs(qd) > (1.0e-12)*std::max(qa,qb)) {
            // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
            rho = qe/qd;
          }

          Real tmp_m = q(n,i) - rho*dqf_minus(i);
          Real tmp_p = q(n,i) + rho*dqf_plus(i);
          Real tmp2_m = q(n,i) - 2.0*dqf_plus(i);
          Real tmp2_p = q(n,i) + 2.0*dqf_minus(i);

          // Check for local extrema
          if (qa_tmp <= 0.0 || qb_tmp <= 0.0 ) {
            // Check if relative change in limited 2nd deriv is > roundoff
            if (rho <= (1.0 - (1.0e-12))) {
              // Limit smooth extrema
              qminus(i) = tmp_m; // (CS eq 23)
              qplus(i) = tmp_p;
            }

          // No extrema detected
          } else {
            // Overshoot k-1/2,R / k,(-) state
            if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
              qminus(i) = tmp2_m;
            }
            // Overshoot k+1/2,L / k,(+) state
            if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
              qplus(i) = tmp2_p;
            }
          }
        }

//--- Step 4b. ---------------------------------------------------------------------------
      // Non-uniform/curvilinear mesh: apply Mignone limiters to parabolic interpolant
      // Note, the Mignone limiter does not check for cell-averaged extrema:
      } else {
        for (int i=il; i<=iu; ++i) {
          Real qa = dqf_minus(i)*dqf_plus(i);
          if (qa <= 0.0) { // Local extrema detected
            qminus(i) = q(n,i);
            qplus(i) = q(n,i);
          } else { // No extrema detected
            // could delete hplus_ratio_k() arrays for curvilinear PPMx3
            // Overshoot k-1/2,R / k,(-) state
            if (fabs(dqf_minus(i)) >= prec->hplus_ratio_k(k)*fabs(dqf_plus(i))) {
              qminus(i) = q(n,i) - prec->hplus_ratio_k(k)*dqf_plus(i);
            }
            // Overshoot k+1/2,L / k,(+) state
            if (fabs(dqf_plus(i)) >= prec->hminus_ratio_k(k)*fabs(dqf_minus(i))) {
              qplus(i) = q(n,i) + prec->hminus_ratio_k(k)*dqf_minus(i);
            }
          }
        }
      }

//--- Step 5. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [il,iu]
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        ql_kph(n,i ) = qplus(i);
        qr_kmh(n,i ) = qminus(i);
      }
    } // end char PPM loop over NWAVE

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      RightEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,ql_kph);
      RightEigenmatrixDotVector(pmb,IVZ,il,iu,bx,wc,qr_kmh);
    }

    // compute ql_(k+1/2) and qr_(k-1/2)
    for (int n=0; n<(NWAVE); ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,k+1,j,i) = ql_kph(n,i);
        wr(n,k  ,j,i) = qr_kmh(n,i);
        // Reapply EOS floors to both L/R reconstructed primitive states
        pmb->peos->ApplyPrimitiveFloors(wl, k+1, j, i);
        pmb->peos->ApplyPrimitiveFloors(wr, k, j, i);
      }
    }
  }}

  return;
}
