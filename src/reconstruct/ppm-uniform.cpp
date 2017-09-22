//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppm-uniform.cpp
//  \brief piecewise parabolic reconstruction with modified McCorquodale/Colella limiter
//         for a uniform Cartesian mesh.
//
// REFERENCES:(MC) P. McCorquodale & P. Colella,  "A high-order finite-volume method for
// conservation laws on locally refined grids", CAMCoS, 6, 1 (2011)
//
// (CW) P. Colella & P. Woodward, "The Piecewise Parabolic Method (PPM) for Gas-Dynamical
// Simulations", JCP, 54, 174 (1984)
//
// (CS) P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at smooth
// extrema", JCP, 227, 7069 (2008)
//
// (CD) P. Colella, M.R. Dorr, J. Hittinger, D. Martin, "High-order, finite-volume methods
// in mapped coordinates", JCP, 230, 2952 (2011)
//========================================================================================

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX1()
//  \brief Returns L/R interface values in X1-dir constructed using fourth-order PPM and
//         Colella-Sekora extremum-preserving limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PPMUniformX1(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;
  Real qa,qb,qc,qd,qe,rho;

  // 1D scratch arrays
  AthenaArray<Real> dph, qplus, qminus, dqf_plus, dqf_minus, d2qf, d2qc;
  dph.InitWithShallowCopy(dph_);
  qplus.InitWithShallowCopy(qplus_);
  qminus.InitWithShallowCopy(qminus_);
  dqf_plus.InitWithShallowCopy(dqf_plus_);
  dqf_plus.InitWithShallowCopy(dqf_minus_);
  d2qf.InitWithShallowCopy(d2qf_);
  d2qc.InitWithShallowCopy(d2qc_);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {

//--- Step 1. ----------------------------------------------------------------------------
// Reconstruct interface averages <a>_{i-1/2} using PPM for uniform mesh (CW eq 1.6)

    for (int i=il-1; i<=(iu+1); ++i) {
      dph(i) = (7.0*(q(nin,k,j,i-1)+q(nin,k,j,i)) - (q(nin,k,j,i+1)+q(nin,k,j,i-2)))/12.0;
      d2qc(i) = q(nin,k,j,i-1) - 2.0*q(nin,k,j,i) + q(nin,k,j,i+1); //(CD eq 85a) (no 1/2)
    }
    d2qc(il-2) = q(nin,k,j,il-3) - 2.0*q(nin,k,j,il-2) + q(nin,k,j,il-1);

    // Limit interpolated interface states as in CD section 4.3.1
    for (int i=il-1; i<=(iu+1); ++i) {
      qa = dph(i) - q(nin,k,j,i-1); // (CD eq 84a)
      qb = q(nin,k,j,i) - dph(i);   // (CD eq 84b)
      if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
        qa = 3.0*(q(nin,k,j,i-1) - 2.0*dph(i) + q(nin,k,j,i));  // (CD eq 85b)
        qb = d2qc(i-1); // (CD eq 85a) (no 1/2)
        qc = d2qc(i);   // (CD eq 85c) (no 1/2)
        qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
          qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
        }
        dph(i) = 0.5*(q(nin,k,j,i-1)+q(nin,k,j,i)) - qd/6.0;
      } 
    }

    for (int i=il-1; i<=iu; ++i) {
      qminus(i) = dph(i  );
      qplus(i) =  dph(i+1 );
    }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)

    for (int i=il-1; i<=iu; ++i) {
      dqf_minus(i) = q(nin,k,j,i) - qminus(i); // (CS eq 25)
      dqf_plus(i)  = qplus(i) - q(nin,k,j,i);
      d2qf(i) = 6.0*(dph(i) - 2.0*q(nin,k,j,i) + dph(i+1)); // a6 coefficient * -2
    }

//--- Step 3. ----------------------------------------------------------------------------
// Apply CS limiters to parabolic interpolant

    for (int i=il-1; i<=iu; ++i) {
      qa = dqf_minus(i)*dqf_plus(i);
      qb = (q(nin,k,j,i+1) - q(nin,k,j,i))*(q(nin,k,j,i) - q(nin,k,j,i-1));

      // Check for local extrema
      if (qa <= 0.0 || qb <= 0.0 ) {
        // Check if extrema is smooth
        qa = d2qc(i-1);
        qb = d2qc(i);
        qc = d2qc(i+1);
        qd = d2qf(i);
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
          // Extrema is smooth
          qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                  std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
        } else {
          // Extrema is at interface adjacent to a discontinuity: flatten derivative
          qe =0.0;
        }

        // Check if 2nd derivative is close to roundoff error
        qa = std::max(fabs(q(nin,k,j,i-1)),fabs(q(nin,k,j,i-2)));
        qb = std::max(std::max(fabs(q(nin,k,j,i)),fabs(q(nin,k,j,i+1))),
                        fabs(q(nin,k,j,i+2)));

        if (fabs(qd) <= (1.0e-12)*std::max(qa,qb)) {
          // 2nd derivative magnitude is too small: flatten
          rho = 0.0;
        } else {
          // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
          rho = qe/qd; 
        }

        // Check if relative change in limited 2nd deriv is > roundoff
        if (rho <= (1.0 - (1.0e-12))) {
          // Limit smooth extrema
          qminus(i) = q(nin,k,j,i) - rho*dqf_minus(i); // (CS eq 23)
          qplus(i) = q(nin,k,j,i) + rho*dqf_plus(i);
        }

      // No extrema detected
      } else {
        // Overshoot i-1/2,R / i,(-) state
        if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
          qminus(i) = q(nin,k,j,i) - 2.0*dqf_plus(i);
        }
        // Overshoot i+1/2,L / i,(+) state
        if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
          qplus(i) = q(nin,k,j,i) + 2.0*dqf_minus(i);
        }
      }
    } 

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [il,iu]

    for (int i=il-1; i<=iu; ++i) {
      ql(nout,k,j,i+1) = qplus(i);
      qr(nout,k,j,i  ) = qminus(i);
    }

  }}

  return;
}

//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief Returns L/R interface values in X2-dir constructed using fourth-order PPM and
//         Colella-Sekora extremum-preserving limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PPMUniformX2(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25; 
  Real qa,qb,qc,qd,qe,rho;

  // 1D scratch arrays
  AthenaArray<Real> dph, dph_jp1, qplus, qminus, dqf_plus, dqf_minus, d2qf;
  AthenaArray<Real> d2qc_jm1, d2qc, d2qc_jp1;
  dph.InitWithShallowCopy(dph_);
  dph_jp1.InitWithShallowCopy(dph_p1_);
  qplus.InitWithShallowCopy(qplus_);
  qminus.InitWithShallowCopy(qminus_);
  dqf_plus.InitWithShallowCopy(dqf_plus_);
  dqf_plus.InitWithShallowCopy(dqf_minus_);
  d2qf.InitWithShallowCopy(d2qf_);
  d2qc_jm1.InitWithShallowCopy(d2qc_m1_);
  d2qc.InitWithShallowCopy(d2qc_);
  d2qc_jp1.InitWithShallowCopy(d2qc_p1_);

  for (int k=kl; k<=ku; ++k) {

//--- Step 1a. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} at j=jl-1 (CW eq 1.6)

    // initialize interface states along 1-D vector at j=jl-1
    for (int i=il; i<=iu; ++i) {
      dph(i) = ( 7.0*(q(nin,k,jl-2,i) + q(nin,k,jl-1,i)) - 
                     (q(nin,k,jl  ,i) + q(nin,k,jl-3,i)) )/12.0;
      // Approximate second-derivative at cell center using cell averages +/-
      d2qc_jm1(i) = q(nin,k,jl-3,i) - 2.0*q(nin,k,jl-2,i) + q(nin,k,jl-1,i);
      d2qc(i)     = q(nin,k,jl-2,i) - 2.0*q(nin,k,jl-1,i) + q(nin,k,jl  ,i);
    }

    // Limit interpolated interface states at j=jl-1 as in CD section 4.3.1
    for (int i=il; i<=iu; ++i) {
      qa = dph(i) - q(nin,k,jl-2,i); // (CD eq 84a)
      qb = q(nin,k,jl-1,i) - dph(i);   // (CD eq 84b)
      if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
        qa = 3.0*(q(nin,k,jl-2,i) - 2.0*dph(i) + q(nin,k,jl-1,i)); // (CD eq 85b)
        qb = d2qc_jm1(i); // (CD eq 85a)(no 1/2)
        qc = d2qc(i);     // (CD eq 85c)(no 1/2)
        qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
          qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
        }
        dph(i) = 0.5*(q(nin,k,jl-2,i)+q(nin,k,jl-1,i)) - qd/6.0;
      }
    }

// start loop over all j
    for (int j=jl-1; j<=ju; ++j) {

//--- Step 1b. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} at j=j+1 (CW eq 1.6)

      for (int i=il; i<=iu; ++i) {
        dph_jp1(i) = ( 7.0*(q(nin,k,j  ,i) + q(nin,k,j+1,i)) - 
                           (q(nin,k,j+2,i) + q(nin,k,j-1,i)) )/12.0;
        d2qc_jp1(i) = q(nin,k,j,i) - 2.0*q(nin,k,j+1,i) + q(nin,k,j+2,i);
      }

      // Limit interpolated interface states at j=j+1 as in CD section 4.3.1
      for (int i=il; i<=iu; ++i) {
        qa = dph_jp1(i) - q(nin,k,j,i); // (CD eq 84a)
        qb = q(nin,k,j+1,i) - dph_jp1(i); // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
          qa = 3.0*(q(nin,k,j,i) - 2.0*dph_jp1(i) + q(nin,k,j+1,i)); // (CD eq 85b)
          qb = d2qc(i);      // (CD eq 85a)(no 1/2)
          qc = d2qc_jp1(i);  // (CD eq 85c)(no 1/2)
          qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          dph_jp1(i) = 0.5*(q(nin,k,j,i)+q(nin,k,j+1,i)) - qd/6.0;
        }
      }

      // Initialize cell-indexed interface states / parabolic coefficients
      for (int i=il; i<=iu; ++i) {
        qminus(i) = dph(i);       // value at j
        qplus(i)  = dph_jp1(i);   // value at j+1
      }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)

      for (int i=il; i<=iu; ++i) {
        dqf_minus(i) = q(nin,k,j,i) - qminus(i);
        dqf_plus(i) = qplus(i) - q(nin,k,j,i);
        d2qf(i) = 6.0*(dph(i) -2.0*q(nin,k,j,i) + dph_jp1(i)); // a6 coefficient * -2
      }

//--- Step 3. ----------------------------------------------------------------------------
// Apply CS limiters to parabolic interpolant

      for (int i=il; i<=iu; ++i) {
        qa = dqf_minus(i)*dqf_plus(i);
        qb = (q(nin,k,j+1,i) - q(nin,k,j,i))*(q(nin,k,j,i) - q(nin,k,j-1,i));

        // check for local extrema
        if (qa <= 0.0 || qb <= 0.0 ) {
          // Check if extrema is smooth
          qa = d2qc_jm1(i);
          qb = d2qc(i);
          qc = d2qc_jp1(i);
          qd = d2qf(i);
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          } else {
            // Extrema is at interface adjacent to a discontinuity: flatten derivative
            qe = 0.0;
          }

          // Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q(nin,k,j-1,i)),fabs(q(nin,k,j-2,i)));
          qb = std::max(std::max(fabs(q(nin,k,j,i)),fabs(q(nin,k,j+1,i))),
                        fabs(q(nin,k,j+2,i)));

          if (fabs(qd) <= (1.0e-12)*std::max(qa,qb)) {
            // 2nd derivative magnitude is too small: flatten
            rho = 0.0;
          } else {
            // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
            rho = qe/qd;
          }

          // Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus(i) = q(nin,k,j,i) - rho*dqf_minus(i); // (CS eq 23)
            qplus(i)  = q(nin,k,j,i) + rho*dqf_plus(i);
          }

        // No extrema detected
        } else {
          // Overshoot j-1/2,R / j,(-) state
          if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
            qminus(i) = q(nin,k,j,i) - 2.0*dqf_plus(i);
          }
          // Overshoot j+1/2,L / j,(+) state
          if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
            qplus(i) = q(nin,k,j,i) + 2.0*dqf_minus(i);
          }
        } 
      }

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [jl,ju]

      for (int i=il; i<=iu; ++i) {
        ql(nout,k,j+1,i) = qplus(i);
        qr(nout,k,j  ,i) = qminus(i);
      } 

      // Copy 1D temporary arrays for next value of j unless j-loop finished
      if (j < ju) {
        for (int i=il; i<=iu; ++i) {
          dph(i) = dph_jp1(i);
          d2qc_jm1(i) = d2qc    (i);
          d2qc    (i) = d2qc_jp1(i);
        } 
      }

    } // end loop over [jl-1,ju]
  }   // end loop over [kl,ku]

  return;
}

//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//         Colella-Sekora extremum-preserving limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PPMUniformX3(Coordinates *pco, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  // CS08 constant used in second derivative limiter
  Real C2 = 1.25; // >1 , independent of h
  Real qa,qb,qc,qd,qe,rho;

  // 1D scratch arrays
  AthenaArray<Real> dph, dph_kp1, qplus, qminus, dqf_plus, dqf_minus, d2qf;
  AthenaArray<Real> d2qc_km1, d2qc, d2qc_kp1;
  dph.InitWithShallowCopy(dph_);
  dph_kp1.InitWithShallowCopy(dph_p1_);
  qplus.InitWithShallowCopy(qplus_);
  qminus.InitWithShallowCopy(qminus_);
  dqf_plus.InitWithShallowCopy(dqf_plus_);
  dqf_plus.InitWithShallowCopy(dqf_minus_);
  d2qf.InitWithShallowCopy(d2qf_);
  d2qc_km1.InitWithShallowCopy(d2qc_m1_);
  d2qc.InitWithShallowCopy(d2qc_);
  d2qc_kp1.InitWithShallowCopy(d2qc_p1_);

  for (int j=jl; j<=ju; ++j) {

//--- Step 1a. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} at k=kl-1 (CW eq 1.6)

    // initialize interface states along 1-D vector at k=kl-1
    for (int i=il; i<=iu; ++i) {
      dph(i) = ( 7.0*(q(nin,kl-2,j,i) + q(nin,kl-1,j,i)) -
                     (q(nin,kl  ,j,i) + q(nin,kl-3,j,i)) )/12.0;
      // Approximate second-derivative at cell center using cell averages +/-
      d2qc_km1(i) = q(nin,kl-3,j,i) - 2.0*q(nin,kl-2,j,i) + q(nin,kl-1,j,i);
      d2qc(i)     = q(nin,kl-2,j,i) - 2.0*q(nin,kl-1,j,i) + q(nin,kl  ,j,i);
    }

    // Limit interpolated interface states at k=kl-1 as in CD section 4.3.1
    for (int i=il; i<=iu; ++i) {
      qa = dph(i) - q(nin,kl-2,j,i);   // (CD eq 84a)
      qb = q(nin,kl-1,j,i) - dph(i);   // (CD eq 84b)
      if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
        qa = 3.0*(q(nin,kl-2,j,i) - 2.0*dph(i) + q(nin,kl-1,j,i)); // (CD eq 85b)
        qb = d2qc_km1(i); // (CD eq 85a)(no 1/2)
        qc = d2qc(i);     // (CD eq 85c)(no 1/2)
        qd = 0.0;
        if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
          qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
        }
        dph(i) = 0.5*(q(nin,kl-2,j,i)+q(nin,kl-1,j,i)) - qd/6.0;
      }
    }

// start loop over all k
    for (int k=kl-1; k<=ku; ++k) {

//--- Step 1b. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} at k=k+1 (CW eq 1.6)

      for (int i=il; i<=iu; ++i) {
        dph_kp1(i) = ( 7.0*(q(nin,k  ,j,i) + q(nin,k+1,j,i)) -
                           (q(nin,k+2,j,i) + q(nin,k-1,j,i)) )/12.0;
        d2qc_kp1(i) = q(nin,k,j,i) - 2.0*q(nin,k+1,j,i) + q(nin,k+2,j,i);
      }

      // Limit interpolated interface states at k=k+1 as in CD section 4.3.1
      for (int i=il; i<=iu; ++i) {
        qa = dph_kp1(i) - q(nin,k,j,i);   // (CD eq 84a)
        qb = q(nin,k+1,j,i) - dph_kp1(i); // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
          qa = 3.0*(q(nin,k,j,i) - 2.0*dph_kp1(i) + q(nin,k+1,j,i)); // (CD eq 85b)
          qb = d2qc(i);      // (CD eq 85a)(no 1/2)
          qc = d2qc_kp1(i);  // (CD eq 85c)(no 1/2)
          qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          dph_kp1(i) = 0.5*(q(nin,k,j,i)+q(nin,k+1,j,i)) - qd/6.0;
        }
      }

      // Initialize cell-indexed interface states / parabolic coefficients
      for (int i=il; i<=iu; ++i) {
        qminus(i) = dph(i);       // value at k
        qplus(i)  = dph_kp1(i);   // value at k+1
      }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)

      for (int i=il; i<=iu; ++i) {
        dqf_minus(i) = q(nin,k,j,i) - qminus(i);
        dqf_plus(i) = qplus(i) - q(nin,k,j,i);
        d2qf(i) = 6.0*(dph(i) -2.0*q(nin,k,j,i) + dph_kp1(i)); // a6 coefficient * -2
      }

//--- Step 3. ----------------------------------------------------------------------------
// Apply CS limiters to parabolic interpolant

      for (int i=il; i<=iu; ++i) {
        qa = dqf_minus(i)*dqf_plus(i);
        qb = (q(nin,k+1,j,i) - q(nin,k,j,i))*(q(nin,k,j,i) - q(nin,k-1,j,i));

        // check for local extrema
        if (qa <= 0.0 || qb <= 0.0 ) {
          // Check if extrema is smooth
          qa = d2qc_km1(i);
          qb = d2qc(i);
          qc = d2qc_kp1(i);
          qd = d2qf(i);
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          } else {
            // Extrema is at interface adjacent to a discontinuity: flatten derivative
            qe = 0.0;
          }

          // Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q(nin,k-1,j,i)),fabs(q(nin,k-2,j,i)));
          qb = std::max(std::max(fabs(q(nin,k,j,i)),fabs(q(nin,k+1,j,i))),
                        fabs(q(nin,k+2,j,i)));

          if (fabs(qd) <= (1.0e-12)*std::max(qa,qb)) {
            // 2nd derivative magnitude is too small: flatten
            rho = 0.0;
          } else {
            // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
            rho = qe/qd;
          }

          // Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus(i) = q(nin,k,j,i) - rho*dqf_minus(i); // (CS eq 23)
            qplus(i)  = q(nin,k,j,i) + rho*dqf_plus(i);
          }

        // No extrema detected
        } else {
          // Overshoot k-1/2,R / j,(-) state
          if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
            qminus(i) = q(nin,k,j,i) - 2.0*dqf_plus(i);
          }
          // Overshoot k+1/2,L / j,(+) state
          if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
            qplus(i) = q(nin,k,j,i) + 2.0*dqf_minus(i);
          }
        }
      }

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [jl,ju]

      for (int i=il; i<=iu; ++i) {
        ql(nout,k+1,j,i) = qplus(i);
        qr(nout,k  ,j,i) = qminus(i);
      }

      // Copy 1D temporary arrays for next value of k unless k-loop finished
      if (k < ku) {
        for (int i=il; i<=iu; ++i) {
          dph(i) = dph_kp1(i);
          d2qc_km1(i) = d2qc    (i);
          d2qc    (i) = d2qc_kp1(i);
        }
      }

    } // end loop over [kl-1,ku]
  }   // end loop over [jl,ju]

  return;
}
