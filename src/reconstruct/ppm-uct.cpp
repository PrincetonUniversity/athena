//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppm-uct.cpp
//  \brief piecewise parabolic primitive reconstruction with modified McCorquodale/Colella
//         limiter for uniform coordinates application of limited UCT
//
// REFERENCES:(MC) P. McCorquodale & P. Colella,  "A high-order finite-volume method for
// conservation laws on locally refined grids", CAMCoS, 6, 1 (2011)
//
// (CW)
// P. Colella & P. Woodward, "The Piecewise Parabolic Method (PPM) for Gas-Dynamical
// Simulations", JCP, 54, 174 (1984)
//
// (Mignone)
// A. Mignone, "High-order conservative reconstruction schemes for finite volume methods
// in cylindrical and spherical coordinates", JCP, 270, 784 (2014)
//
// (CS)
// P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at smooth extrema"
// , JCP, 227, 7069 (2008)
//
// (CD)
// P. Colella, M.R. Dorr, J. Hittinger, D. Martin, "High-order, finite-volume methods in
// mapped coordinates", JCP, 230, 2952 (2011)
//========================================================================================

// C++ headers
#include <algorithm>

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

void Reconstruction::PiecewiseParabolicUCTx1(MeshBlock *pmb, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  // CS08 constant used in second derivative limiter
  Real C2 = 1.25; // >1 , independent of h

  int ncells1 = iu - il +2*NGHOST;
  // 1D scratch array for nonlimited reconstructed interface values. Face-indexed j-1/2
  AthenaArray<Real> dph_;
  dph_.NewAthenaArray(ncells1);

  // 1D scratch arrays for limited reconstructed interface values. Cell-indexed j
  // plus/minus match +/- notation as in CS and CMHOG, and R/L notation as in CW
  AthenaArray<Real> qplus_, qminus_;
  qplus_.NewAthenaArray(ncells1); // Upper interface
  qminus_.NewAthenaArray(ncells1); // Lower interface

  // 1D arrays of finite difference approximations to derivatives
  AthenaArray<Real> dqf_plus_, dqf_minus_, d2qf_, d2qc_;
  dqf_plus_.NewAthenaArray(ncells1);
  dqf_minus_.NewAthenaArray(ncells1);
  d2qf_.NewAthenaArray(ncells1);
  d2qc_.NewAthenaArray(ncells1);

  // Assorted temporary scalars used in limiter
  Real qa,qb,qc,qd,qe,rho;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // Reconstruct interface averages <a>_{j-1/2}
      for (int i=il-1; i<=(iu+1); ++i) {
        // The first two ghost interfaces are required for limiting last real interface
        // Fourth-order reconstruction for uniform mesh (CW eq 1.6)
        // No van Leer monotonization of average zone slopes
        dph_(i) = 7.0/12.0*(q(nin,k,j,i-1) + q(nin,k,j,i))
          - 1.0/12.0*(q(nin,k,j,i+1) + q(nin,k,j,i-2));
      } // end reconstruction

      //-- Limit interpolated interface states as in CD section 4.3.1
      for (int i=il-1; i<=(iu+1); ++i) {
        // Either check monotonicity directly (CS eq 13) or equivalently,
        qa = dph_(i) - q(nin,k,j,i-1); // (CD eq 84a)
        qb = q(nin,k,j,i) - dph_(i); // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
          // centered appox. to second-derivative at i-1/2 (must recompute after limiting)
          qa = 3.0*(q(nin,k,j,i-1) - 2.0*dph_(i) + q(nin,k,j,i)); // (CD eq 85b)
          // Typo in CD! Use 1.0, not 0.5 coefficients
          // left appox. to second-derivative at i-1/2
          qb = 1.0*(q(nin,k,j,i-2) - 2.0*q(nin,k,j,i-1) + q(nin,k,j,i)); // (CD eq 85a)
          // right appox. to second-derivative at i-1/2
          qc = 1.0*(q(nin,k,j,i-1) - 2.0*q(nin,k,j,i) + q(nin,k,j,i+1)); // (CD eq 85c);
          // (the above two stencils are wastefully recomputed after this)

          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          } else {
            qd =0.0;
          }
          dph_(i) = 0.5*(q(nin,k,j,i-1)+q(nin,k,j,i)) - qd/6.0;
        } // end check of interpolated interface monotonicity
      }

      // Initialize cell-indexed interface states / parabolic coefficients
      for (int i=il-1; i<=iu; ++i) {
        qminus_(i) = dph_(i  );
        qplus_(i) =  dph_(i+1 );
      }

      //-- Compute cell-centered difference stencils (MC section 2.4.1)
      for (int i=il-1; i<=iu; ++i) {
        // Approximate first derivative *(h), alpha_{+/-} coefficients (CS eq 25)
        dqf_minus_(i) = q(nin,k,j,i) - qminus_(i);
        dqf_plus_(i) = qplus_(i) - q(nin,k,j,i);
        // Approximate second derivative at cell center *(h^2)
        // 1) Centered stencil using reconstructed interfaces +/- 1/2
        d2qf_(i) = 6.0*(dph_(i) -2.0*q(nin,k,j,i) + dph_(i+1)); // a6 coefficient * -2
        // 2) Centered stencil using cell averages +/- 1 quadrature points
        d2qc_(i) = q(nin,k,j,i-1) - 2.0*q(nin,k,j,i) + q(nin,k,j,i+1);
      }
      // Approximate additional second derivatives:
      d2qc_(il-2) = q(nin,k,j,il-3) - 2.0*q(nin,k,j,il-2) + q(nin,k,j,il-1);
      d2qc_(iu+1) = q(nin,k,j,iu) - 2.0*q(nin,k,j,iu+1) + q(nin,k,j,iu+2);

      //--------- Apply CS limiters to parabolic interpolant
      for (int i=il-1; i<=iu; ++i) {
        //--------- Check for local extrema
        // Condition #1: Cell average and interface average differences (CS eq 20)
        qa = dqf_minus_(i)*dqf_plus_(i);
        // Condition #2: Cell average and other cell average differences
        qb = (q(nin,k,j,i+1) - q(nin,k,j,i))*(q(nin,k,j,i) - q(nin,k,j,i-1));
        // Could use stencil width of +/-2 to reduce sensitivity to roundoff error

        //----------- Local extrema detected
        if (qa <= 0.0 || qb <= 0.0 ) { // If either face- or cell-avergage extrema
          //----------- Check if extrema is smooth
          qa = d2qc_(i-1);
          qb = d2qc_(i);
          qc = d2qc_(i+1);
          qd = d2qf_(i);
          // Limiting second derivative as a nonlinear combination of 4 approximations
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          } else {
            // Extrema is at interface adjacent to a discontinuity: flatten derivative
            qe =0.0;
          }
          //------------- Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q(nin,k,j,i-1)),fabs(q(nin,k,j,i-2)));
          qb = std::max(std::max(fabs(q(nin,k,j,i)),fabs(q(nin,k,j,i+1))),
                        fabs(q(nin,k,j,i+2)));
          //--------------- 2nd derivative magnitude is too small: flatten
          // and/or dont divide by zero
          if (fabs(qd) <= (1.0e-12)*std::max(qa,qb))
            rho = 0.0;
          //--------------- Limiter is not sensitive to roundoff. Use limited ratio
          else
            rho = qe/qd; // (MC eq 27)

          //------------- Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus_(i) = q(nin,k,j,i) - rho*dqf_minus_(i); // (CS eq 23)
            qplus_(i) = q(nin,k,j,i) + rho*dqf_plus_(i);
          } // end loop over rho roundoff sensitivity
        } else { // end check of local extrema 2x conditions
          //----------- No extrema detected
          // Construct the parabolic interpolant with a monotonicity preserving limiter
          // Using standard CW overshoot limiters; see Martin memo for diffusivity info

          // Overshoot j-1/2,R / j,(-) state
          if (fabs(dqf_minus_(i)) >= 2.0*fabs(dqf_plus_(i)))
            qminus_(i) = q(nin,k,j,i) - 2.0*dqf_plus_(i);
          // Overshoot j+1/2,L / j,(+) state
          if (fabs(dqf_plus_(i)) >= 2.0*fabs(dqf_minus_(i)))
            qplus_(i) = q(nin,k,j,i) + 2.0*dqf_minus_(i);
        } // end else statement (not a local extrema)
      } // end loop over i=il-1:iu, applying Colella/Sekora limiters

      // Convert limited cell-centered values to interface-centered L/R Riemann states
      for (int i=il; i<=iu; ++i) {
        ql(nout,k,j,i) = qplus_(i-1);
        qr(nout,k,j,i) = qminus_(i);
      }
    } // end loop over j
  }// end loop over k

  // Free 1D scratch arrays
  qplus_.DeleteAthenaArray();
  qminus_.DeleteAthenaArray();
  dph_.DeleteAthenaArray();

  dqf_plus_.DeleteAthenaArray();
  dqf_minus_.DeleteAthenaArray();
  d2qf_.DeleteAthenaArray();
  d2qc_.DeleteAthenaArray();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

void Reconstruction::PiecewiseParabolicUCTx2(MeshBlock *pmb, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  // CS08 constant used in second derivative limiter
  Real C2 = 1.25; // >1 , independent of h

  int ncells1 = iu - il +2*NGHOST;
  // 1D scratch arrays for nonlimited reconstructed interface values. Face-indexed j-1/2
  AthenaArray<Real> dph_x2_;
  AthenaArray<Real> dph_, dph_jp1_;
  dph_x2_.NewAthenaArray(2,ncells1);
  dph_.InitWithShallowSlice(dph_x2_, 2, 0, 0);
  dph_jp1_.InitWithShallowSlice(dph_x2_, 2, 1, 0);

  // 1D scratch arrays for limited reconstructed interface values. Cell-indexed j
  // plus/minus match +/- notation as in CS and CMHOG, and R/L notation as in CW
  AthenaArray<Real> qplus_, qminus_;
  qplus_.NewAthenaArray(ncells1); // Upper interface
  qminus_.NewAthenaArray(ncells1); // Lower interface

  // 1D arrays of finite difference approximations to derivatives
  AthenaArray<Real> dqf_plus_, dqf_minus_, d2qf_;
  dqf_plus_.NewAthenaArray(ncells1);
  dqf_minus_.NewAthenaArray(ncells1);
  d2qf_.NewAthenaArray(ncells1);
  AthenaArray<Real> d2qc_x2_;
  AthenaArray<Real> d2qc_jm1_, d2qc_, d2qc_jp1_;
  d2qc_x2_.NewAthenaArray(3,ncells1);
  d2qc_jm1_.InitWithShallowSlice(d2qc_x2_, 2, 0, 0);
  d2qc_.InitWithShallowSlice(d2qc_x2_, 2, 1, 0);
  d2qc_jp1_.InitWithShallowSlice(d2qc_x2_, 2, 2, 0);
  // Temporary scratch pointer for swapping pointers in x2
  AthenaArray<Real> x2scratch_;

  // Assorted temporary scalars used in limiter
  Real qa,qb,qc,qd,qe,rho;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl-1; j<=ju; ++j) {
      // Startup the x1 sweeps on lowest j
      if (j==jl-1) {
        for (int i=il; i<=iu; ++i) {
          // Reconstruct interface average for lowest x2 cells
          dph_(i) = 7.0/12.0*(q(nin,k,j-1,i) + q(nin,k,j,i))
            - 1.0/12.0*(q(nin,k,j+1,i) + q(nin,k,j-2,i));
          // Approximate second-derivative at cell center using cell averages +/-
          d2qc_jm1_(i) = q(nin,k,j-2,i) - 2.0*q(nin,k,j-1,i) + q(nin,k,j,i);
          d2qc_(i) = q(nin,k,j-1,i) - 2.0*q(nin,k,j,i) + q(nin,k,j+1,i);
          d2qc_jp1_(i) = q(nin,k,j,i) - 2.0*q(nin,k,j+1,i) + q(nin,k,j+2,i);

          //-- Limit interpolated interface states as in CD section 4.3.1
          // Either check monotonicity directly (CS eq 13) or equivalently,
          qa = dph_(i) - q(nin,k,j-1,i); // (CD eq 84a)
          qb = q(nin,k,j,i) - dph_(i); // (CD eq 84b)
          if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
            // centered appox. to second-derivative at i-1/2  (must recompute after limit)
            qa = 3.0*(q(nin,k,j-1,i) - 2.0*dph_(i) + q(nin,k,j,i)); // (CD eq 85b)
            // Typo in CD! Use 1.0, not 0.5 coefficients
            // left appox. to second-derivative at i-1/2
            qb = 1.0*(q(nin,k,j-2,i) - 2.0*q(nin,k,j-1,i) + q(nin,k,j,i)); // (CD eq 85a)
            // right appox. to second-derivative at i-1/2
            qc = 1.0*(q(nin,k,j-1,i) - 2.0*q(nin,k,j,i) + q(nin,k,j+1,i)); // (CD eq 85c)
            // (the above two stencils are wastefully recomputed)

            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
              qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
            } else {
              qd =0.0;
            }
            dph_(i) = 0.5*(q(nin,k,j-1,i)+q(nin,k,j,i)) - qd/6.0;
          }  // end check of interpolated interface monotonicity
        }
      } else {
        // Reuse old x1 slice information from previous j
        //-------- Swap pointers of 1D scratch arrays from previous x1 slice sweep
         // Non-limited interface averages at j-1/2, j+1/2
        x2scratch_.InitWithShallowCopy(dph_);
        dph_.InitWithShallowCopy(dph_jp1_);
        dph_jp1_.InitWithShallowCopy(x2scratch_);
        // Shift (downward) 1D array pointers, dequeue lowest one
        x2scratch_.InitWithShallowCopy(d2qc_jm1_);
        d2qc_jm1_.InitWithShallowCopy(d2qc_);
        d2qc_.InitWithShallowCopy(d2qc_jp1_);
        d2qc_jp1_.InitWithShallowCopy(x2scratch_);
      }

      // Reconstruct x2 interfaces along x1 slices for performance
      for (int i=il; i<=iu; ++i) {
        // Reconstruct interface averages <a>_{j+1/2}
        // Fourth-order reconstruction for uniform mesh (CW eq 1.6)
        // No van Leer monotonization of average zone slopes
        dph_jp1_(i) = 7.0/12.0*(q(nin,k,j,i) + q(nin,k,j+1,i))
          - 1.0/12.0*(q(nin,k,j+2,i) + q(nin,k,j-1,i));

        //-- Limit interpolated interface states as in CD section 4.3.1
        // Either check monotonicity directly (CS eq 13) or equivalently,
        qa = dph_jp1_(i) - q(nin,k,j,i); // (CD eq 84a)
        qb = q(nin,k,j+1,i) - dph_jp1_(i); // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at j+1/2 face
          // centered appox. to second-derivative at j+1/2  (must recompute after limit)
          qa = 3.0*(q(nin,k,j,i) - 2.0*dph_jp1_(i) + q(nin,k,j+1,i)); // (CD eq 85b)
          // Typo in CD! Use 1.0, not 0.5 coefficients
          // left appox. to second-derivative at j+1/2
          qb = 1.0*(q(nin,k,j-1,i) - 2.0*q(nin,k,j,i) + q(nin,k,j+1,i)); // (CD eq 85a)
          // right appox. to second-derivative at j+1/2
          qc = 1.0*(q(nin,k,j,i) - 2.0*q(nin,k,j+1,i) + q(nin,k,j+2,i)); // (CD eq 85c)
          // (the above two stencils are wastefully recomputed)

          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          } else {
            qd =0.0;
          }
          dph_jp1_(i) = 0.5*(q(nin,k,j+1,i)+q(nin,k,j,i)) - qd/6.0;
        } // end check of interpolated interface monotonicity

        //-- Compute cell-centered difference stencils (MC section 2.4.1)
        // Initialize cell-indexed interface states / parabolic coefficients
        qminus_(i) = dph_(i);
        qplus_(i)  =  dph_jp1_(i);

        // Approximate first derivative *(h), alpha_{+/-} coefficients (CS eq 25)
        dqf_minus_(i) = q(nin,k,j,i) - qminus_(i);
        dqf_plus_(i) = qplus_(i) - q(nin,k,j,i);
        // Approximate second derivative at cell center *(h^2)
        // 1) Centered stencil using reconstructed interfaces +/- 1/2
        d2qf_(i) = 6.0*(dph_(i) -2.0*q(nin,k,j,i) + dph_jp1_(i)); // a6 coefficient * -2
        // 2) Centered stencil using cell averages +/- 1
        d2qc_jp1_(i) = q(nin,k,j,i) - 2.0*q(nin,k,j+1,i) + q(nin,k,j+2,i);

        //-- Compute interface-centered difference stencils

        //--------- Apply CS limiters to parabolic interpolant
        //----------- Check for local extrema
        // Condition #1: Cell average and interface average differences (CS eq 20)
        qa = dqf_minus_(i)*dqf_plus_(i);
        // Condition #2: Cell average and other cell average differences
        qb = (q(nin,k,j+1,i) - q(nin,k,j,i))*(q(nin,k,j,i) - q(nin,k,j-1,i));
        // Could use stencil width of +/-2 to reduce sensitivity to roundoff error

        //----------- Local extrema detected
        if (qa <= 0.0 || qb <= 0.0 ) { // If either face- or cell-avergage extrema
          //----------- Check if extrema is smooth
          qa = d2qc_jm1_(i);
          qb = d2qc_(i);
          qc = d2qc_jp1_(i);
          qd = d2qf_(i);
          // Limiting second derivative as a nonlinear combination of 4 approximations
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is possibly smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          } else {
            // Extrema is at interface adjacent to a discontinuity: flatten derivative
            qe =0.0;
          }
          //------------- Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q(nin,k,j-1,i)),fabs(q(nin,k,j-2,i)));
          qb = std::max(std::max(fabs(q(nin,k,j,i)),fabs(q(nin,k,j+1,i))),
                        fabs(q(nin,k,j+2,i)));
          //--------------- 2nd derivative magnitude is too small: flatten
          // and/or dont divide by zero
          if (fabs(qd) <= (1.0e-12)*std::max(qa,qb))
            rho = 0.0;
          //--------------- Limiter is not sensitive to roundoff. Use limited ratio
          else
            rho = qe/qd; // (MC eq 27)

          //------------- Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus_(i) = q(nin,k,j,i) - rho*dqf_minus_(i); // (CS eq 23)
            qplus_(i) = q(nin,k,j,i) + rho*dqf_plus_(i);
          } // end loop over rho roundoff sensitivity
        } else { // end check of local extrema 2x conditions
          //----------- No extrema detected:
          // Construct the parabolic interpolant with a monotonicity preserving limiter
          // Using standard CW overshoot limiters; see Martin memo for diffusivity info

          // Overshoot j-1/2,R / j,(-) state
          if (fabs(dqf_minus_(i)) >= 2.0*fabs(dqf_plus_(i)))
            qminus_(i) = q(nin,k,j,i) - 2.0*dqf_plus_(i);
          // Overshoot j+1/2,L / j,(+) state
          if (fabs(dqf_plus_(i)) >= 2.0*fabs(dqf_minus_(i)))
            qplus_(i) = q(nin,k,j,i) + 2.0*dqf_minus_(i);
        } // end else statement (not a local extrema)

        // Convert limited cell-centered values to interface-centered L/R Riemann states
        ql(nout,k,j+1,i) = qplus_(i);
        qr(nout,k,j,i) = qminus_(i);
      } // end loop over i=il:iu, applying McCorquodale limiters
    } // end loop over j
  }// end loop over k

  // Free 1D scratch arrays
  // Allocated AthenaArray
  qplus_.DeleteAthenaArray();
  qminus_.DeleteAthenaArray();
  dqf_plus_.DeleteAthenaArray();
  dqf_minus_.DeleteAthenaArray();
  d2qc_x2_.DeleteAthenaArray();
  dph_x2_.DeleteAthenaArray();

  // Shallow copied AthenaArray
  dph_.DeleteAthenaArray();
  dph_jp1_.DeleteAthenaArray();

  d2qf_.DeleteAthenaArray();

  d2qc_jm1_.DeleteAthenaArray();
  d2qc_.DeleteAthenaArray();
  d2qc_jp1_.DeleteAthenaArray();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

void Reconstruction::PiecewiseParabolicUCTx3(MeshBlock *pmb, const int kl, const int ku,
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  // CS08 constant used in second derivative limiter
  Real C2 = 1.25; // >1 , independent of h

  int ncells1 = iu - il +2*NGHOST;
  // 1D scratch arrays for nonlimited reconstructed interface values. Face-indexed j-1/2
  AthenaArray<Real> dph_x3_;
  AthenaArray<Real> dph_, dph_kp1_;
  dph_x3_.NewAthenaArray(2,ncells1);
  dph_.InitWithShallowSlice(dph_x3_, 2, 0, 0);
  dph_kp1_.InitWithShallowSlice(dph_x3_, 2, 1, 0);

  // 1D scratch arrays for limited reconstructed interface values. Cell-indexed k
  // plus/minus match +/- notation as in CS and CMHOG, and R/L notation as in CW
  AthenaArray<Real> qplus_, qminus_;
  qplus_.NewAthenaArray(ncells1); // Upper interface
  qminus_.NewAthenaArray(ncells1); // Lower interface

  // 1D arrays of finite difference approximations to derivatives
  AthenaArray<Real> dqf_plus_, dqf_minus_, d2qf_;
  dqf_plus_.NewAthenaArray(ncells1);
  dqf_minus_.NewAthenaArray(ncells1);
  d2qf_.NewAthenaArray(ncells1);
  AthenaArray<Real> d2qc_x3_;
  AthenaArray<Real> d2qc_km1_, d2qc_, d2qc_kp1_;
  d2qc_x3_.NewAthenaArray(3,ncells1);
  d2qc_km1_.InitWithShallowSlice(d2qc_x3_, 2, 0, 0);
  d2qc_.InitWithShallowSlice(d2qc_x3_, 2, 1, 0);
  d2qc_kp1_.InitWithShallowSlice(d2qc_x3_, 2, 2, 0);
  // Temporary scratch pointer for swapping pointers in x3
  AthenaArray<Real> x3scratch_;

  // Assorted temporary scalars used in limiter
  Real qa,qb,qc,qd,qe,rho;

  for (int j=jl; j<=ju; ++j) {
    for (int k=kl-1; k<=ku; ++k) {
      // Startup the x1 sweeps on lowest k
      if (k==kl-1) {
        for (int i=il; i<=iu; ++i) {
          // Reconstruct interface average for lowest x3 cells
          dph_(i) = 7.0/12.0*(q(nin,k-1,j,i) + q(nin,k,j,i))
            - 1.0/12.0*(q(nin,k+1,j,i) + q(nin,k-2,j,i));
          // Approximate second-derivative at cell center using cell averages +/-
          d2qc_km1_(i) = q(nin,k-2,j,i) - 2.0*q(nin,k-1,j,i) + q(nin,k,j,i);
          d2qc_(i) = q(nin,k-1,j,i) - 2.0*q(nin,k,j,i) + q(nin,k+1,j,i);
          d2qc_kp1_(i) = q(nin,k,j,i) - 2.0*q(nin,k+1,j,i) + q(nin,k+2,j,i);

          //-- Limit interpolated interface states as in CD section 4.3.1
          // Either check monotonicity directly (CS eq 13) or equivalently,
          qa = dph_(i) - q(nin,k-1,j,i); // (CD eq 84a)
          qb = q(nin,k,j,i) - dph_(i); // (CD eq 84b)
          if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
            // centered appox. to second-derivative at i-1/2  (must recompute after limit)
            qa = 3.0*(q(nin,k-1,j,i) - 2.0*dph_(i) + q(nin,k,j,i)); // (CD eq 85b)
            // Typo in CD! Use 1.0, not 0.5 coefficients
            // left appox. to second-derivative at i-1/2
            qb = 1.0*(q(nin,k-2,j,i) - 2.0*q(nin,k-1,j,i) + q(nin,k,j,i)); // (CD eq 85a)
            // right appox. to second-derivative at i-1/2
            qc = 1.0*(q(nin,k-1,j,i) - 2.0*q(nin,k,j,i) + q(nin,k+1,j,i)); // (CD eq 85c)
            // (the above two stencils are wastefully recomputed)

            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
              qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
            } else {
              qd = 0.0;
            }
            dph_(i) = 0.5*(q(nin,k-1,j,i)+q(nin,k,j,i)) - qd/6.0;
          }  // end check of interpolated interface monotonicity
        }
      } else { // Reuse old x1 slice information from previous j:
        //-------- Swap pointers of 1D scratch arrays from previous x1 slice sweep
         // Non-limited interface averages at j-1/2, j+1/2
        x3scratch_.InitWithShallowCopy(dph_);
        dph_.InitWithShallowCopy(dph_kp1_);
        dph_kp1_.InitWithShallowCopy(x3scratch_);
        // Shift (downward) 1D array pointers, dequeue lowest one
        x3scratch_.InitWithShallowCopy(d2qc_km1_);
        d2qc_km1_.InitWithShallowCopy(d2qc_);
        d2qc_.InitWithShallowCopy(d2qc_kp1_);
        d2qc_kp1_.InitWithShallowCopy(x3scratch_);
      }

      // Reconstruct x3 interfaces along x1 slices for performance
      for (int i=il; i<=iu; ++i) {
        // Reconstruct interface averages <a>_{j+1/2}
        // Fourth-order reconstruction for uniform mesh (CW eq 1.6)
        // No van Leer monotonization of average zone slopes
        dph_kp1_(i) = 7.0/12.0*(q(nin,k,j,i) + q(nin,k+1,j,i))
          - 1.0/12.0*(q(nin,k+2,j,i) + q(nin,k-1,j,i));

        //-- Limit interpolated interface states as in CD section 4.3.1
        // Either check monotonicity directly (CS eq 13) or equivalently,
        qa = dph_kp1_(i) - q(nin,k,j,i); // (CD eq 84a)
        qb = q(nin,k+1,j,i) - dph_kp1_(i); // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at k+1/2 face
          // centered appox. to second-derivative at k+1/2  (must recompute after limit)
          qa = 3.0*(q(nin,k,j,i) - 2.0*dph_kp1_(i) + q(nin,k+1,j,i)); // (CD eq 85b)
          // Typo in CD! Use 1.0, not 0.5 coefficients
          // left appox. to second-derivative at k+1/2
          qb = 1.0*(q(nin,k-1,j,i) - 2.0*q(nin,k,j,i) + q(nin,k+1,j,i)); // (CD eq 85a)
          // right appox. to second-derivative at k+1/2
          qc = 1.0*(q(nin,k,j,i) - 2.0*q(nin,k+1,j,i) + q(nin,k+2,j,i)); // (CD eq 85c)
          // (the above two stencils are wastefully recomputed)

          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)) {
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }     else {
            qd =0.0;
          }
          dph_kp1_(i) = 0.5*(q(nin,k+1,j,i)+q(nin,k,j,i)) - qd/6.0;
        }  // end check of interpolated interface monotonicity

        //-- Compute cell-centered difference stencils (MC section 2.4.1)
        // Initialize cell-indexed interface states / parabolic coefficients
        qminus_(i) = dph_(i);
        qplus_(i)  =  dph_kp1_(i);

        // Approximate first derivative *(h), alpha_{+/-} coefficients (CS eq 25)
        dqf_minus_(i) = q(nin,k,j,i) - qminus_(i);
        dqf_plus_(i) = qplus_(i) - q(nin,k,j,i);
        // Approximate second derivative at cell center *(h^2)
        // 1) Centered stencil using reconstructed interfaces +/- 1/2
        d2qf_(i) = 6.0*(dph_(i) -2.0*q(nin,k,j,i) + dph_kp1_(i));
        // 2) Centered stencil using cell averages +/- 1 // a6 coefficient * -2
        d2qc_kp1_(i) = q(nin,k,j,i) - 2.0*q(nin,k+1,j,i) + q(nin,k+2,j,i);

        //-- Compute interface-centered difference stencils

        //--------- Apply CS limiters to parabolic interpolant
        //----------- Check for local extrema
        // Condition #1: Cell average and interface average differences (CS eq 20)
        qa = dqf_minus_(i)*dqf_plus_(i);
        // Condition #2: Cell average and other cell average differences
        qb = (q(nin,k+1,j,i) - q(nin,k,j,i))*(q(nin,k,j,i) - q(nin,k-1,j,i));
        // Could use stencil width of +/-2 to reduce sensitivity to roundoff error

        //----------- Local extrema detected
        if (qa <= 0.0 || qb <= 0.0 ) { // If either face- or cell-avergage extrema
          //----------- Check if extrema is smooth
          qa = d2qc_km1_(i);
          qb = d2qc_(i);
          qc = d2qc_kp1_(i);
          qd = d2qf_(i);
          // Limiting second derivative as a nonlinear combination of 4 approximations
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is possibly smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          } else {
            // Extrema is at interface adjacent to a discontinuity: flatten derivative
            qe =0.0;
          }
          //------------- Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q(nin,k-1,j,i)),fabs(q(nin,k-2,j,i)));
          qb = std::max(std::max(fabs(q(nin,k,j,i)),fabs(q(nin,k+1,j,i))),
                        fabs(q(nin,k+2,j,i)));
          //--------------- 2nd derivative magnitude is too small: flatten
          // and/or dont divide by zero
          if (fabs(qd) <= (1.0e-12)*std::max(qa,qb))
            rho = 0.0;
          //--------------- Limiter is not sensitive to roundoff. Use limited ratio
          else
            rho = qe/qd; // (MC eq 27)

          //------------- Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus_(i) = q(nin,k,j,i) - rho*dqf_minus_(i); // (CS eq 23)
            qplus_(i) = q(nin,k,j,i) + rho*dqf_plus_(i);
          } // end loop over rho roundoff sensitivity
        } else { // end check of local extrema 2x conditions
          //----------- No extrema detected
          // Construct the parabolic interpolant with a monotonicity preserving limiter
          // Using standard CW overshoot limiters; see Martin memo for diffusivity info

          // Overshoot j-1/2,R / j,(-) state
          if (fabs(dqf_minus_(i)) >= 2.0*fabs(dqf_plus_(i)))
            qminus_(i) = q(nin,k,j,i) - 2.0*dqf_plus_(i);
          // Overshoot j+1/2,L / j,(+) state
          if (fabs(dqf_plus_(i)) >= 2.0*fabs(dqf_minus_(i)))
            qplus_(i) = q(nin,k,j,i) + 2.0*dqf_minus_(i);
        } // end else statement (not a local extrema)

        // Convert limited cell-centered values to interface-centered L/R Riemann states
        ql(nout,k+1,j,i) = qplus_(i);
        qr(nout,k,j,i) = qminus_(i);
      } // end loop over i=il:iu, applying McCorquodale limiters
    } // end loop over j
  }// end loop over k

  // Free 1D scratch arrays
  // Allocated AthenaArray
  qplus_.DeleteAthenaArray();
  qminus_.DeleteAthenaArray();
  dqf_plus_.DeleteAthenaArray();
  dqf_minus_.DeleteAthenaArray();
  d2qc_x3_.DeleteAthenaArray();
  dph_x3_.DeleteAthenaArray();

  // Shallow copied AthenaArray
  dph_.DeleteAthenaArray();
  dph_kp1_.DeleteAthenaArray();

  d2qf_.DeleteAthenaArray();

  d2qc_km1_.DeleteAthenaArray();
  d2qc_.DeleteAthenaArray();
  d2qc_kp1_.DeleteAthenaArray();

  return;
}
