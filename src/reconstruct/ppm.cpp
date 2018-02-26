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

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX1()
//  \brief Returns L/R interface values in X1-dir constructed using fourth-order PPM and
//         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX1(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
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

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
    // cache the x1-sliced primitive states for eigensystem calculation
    for (int n=0; n<(NHYDRO); ++n) {
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        wc(n,i) = w(n,k,j,i);
        q    (n,i) = w(n,k,j,i  );
        q_im2(n,i) = w(n,k,j,i-2);
        q_im1(n,i) = w(n,k,j,i-1);
        q_ip1(n,i) = w(n,k,j,i+1);
        q_ip2(n,i) = w(n,k,j,i+2);
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        bx(i) = bcc(IB1,k,j,i);
      }
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        wc(IBY,i) = bcc(IB2,k,j,i);
        q    (IBY,i) = bcc(IB2,k,j,i  );
        q_im2(IBY,i) = bcc(IB2,k,j,i-2);
        q_im1(IBY,i) = bcc(IB2,k,j,i-1);
        q_ip1(IBY,i) = bcc(IB2,k,j,i+1);
        q_ip2(IBY,i) = bcc(IB2,k,j,i+2);
      }
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        wc(IBZ,i) = bcc(IB3,k,j,i);
        q    (IBZ,i) = bcc(IB3,k,j,i  );
        q_im2(IBZ,i) = bcc(IB3,k,j,i-2);
        q_im1(IBZ,i) = bcc(IB3,k,j,i-1);
        q_ip1(IBZ,i) = bcc(IB3,k,j,i+1);
        q_ip2(IBZ,i) = bcc(IB3,k,j,i+2);
      }
    }

    // Project cell-averages to characteristic variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      // equivalent to looping from il-3 to iu+2 for primitive reconstruction
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
      for (int i=il-1; i<=iu; ++i){
        Real qa = (q(n,i) - q_im1(n,i));
        Real qb = (q_ip1(n,i) - q(n,i));
        dd_im1(i) = prec->c1i(i-1)*qa + prec->c2i(i-1)*(q_im1(n,i) - q_im2(n,i));
        dd    (i) = prec->c1i(i  )*qb + prec->c2i(i  )*qa;
        dd_ip1(i) = prec->c1i(i+1)*(q_ip2(n,i) - q_ip1(n,i)) + prec->c2i(i+1)*qb;
      }
      // Approximate interface average at i-1/2 and i+1/2 using PPM (CW eq 1.6)
      for (int i=il-1; i<=iu; ++i) {
        dph(i)= prec->c3i(i)*q_im1(n,i) + prec->c4i(i)*q(n,i) +
                prec->c5i(i)*dd_im1(i) + prec->c6i(i)*dd(i);
        dph_ip1(i)= prec->c3i(i+1)*q(n,i) + prec->c4i(i+1)*q_ip1(n,i) +
                    prec->c5i(i+1)*dd(i) + prec->c6i(i+1)*dd_ip1(i);
      }

//--- Step 2a. ---------------------------------------------------------------------------
      // For a uniform grid, limit interpolated interface states as in CD section 4.3.1
      if (pmb->block_size.x1rat == 1.0) {
        // approximate second derivative at interfaces for smooth extrema preservation
#pragma simd
        for (int i=il-1; i<=iu+1; ++i) {
          d2qc_im1(i) = q_im2(n,i) - 2.0*q_im1(n,i) + q    (n,i);
          d2qc    (i) = q_im1(n,i) - 2.0*q    (n,i) + q_ip1(n,i); //(CD eq 85a) (no 1/2)
          d2qc_ip1(i) = q    (n,i) - 2.0*q_ip1(n,i) + q_ip2(n,i);
        }

        // i-1/2
        // #pragma simd // poor vectorization efficiency
        for (int i=il-1; i<=(iu+1); ++i) {
          Real qa = dph(i) - q_im1(n,i); // (CD eq 84a)
          Real qb = q(n,i) - dph(i);     // (CD eq 84b)
          if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
            qa = 3.0*(q_im1(n,i) - 2.0*dph(i) + q(n,i));  // (CD eq 85b)
            qb = d2qc_im1(i);    // (CD eq 85a) (no 1/2)
            Real qc = d2qc(i);   // (CD eq 85c) (no 1/2)
            Real qd = 0.0;
            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
              qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
            }
            dph(i) = 0.5*(q_im1(n,i)+q(n,i)) - qd/6.0;
          }
        }
        // i+1/2
        // #pragma simd // poor vectorization efficiency
        for (int i=il-1; i<=(iu+1); ++i) {
          Real qa = dph_ip1(i) - q(n,i);       // (CD eq 84a)
          Real qb = q_ip1(n,i) - dph_ip1(i);   // (CD eq 84b)
          if (qa*qb < 0.0) { // Local extrema detected at i+1/2 face
            qa = 3.0*(q(n,i) - 2.0*dph_ip1(i) + q_ip1(n,i));  // (CD eq 85b)
            qb = d2qc(i);            // (CD eq 85a) (no 1/2)
            Real qc = d2qc_ip1(i);   // (CD eq 85c) (no 1/2)
            Real qd = 0.0;
            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
              qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
            }
            dph_ip1(i) = 0.5*(q(n,i)+q_ip1(n,i)) - qd/6.0;
          }
        }

#pragma simd
        for (int i=il-1; i<=iu; ++i) {
          d2qf(i) = 6.0*(dph(i) - 2.0*q(n,i) + dph_ip1(i)); // a6 coefficient * -2
        }

//--- Step 2b. ---------------------------------------------------------------------------
      // For a non-uniform grid, apply strict monotonicity constraints (Mignone eq 45)
      } else {
        for (int i=il-1; i<=iu; ++i) {
          dph    (i) = std::min(dph    (i), std::max(q(n,i),q_im1(n,i)));
          dph_ip1(i) = std::min(dph_ip1(i), std::max(q(n,i),q_ip1(n,i)));

          dph    (i) = std::max(dph    (i), std::min(q(n,i),q_im1(n,i)));
          dph_ip1(i) = std::max(dph_ip1(i), std::min(q(n,i),q_ip1(n,i)));
        }
      }

      // Cache Riemann states for both non-/uniform limiters
#pragma simd
      for (int i=il-1; i<=iu; ++i) {
        qminus(i) = dph(i  );
        qplus(i) =  dph_ip1(i );
      }

//--- Step 3. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
#pragma simd
      for (int i=il-1; i<=iu; ++i) {
        dqf_minus(i) = q(n,i) - qminus(i); // (CS eq 25)
        dqf_plus(i)  = qplus(i) - q(n,i);
      }

//--- Step 4a. ---------------------------------------------------------------------------
      // For uniform mesh: apply CS limiters to parabolic interpolant
      if (pmb->block_size.x1rat == 1.0) {
        // #pragma simd // poor vectorization efficiency
        for (int i=il-1; i<=iu; ++i) {
          Real qa = dqf_minus(i)*dqf_plus(i);
          Real qb = (q_ip1(n,i) - q(n,i))*(q(n,i) - q_im1(n,i));

          // Check for local extrema
          if (qa <= 0.0 || qb <= 0.0 ) {
            // Check if extrema is smooth
            qa = d2qc_im1(i);
            qb = d2qc(i);
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

            // Check if relative change in limited 2nd deriv is > roundoff
            if (rho <= (1.0 - (1.0e-12))) {
              // Limit smooth extrema
              qminus(i) = q(n,i) - rho*dqf_minus(i); // (CS eq 23)
              qplus(i) = q(n,i) + rho*dqf_plus(i);
            }

          // No extrema detected
          } else {
            // Overshoot i-1/2,R / i,(-) state
            if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
              qminus(i) = q(n,i) - 2.0*dqf_plus(i);
            }
            // Overshoot i+1/2,L / i,(+) state
            if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
              qplus(i) = q(n,i) + 2.0*dqf_minus(i);
            }
          }
        }

//--- Step 4b. ---------------------------------------------------------------------------
      // For non-uniform mesh: apply Mignone limiters to parabolic interpolant
      // Note Mignone limiter does not check for cell-averaged extrema:
      } else {
        for (int i=il-1; i<=iu; ++i) {
          Real qa = dqf_minus(i)*dqf_plus(i);
          if (qa <= 0.0) { // Local extrema detected
            qminus(i) = q(n,i);
            qplus(i) = q(n,i);
          } else { // No extrema detected
            // Overshoot i-1/2,R / i,(-) state
            if (fabs(dqf_minus(i)) >= 2.0*fabs(dqf_plus(i))) {
              qminus(i) = q(n,i) - 2.0*dqf_plus(i);
            }
            // Overshoot i+1/2,L / i,(+) state
            if (fabs(dqf_plus(i)) >= 2.0*fabs(dqf_minus(i))) {
              qplus(i) = q(n,i) + 2.0*dqf_minus(i);
            }
          }
        }
      }

//--- Step 5. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [il,iu]
#pragma simd
      for (int i=il-1; i<=iu; ++i) {
        ql_iph(n,i ) = qplus(i);
        qr_imh(n,i ) = qminus(i);
      }
    } // end char PPM loop over NWAVE

    // Project limited slope back to primitive variables, if necessary
    if (pmb->precon->characteristic_reconstruction) {
      VectorDotRightEigenmatrix(pmb,IVX,il-1,iu,bx,wc,ql_iph);
      VectorDotRightEigenmatrix(pmb,IVX,il-1,iu,bx,wc,qr_imh);
    }

    // compute ql_(i+1/2) and qr_(i-1/2) using monotonized slopes
    for (int n=0; n<(NWAVE); ++n) {
#pragma simd
      for (int i=il-1; i<=iu; ++i){
        wl(n,k,j,i+1) = ql_iph(n,i);
        wr(n,k,j,i  ) = qr_imh(n,i);
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
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
#ifdef NO_COMPILE
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

    //  Apply van Leer limiter
    for (int n=0; n<(NWAVE); ++n) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        dw2(n,i) = dwl(n,i)*dwr(n,i);
        dwm(n,i) = 0.0;
        if (dw2(n,i) > 0.0) {
          dwm(n,i) = dw2(n,i)/(dwl(n,i) + dwr(n,i));
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
        wl(n,k,j+1,i) = wc(n,i) + dwm(n,i);
        wr(n,k,j  ,i) = wc(n,i) - dwm(n,i);
      }
    }
  }}

  dwl.DeleteAthenaArray();
  dwr.DeleteAthenaArray();
  dw2.DeleteAthenaArray();
  dwm.DeleteAthenaArray();
  wc.DeleteAthenaArray();
  bx.DeleteAthenaArray();
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::PiecewiseParabolicX3()
//  \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//         Colella-Sekora or Mignone limiting over [kl,ku][jl,ju][il,iu]

void Reconstruction::PiecewiseParabolicX3(MeshBlock *pmb,
  const int kl, const int ku, const int jl, const int ju, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
#ifdef NO_COMPLIE
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


    //  Apply van Leer limiter
    for (int n=0; n<(NWAVE); ++n) {
#pragma simd
      for (int i=il; i<=iu; ++i){
        dw2(n,i) = dwl(n,i)*dwr(n,i);
        dwm(n,i) = 0.0;
        if (dw2(n,i) > 0.0) {
          dwm(n,i) = dw2(n,i)/(dwl(n,i) + dwr(n,i));
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
        wl(n,k+1,j,i) = wc(n,i) + dwm(n,i);
        wr(n,k  ,j,i) = wc(n,i) - dwm(n,i);
      }
    }
  }}

  dwl.DeleteAthenaArray();
  dwr.DeleteAthenaArray();
  dw2.DeleteAthenaArray();
  dwm.DeleteAthenaArray();
  wc.DeleteAthenaArray();
  bx.DeleteAthenaArray();
#endif
  return;
}
