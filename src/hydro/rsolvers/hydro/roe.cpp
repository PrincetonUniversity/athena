//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file  roe.cpp
//! \brief Roe's linearized Riemann solver.
//!
//! Computes 1D fluxes using Roe's linearization.  When Roe's method fails because of
//! negative density in the intermediate states, LLF fluxes are used instead
//! (only density, not pressure, is checked in this version).
//!
//! REFERENCES:
//! - P. Roe, "Approximate Riemann solvers, parameter vectors, and difference schemes",
//!   JCP, 43, 357 (1981).

// C headers

// C++ headers
#include <algorithm>  // max()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"
#include "../../hydro.hpp"

namespace {
// prototype for function to compute Roe fluxes from eigenmatrices
inline void RoeFlux(const Real wroe[], const Real du[], const Real wli[], Real flx[],
                    Real eigenvalues[], int &flag);

// (gamma-1) and isothermal sound speed made global so can be shared with flux fn
Real gm1, iso_cs;
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//! \brief The Roe Riemann solver for hydrodynamics (both adiabatic and isothermal)

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr,
                          AthenaArray<Real> &flx, const AthenaArray<Real> &dxw) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NHYDRO)],wri[(NHYDRO)],wroe[(NHYDRO)];
  Real flxi[(NHYDRO)],fl[(NHYDRO)],fr[(NHYDRO)];
  gm1 = pmy_block->peos->GetGamma() - 1.0;
  iso_cs = pmy_block->peos->GetIsoSoundSpeed();

  Real ev[(NHYDRO)],du[(NHYDRO)];

#pragma omp simd private(wli,wri,wroe,flxi,fl,fr,ev,du)
  for (int i=il; i<=iu; ++i) {
    //--- Step 1.  Load L/R states into local variables
    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    if (NON_BAROTROPIC_EOS) wli[IPR]=wl(IPR,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    if (NON_BAROTROPIC_EOS) wri[IPR]=wr(IPR,i);

    //--- Step 2.  Compute Roe-averaged data from left- and right-states

    Real sqrtdl = std::sqrt(wli[IDN]);
    Real sqrtdr = std::sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN]  = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;

    // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
    // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    Real el,er;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IPR]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
      er = wri[IPR]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
      wroe[IPR] = ((el + wli[IPR])/sqrtdl + (er + wri[IPR])/sqrtdr)*isdlpdr;
    }

    //--- Step 3.  Compute L/R fluxes

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];

    fl[IDN] = mxl;
    fr[IDN] = mxr;

    fl[IVX] = mxl*wli[IVX];
    fr[IVX] = mxr*wri[IVX];

    fl[IVY] = mxl*wli[IVY];
    fr[IVY] = mxr*wri[IVY];

    fl[IVZ] = mxl*wli[IVZ];
    fr[IVZ] = mxr*wri[IVZ];

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IPR];
      fr[IVX] += wri[IPR];
      fl[IEN] = (el + wli[IPR])*wli[IVX];
      fr[IEN] = (er + wri[IPR])*wri[IVX];
    } else {
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

    //--- Step 4.  Compute Roe fluxes.

    du[IDN] = wri[IDN]          - wli[IDN];
    du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
    du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
    du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) du[IEN] = er - el;

    flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]);
    flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]);
    flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]);
    flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]);
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]);

    int llf_flag = 0;
    RoeFlux(wroe,du,wli,flxi,ev,llf_flag);

    //--- Step 5.  Overwrite with upwind flux if flow is supersonic

    if (ev[0] >= 0.0) {
      flxi[IDN] = fl[IDN];
      flxi[IVX] = fl[IVX];
      flxi[IVY] = fl[IVY];
      flxi[IVZ] = fl[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fl[IEN];
    }
    if (ev[NWAVE-1] <= 0.0) {
      flxi[IDN] = fr[IDN];
      flxi[IVX] = fr[IVX];
      flxi[IVY] = fr[IVY];
      flxi[IVZ] = fr[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fr[IEN];
    }

    //--- Step 6.  Overwrite with LLF flux if any of intermediate states are negative

    if (llf_flag != 0) {
      Real cl = pmy_block->peos->SoundSpeed(wli);
      Real cr = pmy_block->peos->SoundSpeed(wri);
      Real a  = 0.5*std::max( (std::abs(wli[IVX]) + cl), (std::abs(wri[IVX]) + cr) );

      flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]) - a*du[IDN];
      flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]) - a*du[IVX];
      flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]) - a*du[IVY];
      flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]) - a*du[IVZ];
      if (NON_BAROTROPIC_EOS) {
        flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]) - a*du[IEN];
      }
    }

    //--- Step 7. Store results into 3D array of fluxes

    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,k,j,i) = flxi[IEN];
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn RoeFlux()
//! \brief Computes Roe fluxes for the conserved variables.
//!
//! Computes Roe fluxes for the conserved variables, that is
//!           F[n] = 0.5*(F_l + F_r) - SUM_m(coeff[m]*rem[n][m])
//! where     coeff[n] = 0.5*ev[n]*SUM_m(dU[m]*lem[n][m])
//! and the rem[n][m] and lem[n][m] are matrices of the L- and R-eigenvectors of Roe's
//! matrix "A". Also returns the eigenvalues through the argument list.
//!
//! INPUT:
//!  - wroe: vector of Roe averaged primitive variables
//!  - du: Ur - Ul, difference in L/R-states in conserved variables
//!  - wli: Wl, left state in primitive variables
//!  - flx: (F_l + F_r)/2
//! OUTPUT:
//!  - flx: final Roe flux
//!  - ev: vector of eingenvalues
//!  - llf_flag: flag set to 1 if d<0 in any intermediate state
//! \note
//!  The order of the components in the input vectors should be:
//!     (IDN,IVX,IVY,IVZ,[IPR])
//!
//! REFERENCES:
//! - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//!   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.
#pragma omp declare simd simdlen(SIMD_WIDTH) notinbranch
inline void RoeFlux(const Real wroe[], const Real du[], const Real wli[], Real flx[],
                    Real ev[], int &llf_flag) {
  Real d  = wroe[IDN];
  Real v1 = wroe[IVX];
  Real v2 = wroe[IVY];
  Real v3 = wroe[IVZ];

  //--- Adiabatic hydrodynamics
  if (NON_BAROTROPIC_EOS) {
    Real h = wroe[IPR];
    Real vsq = v1*v1 + v2*v2 + v3*v3;
    Real q = h - 0.5*vsq;
    Real cs_sq = (q < 0.0) ? (TINY_NUMBER) : gm1*q;
    Real cs = std::sqrt(cs_sq);

    // Compute eigenvalues (eq. B2)
    ev[0] = v1 - cs;
    ev[1] = v1;
    ev[2] = v1;
    ev[3] = v1;
    ev[4] = v1 + cs;

    // Compute projection of dU onto L-eigenvectors using matrix elements from eq. B4
    Real a[(NHYDRO)];
    Real na = 0.5/cs_sq;
    a[0]  = du[0]*(0.5*gm1*vsq + v1*cs);
    a[0] -= du[1]*(gm1*v1 + cs);
    a[0] -= du[2]*gm1*v2;
    a[0] -= du[3]*gm1*v3;
    a[0] += du[4]*gm1;
    a[0] *= na;

    a[1]  = du[0]*(-v2);
    a[1] += du[2];

    a[2]  = du[0]*(-v3);
    a[2] += du[3];

    Real qa = gm1/cs_sq;
    a[3]  = du[0]*(1.0 - na*gm1*vsq);
    a[3] += du[1]*qa*v1;
    a[3] += du[2]*qa*v2;
    a[3] += du[3]*qa*v3;
    a[3] -= du[4]*qa;

    a[4]  = du[0]*(0.5*gm1*vsq - v1*cs);
    a[4] -= du[1]*(gm1*v1 - cs);
    a[4] -= du[2]*gm1*v2;
    a[4] -= du[3]*gm1*v3;
    a[4] += du[4]*gm1;
    a[4] *= na;

    Real coeff[(NHYDRO)];
    coeff[0] = -0.5*std::abs(ev[0])*a[0];
    coeff[1] = -0.5*std::abs(ev[1])*a[1];
    coeff[2] = -0.5*std::abs(ev[2])*a[2];
    coeff[3] = -0.5*std::abs(ev[3])*a[3];
    coeff[4] = -0.5*std::abs(ev[4])*a[4];

    // compute density in intermediate states and check that it is positive, set flag
    // This requires computing the [0][*] components of the right-eigenmatrix
    Real dens = wli[IDN] + a[0];  // rem[0][0]=1, so don't bother to compute or store
    if (dens < 0.0) llf_flag=1;

    dens += a[3];  // rem[0][3]=1, so don't bother to compute or store
    if (dens < 0.0) llf_flag=1;

    // Now multiply projection with R-eigenvectors from eq. B3 and SUM into output fluxes
    flx[0] += coeff[0];
    flx[0] += coeff[3];
    flx[0] += coeff[4];

    flx[1] += coeff[0]*(v1 - cs);
    flx[1] += coeff[3]*v1;
    flx[1] += coeff[4]*(v1 + cs);

    flx[2] += coeff[0]*v2;
    flx[2] += coeff[1];
    flx[2] += coeff[3]*v2;
    flx[2] += coeff[4]*v2;

    flx[3] += coeff[0]*v3;
    flx[3] += coeff[2];
    flx[3] += coeff[3]*v3;
    flx[3] += coeff[4]*v3;

    flx[4] += coeff[0]*(h - v1*cs);
    flx[4] += coeff[1]*v2;
    flx[4] += coeff[2]*v3;
    flx[4] += coeff[3]*0.5*vsq;
    flx[4] += coeff[4]*(h + v1*cs);

    //--- Isothermal hydrodynamics

  } else {
    // Compute eigenvalues (eq. B6)
    ev[0] = v1 - iso_cs;
    ev[1] = v1;
    ev[2] = v1;
    ev[3] = v1 + iso_cs;

    // Compute projection of dU onto L-eigenvectors using matrix elements from eq. B7
    Real a[(NHYDRO)];
    a[0]  = du[0]*(0.5 + 0.5*v1/iso_cs);
    a[0] -= du[1]*0.5/iso_cs;

    a[1]  = du[0]*(-v2);
    a[1] += du[2];

    a[2]  = du[0]*(-v3);
    a[2] += du[3];

    a[3]  = du[0]*(0.5 - 0.5*v1/iso_cs);
    a[3] += du[1]*0.5/iso_cs;

    Real coeff[(NHYDRO)];
    coeff[0] = -0.5*std::abs(ev[0])*a[0];
    coeff[1] = -0.5*std::abs(ev[1])*a[1];
    coeff[2] = -0.5*std::abs(ev[2])*a[2];
    coeff[3] = -0.5*std::abs(ev[3])*a[3];

    // compute density in intermediate states and check that it is positive, set flag
    // This requires computing the [0][*] components of the right-eigenmatrix
    Real dens = wli[IDN] + a[0];  // rem[0][0]=1, so don't bother to compute or store
    if (dens < 0.0) llf_flag=1;

    dens += a[3];  // rem[0][3]=1, so don't bother to compute or store
    if (dens < 0.0) llf_flag=1;

    // Now multiply projection with R-eigenvectors from eq. B3 and SUM into output fluxes
    flx[0] += coeff[0];
    flx[0] += coeff[3];

    flx[1] += coeff[0]*(v1 - iso_cs);
    flx[1] += coeff[3]*(v1 + iso_cs);

    flx[2] += coeff[0]*v2;
    flx[2] += coeff[1];
    flx[2] += coeff[3]*v2;

    flx[3] += coeff[0]*v3;
    flx[3] += coeff[2];
    flx[3] += coeff[3]*v3;
  }
  return;
}
} // namespace
