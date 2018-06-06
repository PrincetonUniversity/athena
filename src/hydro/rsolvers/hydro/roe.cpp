//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file  roe.cpp
//  \brief Roe's linearized Riemann solver.
//
// Computes 1D fluxes using Roe's linearization.  When Roe's method fails because of
// negative density or pressure in the intermediate states, LLF fluxes are used instead.
//
// REFERENCES:
// - P. Roe, "Approximate Riemann solvers, parameter vectors, and difference schemes",
//   JCP, 43, 357 (1981).

// C/C++ headers
#include <algorithm>  // max()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"

// prototype for functions to compute inner produces with eigenmatrices
inline void LeftRoeEigenmatrixDotVector(const Real wroe[], const Real in[], Real out[],
  Real eigenvalues[]);
inline void SumRightRoeEigenmatrixDotVector(const Real wroe[],const Real in[],Real out[],
  int &flag);

// (gamma-1) and isothermal sound speed made global so can be shared with eigensystem fns
static Real gm1, iso_cs;

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//  \brief The Roe Riemann solver for hydrodynamics (both adiabatic and isothermal)

void Hydro::RiemannSolver(const int kl, const int ku, const int jl, const int ju,
  const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
  AthenaArray<Real> &ey, AthenaArray<Real> &ez) {

  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NHYDRO)],wri[(NHYDRO)],wroe[(NHYDRO)];
  Real flxi[(NHYDRO)],fl[(NHYDRO)],fr[(NHYDRO)];
  gm1 = pmy_block->peos->GetGamma() - 1.0;
  iso_cs = pmy_block->peos->GetIsoSoundSpeed();

  Real coeff[(NHYDRO)];
  Real ev[(NHYDRO)],du[(NHYDRO)],a[(NHYDRO)],u[(NHYDRO)];

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {

//--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,k,j,i);
    wli[IVX]=wl(ivx,k,j,i);
    wli[IVY]=wl(ivy,k,j,i);
    wli[IVZ]=wl(ivz,k,j,i);
    if (NON_BAROTROPIC_EOS) wli[IPR]=wl(IPR,k,j,i);

    wri[IDN]=wr(IDN,k,j,i);
    wri[IVX]=wr(ivx,k,j,i);
    wri[IVY]=wr(ivy,k,j,i);
    wri[IVZ]=wr(ivz,k,j,i);
    if (NON_BAROTROPIC_EOS) wri[IPR]=wr(IPR,k,j,i);

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

//--- Step 4.  Compute projection of dU onto L eigenvectors ("vector A"), and eigenvectors

    du[IDN] = wri[IDN]          - wli[IDN];
    du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
    du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
    du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) du[IEN] = er - el;

    LeftRoeEigenmatrixDotVector(wroe,du,a,ev);

//--- Step 5.  Check that the density and pressure in the intermediate states are
// positive.  If not, set a flag that will be checked below.

    int llf_flag = 0;
    u[IDN] = wli[IDN];
    u[IVX] = wli[IDN]*wli[IVX];
    u[IVY] = wli[IDN]*wli[IVY];
    u[IVZ] = wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) u[IEN] = el;

    SumRightRoeEigenmatrixDotVector(wroe,a,u,llf_flag);

//--- Step 6.  Compute Roe flux

    coeff[IDN] = -0.5*fabs(ev[IDN])*a[IDN];
    coeff[IVX] = -0.5*fabs(ev[IVX])*a[IVX];
    coeff[IVY] = -0.5*fabs(ev[IVY])*a[IVY];
    coeff[IVZ] = -0.5*fabs(ev[IVZ])*a[IVZ];
    if (NON_BAROTROPIC_EOS) coeff[IEN] = -0.5*fabs(ev[IEN])*a[IEN];

    flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]);
    flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]);
    flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]);
    flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]);
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]);

    int dummy_flag=0;
    SumRightRoeEigenmatrixDotVector(wroe,coeff,flxi,dummy_flag);

//--- Step 7.  Overwrite with upwind flux if flow is supersonic

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

//--- Step 8.  Overwrite with LLF flux if any of intermediate states are negative

    if (llf_flag != 0) {
      Real a = 0.5*std::max(fabs(ev[0]), fabs(ev[NWAVE-1]));

      flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]) - a*du[IDN];
      flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]) - a*du[IVX];
      flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]) - a*du[IVY];
      flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]) - a*du[IVZ];
      if (NON_BAROTROPIC_EOS) {
        flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]) - a*du[IEN];
      }
    }

    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,k,j,i) = flxi[IEN];
  }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn LeftRoeEigenmatrixDotVector()
//  \brief Computes inner-product of left-eigenmatrix of Roe's matrix A in the primitive
//  variables and an input vector, and returns in out[].  The input vector is UNCHANGED.
//
//  The order of the components in the input vector should be:
//     (IDN,IVX,IVY,IVZ,[IPR],[IBY,IBZ])
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.

inline void LeftRoeEigenmatrixDotVector(const Real wroe[], const Real in[], Real out[],
                                        Real eigenvalues[]) {
  Real d  = wroe[IDN];
  Real v1 = wroe[IVX];
  Real v2 = wroe[IVY];
  Real v3 = wroe[IVZ];

//--- Adiabatic hydrodynamics

  if (NON_BAROTROPIC_EOS) {
    Real h = wroe[IPR];
    Real vsq = v1*v1 + v2*v2 + v3*v3;
    Real q = h - 0.5*vsq;
    Real asq = (q < 0.0) ? 0.0 : gm1*q;
    Real a = std::sqrt(asq);

    // Compute eigenvalues (eq. B2)
    eigenvalues[0] = v1 - a;
    eigenvalues[1] = v1;
    eigenvalues[2] = v1;
    eigenvalues[3] = v1;
    eigenvalues[4] = v1 + a;

    // Multiply row of L-eigenmatrix with vector using matrix elements from eq. B4
    Real na = 0.5/asq;
    out[IDN] = (0.5*gm1*vsq + v1*a)*in[IDN];
    out[IDN] -= gm1*(v1*in[IVX] + v2*in[IVY] + v3*in[IVZ] - in[IPR]) + a*in[IVX];
    out[IDN] *= na;

    out[IVX] = -v2*in[IDN] + in[IVY];
    out[IVY] = -v3*in[IDN] + in[IVZ];

    Real qa = gm1/asq;
    out[IVZ] = (1.0 - na*gm1*vsq)*in[IDN];
    out[IVZ] += qa*(v1*in[IVX] + v2*in[IVY] + v3*in[IVZ] - in[IPR]);

    out[IPR] = (0.5*gm1*vsq - v1*a)*in[IDN];
    out[IPR] -= gm1*(v1*in[IVX] + v2*in[IVY] + v3*in[IVZ] - in[IPR]) - a*in[IVX];
    out[IPR] *= na;

//--- Isothermal hydrodynamics

  } else {
    // Compute eigenvalues (eq. B6)
    eigenvalues[0] = v1 - iso_cs;
    eigenvalues[1] = v1;
    eigenvalues[2] = v1;
    eigenvalues[3] = v1 + iso_cs;

    // Multiply row of L-eigenmatrix with vector using matrix elements from eq. B7
    out[IDN] = 0.5*(in[IDN] + (v1*in[IDN] - in[IVX])/iso_cs);
    out[IVX] = -v2*in[IDN] + in[IVY];
    out[IVY] = -v3*in[IDN] + in[IVZ];
    out[IVZ] = 0.5*(in[IDN] - (v1*in[IDN] - in[IVX])/iso_cs);
  }
}

//----------------------------------------------------------------------------------------
//! \fn SumRightRoeEigenmatrixDotVector()
//  \brief Computes inner-product of right-eigenmatrix of Roe's matrix A in the primitive
//  variables and an input vector, and sums into output vector.  Also checks that the
//  density and pressure remain positive in the intermediate states, else an output flag
//  is set.  This only works if the input vector is L-state of conserved variables.
//  The input vector is UNCHANGED.
//
//  The order of the components in the input vector should be:
//     (IDN,IVX,IVY,IVZ,[IPR],[IBY,IBZ])
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.

inline void SumRightRoeEigenmatrixDotVector(const Real wroe[],const Real in[],Real out[],
                                            int &flag) {
  Real d  = wroe[IDN];
  Real v1 = wroe[IVX];
  Real v2 = wroe[IVY];
  Real v3 = wroe[IVZ];

//--- Adiabatic hydrodynamics

  if (NON_BAROTROPIC_EOS) {
    Real h = wroe[IPR];
    Real vsq = v1*v1 + v2*v2 + v3*v3;
    Real q = h - 0.5*vsq;
    Real asq = (q < 0.0) ? 0.0 : gm1*q;
    Real a = std::sqrt(asq);

    // Multiply row of R-eigenmatrix from eq. B3 with vector and sum into output
    // Also check density and pressure are positive, and set flag if not
    // jump across wave[0]
    out[IDN] += in[IDN];
    out[IVX] += in[IDN]*(v1 - a);
    out[IVY] += in[IDN]*v2;
    out[IVZ] += in[IDN]*v3;
    out[IEN] += in[IDN]*(h - v1*a);
    Real p = out[IEN] - 0.5*(SQR(out[IVX])+SQR(out[IVY])+SQR(out[IVZ]))/out[IDN];
    if (out[IDN] < 0.0) flag=1;
    if (p < 0.0) flag=2;

    // jump across wave[1]
    out[IVY] += in[IVX];
    out[IEN] += in[IVX]*v2;
    p = out[IEN] - 0.5*(SQR(out[IVX])+SQR(out[IVY])+SQR(out[IVZ]))/out[IDN];
    if (p < 0.0) flag=2;

    // jump across wave[2]
    out[IVZ] += in[IVY];
    out[IEN] += in[IVY]*v3;
    p = out[IEN] - 0.5*(SQR(out[IVX])+SQR(out[IVY])+SQR(out[IVZ]))/out[IDN];
    if (p < 0.0) flag=2;

    // jump across wave[3]
    out[IDN] += in[IVZ];
    out[IVX] += in[IVZ]*v1;
    out[IVY] += in[IVZ]*v2;
    out[IVZ] += in[IVZ]*v3;
    out[IEN] += in[IVZ]*0.5*asq;
    p = out[IEN] - 0.5*(SQR(out[IVX])+SQR(out[IVY])+SQR(out[IVZ]))/out[IDN];
    if (out[IDN] < 0.0) flag=1;
    if (p < 0.0) flag=2;

    // jump across wave[4]
    out[IDN] += in[IEN];
    out[IVX] += in[IEN]*(v1+a);
    out[IVY] += in[IEN]*v2;
    out[IVZ] += in[IEN]*v3;
    out[IEN] += in[IEN]*(h + v1*a);
    p = out[IEN] - 0.5*(SQR(out[IVX])+SQR(out[IVY])+SQR(out[IVZ]))/out[IDN];
    if (out[IDN] < 0.0) flag=1;
    if (p < 0.0) flag=2;

//--- Isothermal hydrodynamics

  } else {

    // Multiply row of R-eigenmatrix from eq. B3 with vector and sum into output
    // Also check density and pressure are positive, and set flag if not
    // jump across wave[0]
    out[IDN] += in[IDN];
    out[IVX] += in[IDN]*(v1 - iso_cs);
    out[IVY] += in[IDN]*v2;
    out[IVZ] += in[IDN]*v3;
    if (out[IDN] < 0.0) flag=1;

    // jump across wave[1]
    out[IVY] += in[IVX];

    // jump across wave[2]
    out[IVZ] += in[IVY];

    // jump across wave[3]
    out[IDN] += in[IEN];
    out[IVX] += in[IEN]*(v1 + iso_cs);
    out[IVY] += in[IEN]*v2;
    out[IVZ] += in[IEN]*v3;
    if (out[IDN] < 0.0) flag=1;
  }
}
