//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file  roe_mhd.cpp
//  \brief Roe's linearized Riemann solver for MHD.
//
// Computes 1D fluxes using Roe's linearization.  When Roe's method fails because of
// negative density in the intermediate states, LLF fluxes are used instead (only density,
// not pressure, is checked in this version).
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

// prototype for functions to compute inner product with eigenmatrices
inline void RoeFlux(const Real wroe[], const Real b1, const Real x, const Real y,
  const Real du[], const Real wli[], Real flx[], Real eigenvalues[], int &flag);

// (gamma-1) and isothermal sound speed made global so can be shared with eigensystem fns
static Real gm1, iso_cs;

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//  \brief The Roe Riemann solver for MHD (both adiabatic and isothermal)

void Hydro::RiemannSolver(const int kl, const int ku, const int jl, const int ju,
  const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
  AthenaArray<Real> &ey, AthenaArray<Real> &ez) {

  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NWAVE)],wri[(NWAVE)],wroe[(NWAVE)];
  Real flxi[(NWAVE)],fl[(NWAVE)],fr[(NWAVE)];
  gm1 = pmy_block->peos->GetGamma() - 1.0;
  iso_cs = pmy_block->peos->GetIsoSoundSpeed();

  Real ev[(NWAVE)],du[(NWAVE)];

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
#pragma omp simd private(wli,wri,wroe,flxi,fl,fr,ev,du)
  for (int i=il; i<=iu; ++i) {

//--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,k,j,i);
    wli[IVX]=wl(ivx,k,j,i);
    wli[IVY]=wl(ivy,k,j,i);
    wli[IVZ]=wl(ivz,k,j,i);
    if (NON_BAROTROPIC_EOS) wli[IPR]=wl(IPR,k,j,i);
    wli[IBY]=wl(IBY,k,j,i);
    wli[IBZ]=wl(IBZ,k,j,i);

    wri[IDN]=wr(IDN,k,j,i);
    wri[IVX]=wr(ivx,k,j,i);
    wri[IVY]=wr(ivy,k,j,i);
    wri[IVZ]=wr(ivz,k,j,i);
    if (NON_BAROTROPIC_EOS) wri[IPR]=wr(IPR,k,j,i);
    wri[IBY]=wr(IBY,k,j,i);
    wri[IBZ]=wr(IBZ,k,j,i);

    Real bxi = bx(k,j,i);

//--- Step 2.  Compute Roe-averaged data from left- and right-states

    Real sqrtdl = std::sqrt(wli[IDN]);
    Real sqrtdr = std::sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN] = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;
    // Note Roe average of magnetic field is different
    wroe[IBY] = (sqrtdr*wli[IBY] + sqrtdl*wri[IBY])*isdlpdr;
    wroe[IBZ] = (sqrtdr*wli[IBZ] + sqrtdl*wri[IBZ])*isdlpdr;
    Real x = 0.5*(SQR(wli[IBY]-wri[IBY]) + SQR(wli[IBZ]-wri[IBZ]))/(SQR(sqrtdl+sqrtdr));
    Real y = 0.5*(wli[IDN] + wri[IDN])/wroe[IDN];

    // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
    // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    Real pbl = 0.5*(bxi*bxi + SQR(wli[IBY]) + SQR(wli[IBZ]));
    Real pbr = 0.5*(bxi*bxi + SQR(wri[IBY]) + SQR(wri[IBZ]));
    Real el,er;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IPR]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX])+SQR(wli[IVY])+SQR(wli[IVZ])) +pbl;
      er = wri[IPR]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX])+SQR(wri[IVY])+SQR(wri[IVZ])) +pbr;
      wroe[IPR] = ((el + wli[IPR] + pbl)/sqrtdl + (er + wri[IPR] + pbr)/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute L/R fluxes

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];

    fl[IDN] = mxl;
    fr[IDN] = mxr;

    fl[IVX] = mxl*wli[IVX] + pbl - SQR(bxi);
    fr[IVX] = mxr*wri[IVX] + pbr - SQR(bxi);

    fl[IVY] = mxl*wli[IVY] - bxi*wli[IBY];
    fr[IVY] = mxr*wri[IVY] - bxi*wri[IBY];

    fl[IVZ] = mxl*wli[IVZ] - bxi*wli[IBZ];
    fr[IVZ] = mxr*wri[IVZ] - bxi*wri[IBZ];

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IPR];
      fr[IVX] += wri[IPR];
      fl[IEN] = (el + wli[IPR] + pbl - bxi*bxi)*wli[IVX];
      fr[IEN] = (er + wri[IPR] + pbr - bxi*bxi)*wri[IVX];
      fl[IEN] -= bxi*(wli[IBY]*wli[IVY] + wli[IBZ]*wli[IVZ]);
      fr[IEN] -= bxi*(wri[IBY]*wri[IVY] + wri[IBZ]*wri[IVZ]);
    } else {
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

    fl[IBY] = wli[IBY]*wli[IVX] - bxi*wli[IVY];
    fr[IBY] = wri[IBY]*wri[IVX] - bxi*wri[IVY];

    fl[IBZ] = wli[IBZ]*wli[IVX] - bxi*wli[IVZ];
    fr[IBZ] = wri[IBZ]*wri[IVX] - bxi*wri[IVZ];

//--- Step 4.  Compute Roe fluxes

    du[IDN] = wri[IDN]          - wli[IDN];
    du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
    du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
    du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) du[IEN] = er - el;
    du[IBY] = wri[IBY] - wli[IBY];
    du[IBZ] = wri[IBZ] - wli[IBZ];

    flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]);
    flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]);
    flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]);
    flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]);
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]);
    flxi[IBY] = 0.5*(fl[IBY] + fr[IBY]);
    flxi[IBZ] = 0.5*(fl[IBZ] + fr[IBZ]);

    int llf_flag = 0;
    RoeFlux(wroe,bxi,x,y,du,wli,flxi,ev,llf_flag);

//--- Step 5.  Overwrite with upwind flux if flow is supersonic

    if (ev[0] >= 0.0) {
      flxi[IDN] = fl[IDN];
      flxi[IVX] = fl[IVX];
      flxi[IVY] = fl[IVY];
      flxi[IVZ] = fl[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fl[IEN];
      flxi[IBY] = fl[IBY];
      flxi[IBZ] = fl[IBZ];
    }
    if (ev[NWAVE-1] <= 0.0) {
      flxi[IDN] = fr[IDN];
      flxi[IVX] = fr[IVX];
      flxi[IVY] = fr[IVY];
      flxi[IVZ] = fr[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fr[IEN];
      flxi[IBY] = fr[IBY];
      flxi[IBZ] = fr[IBZ];
    }

//--- Step 6.  Overwrite with LLF flux if any of intermediate states are negative

    if (llf_flag != 0) {
      Real cfl = pmy_block->peos->FastMagnetosonicSpeed(wli,bxi);
      Real cfr = pmy_block->peos->FastMagnetosonicSpeed(wri,bxi);
      Real a = 0.5*std::max( (fabs(wli[IVX]) + cfl), (fabs(wri[IVX]) + cfr) );

      flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]) - a*du[IDN];
      flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]) - a*du[IVX];
      flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]) - a*du[IVY];
      flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]) - a*du[IVZ];
      if (NON_BAROTROPIC_EOS) {
        flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]) - a*du[IEN];
      }
      flxi[IBY] = 0.5*(fl[IBY] + fr[IBY]) - a*du[IBY];
      flxi[IBZ] = 0.5*(fl[IBZ] + fr[IBZ]) - a*du[IBZ];
    }

//--- Step 7. Store results into 3D array of fluxes

    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,k,j,i) = flxi[IEN];
    ey(k,j,i) = -flxi[IBY];
    ez(k,j,i) =  flxi[IBZ];
  }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn RoeFlux()
//  \brief Computes Roe fluxes for the conserved variables, that is
//            F[n] = 0.5*(F_l + F_r) - SUM_m(coeff[m]*rem[n][m])
//  where     coeff[n] = 0.5*ev[n]*SUM_m(dU[m]*lem[n][m])
//  and the rem[n][m] and lem[n][m] are matrices of the L- and R-eigenvectors of Roe's
//  matrix "A". Also returns the eigenvalues through the argument list.
//
// INPUT:
//   wroe: vector of Roe averaged primitive variables
//   du: Ur - Ul, difference in L/R-states in conserved variables
//   wli: Wl, left state in primitive variables
//   flx: (F_l + F_r)/2
//
// OUTPUT:
//   flx: final Roe flux
//   ev: vector of eingenvalues
//   llf_flag: flag set to 1 if d<0 in any intermediate state
//
//  The order of the components in the input vectors should be:
//     (IDN,IVX,IVY,IVZ,[IPR])
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper.
#pragma omp declare simd simdlen(SIMD_WIDTH) notinbranch
inline void RoeFlux(const Real wroe[], const Real b1, const Real x, const Real y,
  const Real du[], const Real wli[], Real flx[], Real ev[], int &llf_flag) {

  Real d  = wroe[IDN];
  Real v1 = wroe[IVX];
  Real v2 = wroe[IVY];
  Real v3 = wroe[IVZ];
  Real b2 = wroe[IBY];
  Real b3 = wroe[IBZ];

  // compute sound and Alfven speeds
  Real di = 1.0/d;
  Real btsq = b2*b2 + b3*b3;
  Real vaxsq = b1*b1*di;

  Real bt_starsq, twid_csq, vsq, hp;
  if (NON_BAROTROPIC_EOS) {
    vsq = v1*v1 + v2*v2 + v3*v3;
    hp = wroe[IPR] - (vaxsq + btsq*di);
    bt_starsq = (gm1 - (gm1 - 1.0)*y)*btsq;
    twid_csq = std::max((gm1*(hp-0.5*vsq)-(gm1-1.0)*x), TINY_NUMBER);
  } else {
    bt_starsq = btsq*y;
    twid_csq = (iso_cs*iso_cs) + x;
  }

  // Compute fast- and slow-magnetosonic speeds (eq. B18)
  Real ct2 = bt_starsq*di;
  Real tsum = vaxsq + ct2 + twid_csq;
  Real tdif = vaxsq + ct2 - twid_csq;
  Real cf2_cs2 = std::sqrt(tdif*tdif + 4.0*twid_csq*ct2);

  Real cfsq = 0.5*(tsum + cf2_cs2);
  Real cf = std::sqrt(cfsq);

  Real cssq = twid_csq*vaxsq/cfsq;
  Real cs = std::sqrt(cssq);

  // Compute beta(s) (eqs. A17, B20, B28)
  Real bt = std::sqrt(btsq);
  Real bt_star = std::sqrt(bt_starsq);
  Real bet2 = 0.0;
  Real bet3 = 0.0;
  if (bt != 0.0) {
    bet2 = b2/bt;
    bet3 = b3/bt;
  }

  Real bet2_star,bet3_star;
  if (NON_BAROTROPIC_EOS) {
    bet2_star = bet2/std::sqrt(gm1 - (gm1-1.0)*y);
    bet3_star = bet3/std::sqrt(gm1 - (gm1-1.0)*y);
  } else {
    bet2_star = bet2/std::sqrt(y);
    bet3_star = bet3/std::sqrt(y);
  }
  Real bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
  Real vbet = v2*bet2_star + v3*bet3_star;
  Real q2_star = 0.0;
  Real q3_star = 0.0;
  if (bet_starsq != 0.0) {
    q2_star = bet2_star/bet_starsq;
    q3_star = bet3_star/bet_starsq;
  }

  // Compute alpha(s) (eq. A16)
  Real alpha_f, alpha_s;
  if ((cfsq - cssq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ((twid_csq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ((cfsq - twid_csq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = std::sqrt((twid_csq - cssq)/(cfsq - cssq));
    alpha_s = std::sqrt((cfsq - twid_csq)/(cfsq - cssq));
  }

  // Compute Q(s) and A(s) (eq. A14-15), etc.
  Real sqrtd = std::sqrt(d);
  Real isqrtd = 1.0/sqrtd;
  Real s = SIGN(b1);
  Real twid_c = std::sqrt(twid_csq);
  Real qf = cf*alpha_f*s;
  Real qs = cs*alpha_s*s;
  Real af_prime = twid_c*alpha_f*isqrtd;
  Real as_prime = twid_c*alpha_s*isqrtd;
  Real afpbb = af_prime*bt_star*bet_starsq;
  Real aspbb = as_prime*bt_star*bet_starsq;
  Real vqstr = (v2*q2_star + v3*q3_star);
  Real vax = std::sqrt(vaxsq);

  // Normalize by 1/2a^{2}: quantities denoted by \hat{f}
  Real norm = 0.5/twid_csq;
  Real cff = norm*alpha_f*cf;
  Real css = norm*alpha_s*cs;
  Real qf_hat = qf*norm;
  Real qs_hat = qs*norm;
  Real af = norm*af_prime*d;
  Real as = norm*as_prime*d;
  Real afpb = norm*af_prime*bt_star;
  Real aspb = norm*as_prime*bt_star;

//--- Adiabatic MHD

  if (NON_BAROTROPIC_EOS) {
    // Compute eigenvalues (eq. B17)
    ev[0] = v1 - cf;
    ev[1] = v1 - vax;
    ev[2] = v1 - cs;
    ev[3] = v1;
    ev[4] = v1 + cs;
    ev[5] = v1 + vax;
    ev[6] = v1 + cf;

    // Compute projection of dU onto L-eigenvectors using matrix elements from eq. B29
    // Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f}
    Real a[(NWAVE)];
    Real alpha_f_bar = alpha_f*gm1*norm;
    Real alpha_s_bar = alpha_s*gm1*norm;
    Real gm1a = gm1/twid_csq;

    a[0]  = du[0]*(alpha_f_bar*(vsq-hp) + cff*(cf+v1) - qs_hat*vqstr - aspb);
    a[0] -= du[1]*(alpha_f_bar*v1 + cff);
    a[0] -= du[2]*(alpha_f_bar*v2 - qs_hat*q2_star);
    a[0] -= du[3]*(alpha_f_bar*v3 - qs_hat*q3_star);
    a[0] += du[4]*alpha_f_bar;
    a[0] += du[5]*(as*q2_star - alpha_f_bar*b2);
    a[0] += du[6]*(as*q3_star - alpha_f_bar*b3);

    a[1]  = du[0]*(v2*bet3 - v3*bet2);
    a[1] -= du[2]*bet3;
    a[1] += du[3]*bet2;
    a[1] -= du[5]*sqrtd*bet3*s;
    a[1] += du[6]*sqrtd*bet2*s;
    a[1] *= 0.5;

    a[2]  = du[0]*(alpha_s_bar*(vsq-hp) + css*(cs+v1) + qf_hat*vqstr + afpb);
    a[2] -= du[1]*(alpha_s_bar*v1 + css);
    a[2] -= du[2]*(alpha_s_bar*v2 + qf_hat*q2_star);
    a[2] -= du[3]*(alpha_s_bar*v3 + qf_hat*q3_star);
    a[2] += du[4]*alpha_s_bar;
    a[2] -= du[5]*(af*q2_star + alpha_s_bar*b2);
    a[2] -= du[6]*(af*q3_star + alpha_s_bar*b3);

    a[3]  = du[0]*(1.0 - gm1a*(0.5*vsq - (gm1-1.0)*x/gm1));
    a[3] += du[1]*gm1a*v1;
    a[3] += du[2]*gm1a*v2;
    a[3] += du[3]*gm1a*v3;
    a[3] -= du[4]*gm1a;
    a[3] += du[5]*gm1a*b2;
    a[3] += du[6]*gm1a*b3;

    a[4]  = du[0]*(alpha_s_bar*(vsq-hp) + css*(cs-v1) - qf_hat*vqstr + afpb);
    a[4] -= du[1]*(alpha_s_bar*v1 - css);
    a[4] -= du[2]*(alpha_s_bar*v2 - qf_hat*q2_star);
    a[4] -= du[3]*(alpha_s_bar*v3 - qf_hat*q3_star);
    a[4] += du[4]*alpha_s_bar;
    a[4] -= du[5]*(af*q2_star + alpha_s_bar*b2);
    a[4] -= du[6]*(af*q3_star + alpha_s_bar*b3);

    a[5]  = du[0]*(v3*bet2 - v2*bet3);
    a[5] += du[2]*bet3;
    a[5] -= du[3]*bet2;
    a[5] -= du[5]*sqrtd*bet3*s;
    a[5] += du[6]*sqrtd*bet2*s;
    a[5] *= 0.5;

    a[6]  = du[0]*(alpha_f_bar*(vsq-hp) + cff*(cf-v1) + qs_hat*vqstr - aspb);
    a[6] -= du[1]*(alpha_f_bar*v1 - cff);
    a[6] -= du[2]*(alpha_f_bar*v2 + qs_hat*q2_star);
    a[6] -= du[3]*(alpha_f_bar*v3 + qs_hat*q3_star);
    a[6] += du[4]*alpha_f_bar;
    a[6] += du[5]*(as*q2_star - alpha_f_bar*b2);
    a[6] += du[6]*(as*q3_star - alpha_f_bar*b3);

    Real coeff[(NWAVE)];
    coeff[0] = -0.5*fabs(ev[0])*a[0];
    coeff[1] = -0.5*fabs(ev[1])*a[1];
    coeff[2] = -0.5*fabs(ev[2])*a[2];
    coeff[3] = -0.5*fabs(ev[3])*a[3];
    coeff[4] = -0.5*fabs(ev[4])*a[4];
    coeff[5] = -0.5*fabs(ev[5])*a[5];
    coeff[6] = -0.5*fabs(ev[6])*a[6];

    // compute density in intermediate states and check that it is positive, set flag
    // This uses the [0][*] components of the right-eigenmatrix
    Real dens = wli[IDN] + a[0]*alpha_f;
    if (dens < 0.0) llf_flag=1;

    dens += a[2]*alpha_s;
    if (dens < 0.0) llf_flag=1;

    dens += a[3];
    if (dens < 0.0) llf_flag=1;

    dens += a[4]*alpha_s;
    if (dens < 0.0) llf_flag=1;

    // Now multiply projection with R-eigenvectors from eq. B21 and SUM into output fluxes
    flx[0] += coeff[0]*alpha_f;
    flx[0] += coeff[2]*alpha_s;
    flx[0] += coeff[3];
    flx[0] += coeff[4]*alpha_s;
    flx[0] += coeff[6]*alpha_f;

    flx[1] += coeff[0]*(alpha_f*(v1 - cf));
    flx[1] += coeff[2]*(alpha_s*(v1 - cs));
    flx[1] += coeff[3]*v1;
    flx[1] += coeff[4]*(alpha_s*(v1 + cs));
    flx[1] += coeff[6]*(alpha_f*(v1 + cf));

    flx[2] += coeff[0]*(alpha_f*v2 + qs*bet2_star);
    flx[2] -= coeff[1]*bet3;
    flx[2] += coeff[2]*(alpha_s*v2 - qf*bet2_star);
    flx[2] += coeff[3]*v2;
    flx[2] += coeff[4]*(alpha_s*v2 + qf*bet2_star);
    flx[2] += coeff[5]*bet3;
    flx[2] += coeff[6]*(alpha_f*v2 - qs*bet2_star);

    flx[3] += coeff[0]*(alpha_f*v3 + qs*bet3_star);
    flx[3] += coeff[1]*bet2;
    flx[3] += coeff[2]*(alpha_s*v3 - qf*bet3_star);
    flx[3] += coeff[3]*v3;
    flx[3] += coeff[4]*(alpha_s*v3 + qf*bet3_star);
    flx[3] -= coeff[5]*bet2;
    flx[3] += coeff[6]*(alpha_f*v3 - qs*bet3_star);

    flx[4] += coeff[0]*(alpha_f*(hp - v1*cf) + qs*vbet + aspbb);
    flx[4] -= coeff[1]*(v2*bet3 - v3*bet2);
    flx[4] += coeff[2]*(alpha_s*(hp - v1*cs) - qf*vbet - afpbb);
    flx[4] += coeff[3]*(0.5*vsq + (gm1-1.0)*x/gm1);
    flx[4] += coeff[4]*(alpha_s*(hp + v1*cs) + qf*vbet - afpbb);
    flx[4] += coeff[5]*(v1*bet3 - v3*bet2);
    flx[4] += coeff[6]*(alpha_f*(hp + v1*cf) - qs*vbet + aspbb);

    flx[5] += coeff[0]*as_prime*bet2_star;
    flx[5] -= coeff[1]*bet3*s*isqrtd;
    flx[5] -= coeff[2]*af_prime*bet2_star;
    flx[5] -= coeff[4]*af_prime*bet2_star;
    flx[5] -= coeff[5]*bet3*s*isqrtd;
    flx[5] += coeff[6]*as_prime*bet2_star;

    flx[6] += coeff[0]*as_prime*bet3_star;
    flx[6] += coeff[1]*bet2*s*isqrtd;
    flx[6] -= coeff[2]*af_prime*bet3_star;
    flx[6] -= coeff[4]*af_prime*bet3_star;
    flx[6] += coeff[5]*bet2*s*isqrtd;
    flx[6] += coeff[6]*as_prime*bet3_star;

//--- Isothermal MHD

  } else {
    // Compute eigenvalues (eq. B38)
    ev[0] = v1 - cf;
    ev[1] = v1 - vax;
    ev[2] = v1 - cs;
    ev[3] = v1 + cs;
    ev[4] = v1 + vax;
    ev[5] = v1 + cf;

    // Compute projection of dU onto L-eigenvectors using matrix elements from eq. B41
    Real a[(NWAVE)];
    a[0]  = du[0]*(cff*(cf+v1) - qs_hat*vqstr - aspb);
    a[0] -= du[1]*cff;
    a[0] += du[2]*qs_hat*q2_star;
    a[0] += du[3]*qs_hat*q3_star;
    a[0] += du[4]*as*q2_star;
    a[0] += du[5]*as*q3_star;

    a[1]  = du[0]*(v2*bet3 - v3*bet2);
    a[1] -= du[2]*bet3;
    a[1] += du[3]*bet2;
    a[1] -= du[4]*sqrtd*bet3*s;
    a[1] += du[5]*sqrtd*bet2*s;
    a[1] *= 0.5;

    a[2]  = du[0]*(css*(cs+v1) + qf_hat*vqstr + afpb);
    a[2] -= du[1]*css;
    a[2] -= du[2]*qf_hat*q2_star;
    a[2] -= du[3]*qf_hat*q3_star;
    a[2] -= du[4]*af*q2_star;
    a[2] -= du[5]*af*q3_star;

    a[3]  = du[0]*(css*(cs-v1) - qf_hat*vqstr + afpb);
    a[3] += du[1]*css;
    a[3] += du[2]*qf_hat*q2_star;
    a[3] += du[3]*qf_hat*q3_star;
    a[3] -= du[4]*af*q2_star;
    a[3] -= du[5]*af*q3_star;

    a[4]  = du[0]*(v3*bet2 - v2*bet3);
    a[4] += du[2]*bet3;
    a[4] -= du[3]*bet2;
    a[4] -= du[4]*sqrtd*bet3*s;
    a[4] += du[5]*sqrtd*bet2*s;
    a[4] *= 0.5;

    a[5]  = du[0]*(cff*(cf-v1) + qs_hat*vqstr - aspb);
    a[5] += du[1]*cff;
    a[5] -= du[2]*qs_hat*q2_star;
    a[5] -= du[3]*qs_hat*q3_star;
    a[5] += du[4]*as*q2_star;
    a[5] += du[5]*as*q3_star;

    Real coeff[(NWAVE)];
    coeff[IDN] = -0.5*fabs(ev[IDN])*a[IDN];
    coeff[IVX] = -0.5*fabs(ev[IVX])*a[IVX];
    coeff[IVY] = -0.5*fabs(ev[IVY])*a[IVY];
    coeff[IVZ] = -0.5*fabs(ev[IVZ])*a[IVZ];
    if (NON_BAROTROPIC_EOS) coeff[IEN] = 0.5*fabs(ev[IEN])*a[IEN];
    coeff[IBY] = -0.5*fabs(ev[IBY])*a[IBY];
    coeff[IBZ] = -0.5*fabs(ev[IBZ])*a[IBZ];

    // compute density in intermediate states and check that it is positive, set flag
    // This uses the [0][*] components of the right-eigenmatrix
    Real dens = wli[IDN] + a[0]*alpha_f;
    if (dens < 0.0) llf_flag=1;

    dens += a[2]*alpha_s;
    if (dens < 0.0) llf_flag=1;

    dens += a[3]*alpha_s;
    if (dens < 0.0) llf_flag=1;

    // Now multiply projection with R-eigenvectors from eq. B21 and SUM into output fluxes
    flx[0] += coeff[0]*alpha_f;
    flx[0] += coeff[2]*alpha_s;
    flx[0] += coeff[3]*alpha_s;
    flx[0] += coeff[5]*alpha_f;

    flx[1] += coeff[0]*alpha_f*(v1 - cf);
    flx[1] += coeff[2]*alpha_s*(v1 - cs);
    flx[1] += coeff[3]*alpha_s*(v1 + cs);
    flx[1] += coeff[5]*alpha_f*(v1 + cf);

    flx[2] += coeff[0]*(alpha_f*v2 + qs*bet2_star);
    flx[2] -= coeff[1]*bet3;
    flx[2] += coeff[2]*(alpha_s*v2 - qf*bet2_star);
    flx[2] += coeff[3]*(alpha_s*v2 + qf*bet2_star);
    flx[2] += coeff[4]*bet3;
    flx[2] += coeff[5]*(alpha_f*v2 - qs*bet2_star);

    flx[3] += coeff[0]*(alpha_f*v3 + qs*bet3_star);
    flx[3] += coeff[1]*bet2;
    flx[3] += coeff[2]*(alpha_s*v3 - qf*bet3_star);
    flx[3] += coeff[3]*(alpha_s*v3 + qf*bet3_star);
    flx[3] -= coeff[4]*bet2;
    flx[3] += coeff[5]*(alpha_f*v3 - qs*bet3_star);

    flx[4] += coeff[0]*as_prime*bet2_star;
    flx[4] -= coeff[1]*bet3*s/sqrtd;
    flx[4] -= coeff[2]*af_prime*bet2_star;
    flx[4] -= coeff[3]*af_prime*bet2_star;
    flx[4] -= coeff[4]*bet3*s/sqrtd;
    flx[4] += coeff[5]*as_prime*bet2_star;

    flx[5] += coeff[0]*as_prime*bet3_star;
    flx[5] += coeff[1]*bet2*s/sqrtd;
    flx[5] -= coeff[2]*af_prime*bet3_star;
    flx[5] -= coeff[3]*af_prime*bet3_star;
    flx[5] += coeff[4]*bet2*s/sqrtd;
    flx[5] += coeff[5]*as_prime*bet3_star;
  }
}
