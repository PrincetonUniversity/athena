//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file llf_mhd.cpp
//  \brief Local Lax Friedrichs (LLF) Riemann solver for MHD
//
//  Computes 1D fluxes using the LLF Riemann solver, also known as Rusanov's method.
//  This flux is very diffusive, even more diffusive than HLLE, and so it is not
//  recommended for use in applications.  However, it is useful for testing, or for
//  problems where other Riemann solvers fail.
//
// REFERENCES:
// - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//   Springer-Verlag, Berlin, (1999) chpt. 10.

// C/C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//  \brief The LLF Riemann solver for MHD (both adiabatic and isothermal)

void Hydro::RiemannSolver(const int kl, const int ku, const int jl, const int ju,
  const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
  AthenaArray<Real> &ey, AthenaArray<Real> &ez) {

  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NWAVE)],wri[(NWAVE)],du[(NWAVE)];
  Real flxi[(NWAVE)],fl[(NWAVE)],fr[(NWAVE)];
  Real gm1 = pmy_block->peos->GetGamma() - 1.0;
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
#pragma omp simd private(wli,wri,du,fl,fr,flxi)
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

//--- Step 2.  Compute wave speeds in L,R states (see Toro eq. 10.43)

    Real cfl = pmy_block->peos->FastMagnetosonicSpeed(wli,bxi);
    Real cfr = pmy_block->peos->FastMagnetosonicSpeed(wri,bxi);
    Real a = 0.5*std::max( (fabs(wli[IVX]) + cfl), (fabs(wri[IVX]) + cfr) );

//--- Step 3.  Compute L/R fluxes

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];
    Real pbl = 0.5*(bxi*bxi + SQR(wli[IBY]) + SQR(wli[IBZ]));
    Real pbr = 0.5*(bxi*bxi + SQR(wri[IBY]) + SQR(wri[IBZ]));

    fl[IDN] = mxl;
    fr[IDN] = mxr;

    fl[IVX] = mxl*wli[IVX] + pbl - SQR(bxi);
    fr[IVX] = mxr*wri[IVX] + pbr - SQR(bxi);

    fl[IVY] = mxl*wli[IVY] - bxi*wli[IBY];
    fr[IVY] = mxr*wri[IVY] - bxi*wri[IBY];

    fl[IVZ] = mxl*wli[IVZ] - bxi*wli[IBZ];
    fr[IVZ] = mxr*wri[IVZ] - bxi*wri[IBZ];

    Real el,er;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IPR]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX])+SQR(wli[IVY])+SQR(wli[IVZ])) + pbl;
      er = wri[IPR]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX])+SQR(wri[IVY])+SQR(wri[IVZ])) + pbr;
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

//--- Step 4.  Compute difference in L/R states dU

    du[IDN] = wri[IDN]          - wli[IDN];
    du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
    du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
    du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) du[IEN] = er - el;
    du[IBY] = wri[IBY] - wli[IBY];
    du[IBZ] = wri[IBZ] - wli[IBZ];

//--- Step 5.  Compute the LLF flux at interface (see Toro eq. 10.42).

    flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]) - a*du[IDN];
    flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]) - a*du[IVX];
    flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]) - a*du[IVY];
    flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]) - a*du[IVZ];
    if (NON_BAROTROPIC_EOS) {
      flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]) - a*du[IEN];
    }
    flxi[IBY] = 0.5*(fl[IBY] + fr[IBY]) - a*du[IBY];
    flxi[IBZ] = 0.5*(fl[IBZ] + fr[IBZ]) - a*du[IBZ];

//--- Step 6. Store results into 3D array of fluxes

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
