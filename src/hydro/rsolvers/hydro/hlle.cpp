//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlle.cpp
//  \brief HLLE Riemann solver for hydrodynamics
//
//  Computes 1D fluxes using the Harten-Lax-van Leer (HLL) Riemann solver.  This flux is
//  very diffusive, especially for contacts, and so it is not recommended for use in
//  applications.  However, as shown by Einfeldt et al.(1991), it is positively
//  conservative (cannot return negative densities or pressure), so it is a useful
//  option when other approximate solvers fail and/or when extra dissipation is needed.
//
// REFERENCES:
// - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//   Springer-Verlag, Berlin, (1999) chpt. 10.
// - Einfeldt et al., "On Godunov-type methods near low densities", JCP, 92, 273 (1991)
// - A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and Godunov-type
//   schemes for hyperbolic conservation laws", SIAM Review 25, 35-61 (1983).

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"
#include "../../hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//  \brief The HLLE Riemann solver for hydrodynamics (both adiabatic and isothermal)

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          const AthenaArray<Real> &dxw) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NHYDRO)],wri[(NHYDRO)],wroe[(NHYDRO)];
  Real fl[(NHYDRO)],fr[(NHYDRO)],flxi[(NHYDRO)];
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
  Real gamma;
  if (GENERAL_EOS) {
    gamma = std::nan("");
  } else {
    gamma = pmy_block->peos->GetGamma();
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;

#pragma omp simd private(wli,wri,wroe,fl,fr,flxi)
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

    //--- Step 2.  Compute middle state estimates with PVRS (Toro 10.5.2)
    Real al, ar, el, er;
    Real cl = pmy_block->peos->SoundSpeed(wli);
    Real cr = pmy_block->peos->SoundSpeed(wri);
    if (NON_BAROTROPIC_EOS) {
      if (GENERAL_EOS) {
        el = pmy_block->peos->EgasFromRhoP(wli[IDN], wli[IPR]) +
             0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
        er = pmy_block->peos->EgasFromRhoP(wri[IDN], wri[IPR]) +
             0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
      } else {
        el = wli[IPR]*igm1 + 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
        er = wri[IPR]*igm1 + 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
      }
      Real rhoa = .5 * (wli[IDN] + wri[IDN]); // average density
      Real ca = .5 * (cl + cr); // average sound speed
      Real pmid = .5 * (wli[IPR] + wri[IPR] + (wli[IVX]-wri[IVX]) * rhoa * ca);
      Real umid = .5 * (wli[IVX] + wri[IVX] + (wli[IPR]-wri[IPR]) / (rhoa * ca));
      Real rhol = wli[IDN] + (wli[IVX] - umid) * rhoa / ca; // mid-left density
      Real rhor = wri[IDN] + (umid - wri[IVX]) * rhoa / ca; // mid-right density

      //--- Step 3.  Compute sound speed in L,R
      Real ql, qr;
      if (GENERAL_EOS) {
        Real gl = pmy_block->peos->AsqFromRhoP(rhol, pmid) * rhol / pmid;
        Real gr = pmy_block->peos->AsqFromRhoP(rhor, pmid) * rhor / pmid;
        ql = (pmid <= wli[IPR]) ? 1.0 :
             (1.0 + (gl + 1) / std::sqrt(2 * gl) * (pmid / wli[IPR]-1.0));
        qr = (pmid <= wri[IPR]) ? 1.0 :
             (1.0 + (gr + 1) / std::sqrt(2 * gr) * (pmid / wri[IPR]-1.0));
      } else {
        ql = (pmid <= wli[IPR]) ? 1.0 :
             (1.0 + (gamma + 1) / std::sqrt(2 * gamma) * (pmid / wli[IPR]-1.0));
        qr = (pmid <= wri[IPR]) ? 1.0 :
             (1.0 + (gamma + 1) / std::sqrt(2 * gamma) * (pmid / wri[IPR]-1.0));
      }

      //--- Step 4. Compute the max/min wave speeds based on L/R states

      al = wli[IVX] - cl*ql;
      ar = wri[IVX] + cr*qr;
    } else { // isothermal; revert to Roe average
      Real a  = iso_cs;
      Real sqrtdl = std::sqrt(wli[IDN]);
      Real sqrtdr = std::sqrt(wri[IDN]);
      Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

      wroe[IDN] = sqrtdl*sqrtdr;
      wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
      //--- Compute the max/min wave speeds based on L/R and Roe-averaged values
      Real al = std::min((wroe[IVX] - a),(wli[IVX] - cl));
      Real ar = std::max((wroe[IVX] + a),(wri[IVX] + cr));
    }

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

    //-- Step 5. Compute L/R fluxes along lines bm/bp: F_L - (S_L)U_L; F_R - (S_R)U_R
    Real vxl = wli[IVX] - bm;
    Real vxr = wri[IVX] - bp;

    fl[IDN] = wli[IDN]*vxl;
    fr[IDN] = wri[IDN]*vxr;

    fl[IVX] = wli[IDN]*wli[IVX]*vxl;
    fr[IVX] = wri[IDN]*wri[IVX]*vxr;

    fl[IVY] = wli[IDN]*wli[IVY]*vxl;
    fr[IVY] = wri[IDN]*wri[IVY]*vxr;

    fl[IVZ] = wli[IDN]*wli[IVZ]*vxl;
    fr[IVZ] = wri[IDN]*wri[IVZ]*vxr;

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IPR];
      fr[IVX] += wri[IPR];
      fl[IEN] = el*vxl + wli[IPR]*wli[IVX];
      fr[IEN] = er*vxr + wri[IPR]*wri[IVX];
    } else {
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

    //--- Step 6. Compute the HLLE flux at interface.
    Real tmp=0.0;
    if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);

    flxi[IDN] = 0.5*(fl[IDN]+fr[IDN]) + (fl[IDN]-fr[IDN])*tmp;
    flxi[IVX] = 0.5*(fl[IVX]+fr[IVX]) + (fl[IVX]-fr[IVX])*tmp;
    flxi[IVY] = 0.5*(fl[IVY]+fr[IVY]) + (fl[IVY]-fr[IVY])*tmp;
    flxi[IVZ] = 0.5*(fl[IVZ]+fr[IVZ]) + (fl[IVZ]-fr[IVZ])*tmp;
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN]+fr[IEN]) + (fl[IEN]-fr[IEN])*tmp;

    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,k,j,i) = flxi[IEN];
  }
  return;
}
