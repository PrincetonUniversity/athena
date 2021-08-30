//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file lhllc.cpp
//! \brief Low-dissipation HLLC (LHLLC) Riemann solver for adiabatic hydrodynamics.
//!        (Minoshima et al. 2021)

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
//! \brief The Low-dissipation HLLC (LHLLC) Riemann solver for adiabatic hydrodynamics

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          const AthenaArray<Real> &dxw) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NHYDRO)],wri[(NHYDRO)];
  Real flxi[(NHYDRO)],fl[(NHYDRO)],fr[(NHYDRO)];
  Real gamma;
  if (GENERAL_EOS) {
    gamma = std::nan("");
  } else {
    gamma = pmy_block->peos->GetGamma();
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;

  CalculateVelocityDifferences(k, j, il, iu, ivx, dvn, dvt);

#pragma distribute_point
#pragma omp simd private(wli,wri,flxi,fl,fr)
  for (int i=il; i<=iu; ++i) {
    //--- Step 1.  Load L/R states into local variables
    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    wli[IPR]=wl(IPR,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    wri[IPR]=wr(IPR,i);

    //--- Step 2.  Compute middle state estimates with PVRS (Toro 10.5.2)

    Real al, ar, el, er;
    Real cl = pmy_block->peos->SoundSpeed(wli);
    Real cr = pmy_block->peos->SoundSpeed(wri);
    Real vsql = SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]);
    Real vsqr = SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]);
    if (GENERAL_EOS) {
      el = pmy_block->peos->EgasFromRhoP(wli[IDN], wli[IPR]) + 0.5*wli[IDN]*vsql;
      er = pmy_block->peos->EgasFromRhoP(wri[IDN], wri[IPR]) + 0.5*wri[IDN]*vsqr;
    } else {
      el = wli[IPR]*igm1 + 0.5*wli[IDN]*vsql;
      er = wri[IPR]*igm1 + 0.5*wri[IDN]*vsqr;
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
           std::sqrt(1.0 + (gl + 1) / (2 * gl) * (pmid / wli[IPR]-1.0));
      qr = (pmid <= wri[IPR]) ? 1.0 :
           std::sqrt(1.0 + (gr + 1) / (2 * gr) * (pmid / wri[IPR]-1.0));
    } else {
      ql = (pmid <= wli[IPR]) ? 1.0 :
           std::sqrt(1.0 + (gamma + 1) / (2 * gamma) * (pmid / wli[IPR]-1.0));
      qr = (pmid <= wri[IPR]) ? 1.0 :
           std::sqrt(1.0 + (gamma + 1) / (2 * gamma) * (pmid / wri[IPR]-1.0));
    }

    //--- Step 4.  Compute the max/min wave speeds based on L/R

    al = wli[IVX] - cl*ql;
    ar = wri[IVX] + cr*qr;

    Real bp = ar > 0.0 ? ar : (TINY_NUMBER);
    Real bm = al < 0.0 ? al : -(TINY_NUMBER);

    //--- Step 5. Compute the contact wave speed and pressure

    Real vxl = al - wli[IVX];
    Real vxr = ar - wri[IVX];

    Real ml = wli[IDN]*vxl;
    Real mr = wri[IDN]*vxr;

    // shock detector
    // Real cmax = std::max(cl*ql, cr*qr);
    Real cmax = std::max(cl, cr);
    Real th1 = std::min(1.0, (cmax-std::min(dvn(i),0.0)) / (cmax-std::min(dvt(i),0.0)));
    Real th = th1 * th1 * th1 * th1; // this 4th power is empirical (see Minoshima+)

    // Determine the contact wave speed...
    Real am = (mr*wri[IVX] - ml*wli[IVX] - th*(wri[IPR]-wli[IPR])) / (mr - ml);

    // ...and the pressure at the contact surface
    Real chi = std::min(1.0, std::sqrt(std::max(vsql, vsqr)) / cmax);
    Real phi = chi * (2.0 - chi);

    Real cp = (mr*wli[IPR] - ml*wri[IPR] + phi*mr*ml*(wri[IVX]-wli[IVX])) / (mr - ml);
    cp = cp > 0.0 ? cp : 0.0;

    // No loop-carried dependencies anywhere in this loop
    //    #pragma distribute_point
    //--- Step 6. Compute L/R fluxes along the line bm, bp

    Real uxl = wli[IVX] - bm;
    Real uxr = wri[IVX] - bp;

    fl[IDN] = wli[IDN]*uxl;
    fr[IDN] = wri[IDN]*uxr;

    fl[IVX] = wli[IDN]*wli[IVX]*uxl + wli[IPR];
    fr[IVX] = wri[IDN]*wri[IVX]*uxr + wri[IPR];

    fl[IVY] = wli[IDN]*wli[IVY]*uxl;
    fr[IVY] = wri[IDN]*wri[IVY]*uxr;

    fl[IVZ] = wli[IDN]*wli[IVZ]*uxl;
    fr[IVZ] = wri[IDN]*wri[IVZ]*uxr;

    fl[IEN] = el*uxl + wli[IPR]*wli[IVX];
    fr[IEN] = er*uxr + wri[IPR]*wri[IVX];

    //--- Step 8. Compute flux weights or scales

    Real sl,sr,sm;
    if (am >= 0.0) {
      sl =  am/(am - bm);
      sr = 0.0;
      sm = -bm/(am - bm);
    } else {
      sl =  0.0;
      sr = -am/(bp - am);
      sm =  bp/(bp - am);
    }

    //--- Step 9. Compute the HLLC flux at interface, including weighted contribution
    // of the flux along the contact

    flxi[IDN] = sl*fl[IDN] + sr*fr[IDN];
    flxi[IVX] = sl*fl[IVX] + sr*fr[IVX] + sm*cp;
    flxi[IVY] = sl*fl[IVY] + sr*fr[IVY];
    flxi[IVZ] = sl*fl[IVZ] + sr*fr[IVZ];
    flxi[IEN] = sl*fl[IEN] + sr*fr[IEN] + sm*cp*am;

    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    flx(IEN,k,j,i) = flxi[IEN];
  }
  return;
}
