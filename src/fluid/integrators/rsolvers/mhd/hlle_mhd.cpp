//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 

//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
#include <iostream>

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../../../athena.hpp"         // enums, macros, Real
#include "../../../../athena_arrays.hpp"  // AthenaArray
#include "../../../fluid.hpp"             // Fluid
#include "../../../eos/eos.hpp"           // GetGamma

//======================================================================================
//! \file hlle_mhd.cpp
//  \brief HLLE Riemann solver for MHD.  See the hydro version for details.
//======================================================================================

void FluidIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl, 
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[NFLUID+NFIELD-1],wri[NFLUID+NFIELD-1],wroe[NFLUID+NFIELD-1];
  Real  fl[NFLUID+NFIELD-1], fr[NFLUID+NFIELD-1],flxi[NFLUID+NFIELD-1];
  Real gm1 = pmy_fluid->pf_eos->GetGamma() - 1.0;

#pragma simd
  for (int i=il; i<=iu; ++i){

//--- Step 1.  Load L/R states into local variables

    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    wli[IDN]=wl(IDN,i);
    if (NON_BAROTROPIC_EOS) wli[IEN]=wl(IEN,i);
    wli[IBY]=wl(IBY,i);
    wli[IBZ]=wl(IBZ,i);

    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    wri[IDN]=wr(IDN,i);
    if (NON_BAROTROPIC_EOS) wri[IEN]=wr(IEN,i);
    wri[IBY]=wr(IBY,i);
    wri[IBZ]=wr(IBZ,i);

    Real bxi = bx(k,j,i);

//--- Step 2.  Compute Roe-averaged state

    Real sqrtdl = sqrt(wli[IDN]);
    Real sqrtdr = sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN] = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;
// Note Roe average of magnetic field is different
    wroe[IBY] = (sqrtdr*wli[IBY] + sqrtdl*wri[IBY])*isdlpdr;
    wroe[IBZ] = (sqrtdr*wli[IBZ] + sqrtdl*wri[IBZ])*isdlpdr;
    Real pbl = 0.5*(bxi*bxi + SQR(wli[IBY]) + SQR(wli[IBZ]));
    Real pbr = 0.5*(bxi*bxi + wri[IBY]*wri[IBY] + wri[IBZ]*wri[IBZ]);
    Real x = 0.5*(SQR(wli[IBY]-wri[IBY]) + SQR(wli[IBZ]-wri[IBZ]))/(SQR(sqrtdl+sqrtdr));
    Real y = 0.5*(wli[IDN] + wri[IDN])/wroe[IDN];

// Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
// rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl

    Real el,er,hroe;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IEN]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX])+SQR(wli[IVY])+SQR(wli[IVZ])) +pbl;
      er = wri[IEN]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX])+SQR(wri[IVY])+SQR(wri[IVZ])) +pbr;
      hroe = ((el + wli[IEN] + pbl)/sqrtdl + (er + wri[IEN] + pbr)/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute fast magnetosonic speed in L,R, and Roe-averaged states

    Real cl = pmy_fluid->pf_eos->FastMagnetosonicSpeed(wli,bxi);
    Real cr = pmy_fluid->pf_eos->FastMagnetosonicSpeed(wri,bxi);

// Compute fast-magnetosonic speed using eq. B18
    Real a;
    if (NON_BAROTROPIC_EOS) {
      Real btsq = SQR(wroe[IBY]) + SQR(wroe[IBZ]);
      Real bt_starsq = (gm1 - (gm1 - 1.0)*y)*btsq;
      Real vaxsq = bxi*bxi/wroe[IDN];
      Real hp = hroe - (vaxsq + btsq/wroe[IDN]);
      Real vsq = SQR(wroe[IVX]) + SQR(wroe[IVY]) + SQR(wroe[IVZ]);
      Real twid_asq = std::max((gm1*(hp-0.5*vsq)-(gm1-1.0)*x), 0.0);
      Real ct2 = bt_starsq/wroe[IDN];
      Real tsum = vaxsq + ct2 + twid_asq;
      Real tdif = vaxsq + ct2 - twid_asq;
      Real cf2_cs2 = sqrt(tdif*tdif + 4.0*twid_asq*ct2);

      Real cfsq = 0.5*(tsum + cf2_cs2);
      a = sqrt(cfsq);
    } else {
      a = pmy_fluid->pf_eos->FastMagnetosonicSpeed(wroe,bxi);
    }

//--- Step 4.  Compute the max/min wave speeds based on L/R and Roe-averaged values

    Real al = std::min((wroe[IVX] - a),(wli[IVX] - cl));
    Real ar = std::max((wroe[IVX] + a),(wri[IVX] + cr));

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

//--- Step 5.  Compute L/R fluxes along the lines bm/bp: F_L - (S_L)U_L; F_R - (S_R)U_R

    fl[IDN] = wli[IDN]*wli[IVX] - bm*wli[IDN];
    fr[IDN] = wri[IDN]*wri[IVX] - bp*wri[IDN];

    fl[IVX] = wli[IDN]*wli[IVX]*(wli[IVX] - bm);
    fr[IVX] = wri[IDN]*wri[IVX]*(wri[IVX] - bp);

    fl[IVY] = wli[IDN]*wli[IVY]*(wli[IVX] - bm);
    fr[IVY] = wri[IDN]*wri[IVY]*(wri[IVX] - bp);

    fl[IVZ] = wli[IDN]*wli[IVZ]*(wli[IVX] - bm);
    fr[IVZ] = wri[IDN]*wri[IVZ]*(wri[IVX] - bp);

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IEN];
      fr[IVX] += wri[IEN];
      fl[IEN] = el*(wli[IVX] - bm) + wli[IEN]*wli[IVX];
      fr[IEN] = er*(wri[IVX] - bp) + wri[IEN]*wri[IVX];
    } else {
      Real iso_cs = pmy_fluid->pf_eos->SoundSpeed(wli);
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

// Add MHD terms

    fl[IVX] -= 0.5*(bxi*bxi - SQR(wli[IBY]) - SQR(wli[IBZ]));
    fr[IVX] -= 0.5*(bxi*bxi - SQR(wri[IBY]) - SQR(wri[IBZ]));

    fl[IVY] -= bxi*wli[IBY];
    fr[IVY] -= bxi*wri[IBY];

    fl[IVZ] -= bxi*wli[IBZ];
    fr[IVZ] -= bxi*wri[IBZ];

    if (NON_BAROTROPIC_EOS) {
      fl[IEN] +=(pbl*wli[IVX] - bxi*(bxi*wli[IVX]+wli[IBY]*wli[IVY]+wli[IBZ]*wli[IVZ]));
      fr[IEN] +=(pbr*wri[IVX] - bxi*(bxi*wri[IVX]+wri[IBY]*wri[IVY]+wri[IBZ]*wri[IVZ]));
    }

    fl[IBY] = wli[IBY]*(wli[IVX] - bm) - bxi*wli[IVY];;
    fr[IBY] = wri[IBY]*(wri[IVX] - bp) - bxi*wri[IVY];;

    fl[IBZ] = wli[IBZ]*(wli[IVX] - bm) - bxi*wli[IVZ];;
    fr[IBZ] = wri[IBZ]*(wri[IVX] - bp) - bxi*wri[IVZ];;

//--- Step 6.  Compute the HLLE flux at interface.

    Real tmp;
    if (bp == bm) {
      tmp = 0.0;
    } else {
      tmp = 0.5*(bp + bm)/(bp - bm);
    }

    flxi[IDN] = 0.5*(fl[IDN]+fr[IDN]) + (fl[IDN]-fr[IDN])*tmp;
    flxi[IVX] = 0.5*(fl[IVX]+fr[IVX]) + (fl[IVX]-fr[IVX])*tmp;
    flxi[IVY] = 0.5*(fl[IVY]+fr[IVY]) + (fl[IVY]-fr[IVY])*tmp;
    flxi[IVZ] = 0.5*(fl[IVZ]+fr[IVZ]) + (fl[IVZ]-fr[IVZ])*tmp;
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN]+fr[IEN]) + (fl[IEN]-fr[IEN])*tmp;
    flxi[IBY] = 0.5*(fl[IBY]+fr[IBY]) + (fl[IBY]-fr[IBY])*tmp;
    flxi[IBZ] = 0.5*(fl[IBZ]+fr[IBZ]) + (fl[IBZ]-fr[IBZ])*tmp;

    flx(IDN,i) = flxi[IDN];
    flx(ivx,i) = flxi[IVX];
    flx(ivy,i) = flxi[IVY];
    flx(ivz,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,i) = flxi[IEN];
    flx(IBY,i) = flxi[IBY];
    flx(IBZ,i) = flxi[IBZ];
  }

  return;
}
