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
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

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
//! \file hllc.cpp
//  \brief HLLC Riemann solver for hydrodynamics, an extension of  the HLLE fluxes to
//    include the contact wave.  Only works for adiabatic hydrodynamics.
//
// REFERENCES:
// - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//   Springer-Verlag, Berlin, (1999) chpt. 10.
//
// - P. Batten, N. Clarke, C. Lambert, and D. M. Causon, "On the Choice of Wavespeeds
//   for the HLLC Riemann Solver", SIAM J. Sci. & Stat. Comp. 18, 6, 1553-1570, (1997).
//======================================================================================

void FluidIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[(NFLUID)],wri[(NFLUID)],wroe[(NFLUID)];
  Real flxi[(NFLUID)],fl[(NFLUID)],fr[(NFLUID)];
  Real gm1 = pmy_fluid->pf_eos->GetGamma() - 1.0;

#pragma simd
  for (int i=il; i<=iu; ++i){

//--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    wli[IEN]=wl(IEN,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    wri[IEN]=wr(IEN,i);

//--- Step2.  Compute Roe-averaged state

    Real sqrtdl = sqrt(wli[IDN]);
    Real sqrtdr = sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN] = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;

    // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
    // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    Real el,er,hroe;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IEN]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
      er = wri[IEN]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
      hroe = ((el + wli[IEN])/sqrtdl + (er + wri[IEN])/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute sound speed in L,R, and Roe-averaged states

    Real cl = pmy_fluid->pf_eos->SoundSpeed(wli);
    Real cr = pmy_fluid->pf_eos->SoundSpeed(wri);
    Real q = hroe - 0.5*(wroe[IVX]*wroe[IVX]+wroe[IVY]*wroe[IVY]+wroe[IVZ]*wroe[IVZ]);
    if (q < 0.0) q=0.0;
    Real a = sqrt(gm1*q);

//--- Step 4.  Compute the max/min wave speeds based on L/R and Roe-averaged values

    Real al = std::min((wroe[IVX] - a),(wli[IVX] - cl));
    Real ar = std::max((wroe[IVX] + a),(wri[IVX] + cr));

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

//--- Step 5.  Compute the contact wave speed and pressure

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];

    Real tl = wli[IEN] + (wli[IVX] - al)*mxl;
    Real tr = wri[IEN] + (wri[IVX] - ar)*mxr;

    Real ml =   mxl - wli[IDN]*al;
    Real mr = -(mxr - wri[IDN]*ar);

    // Determine the contact wave speed...
    Real am = (tl - tr)/(ml + mr);
    // ...and the pressure at the contact surface
    Real cp = (ml*tr + mr*tl)/(ml + mr);
    cp = cp > 0.0 ? cp : 0.0;

//--- Step 6.  Compute L/R fluxes along the line bm, bp

    fl[IDN] = mxl - bm*wli[IDN];
    fr[IDN] = mxr - bp*wri[IDN];

    fl[IVX] = mxl*(wli[IVX] - bm) + wli[IEN];
    fr[IVX] = mxr*(wri[IVX] - bp) + wri[IEN];

    fl[IVY] = wli[IDN]*wli[IVY]*(wli[IVX] - bm);
    fr[IVY] = wri[IDN]*wri[IVY]*(wri[IVX] - bp);

    fl[IVZ] = wli[IDN]*wli[IVZ]*(wli[IVX] - bm);
    fr[IVZ] = wri[IDN]*wri[IVZ]*(wri[IVX] - bp);

    fl[IEN] = wli[IEN]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX])+SQR(wli[IVY])+SQR(wli[IVZ]));
    fr[IEN] = wri[IEN]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX])+SQR(wri[IVY])+SQR(wri[IVZ]));
    fl[IEN] *= (wli[IVX] - bm);
    fr[IEN] *= (wri[IVX] - bp);
    fl[IEN] += wli[IEN]*wli[IVX];
    fr[IEN] += wri[IEN]*wri[IVX];

//--- Step 8.  Compute flux weights or scales

    Real sl,sr,sm;
    if (am >= 0.0) {
      sl =  am/(am - bm);
      sr = 0.0;
      sm = -bm/(am - bm);
    }
    else {
      sl =  0.0;
      sr = -am/(bp - am);
      sm =  bp/(bp - am);
    }

//--- Step 9.  Compute the HLLC flux at interface, including the weighted contribution
// of the flux along the contact

    flxi[IDN] = sl*fl[IDN] + sr*fr[IDN];
    flxi[IVX] = sl*fl[IVX] + sr*fr[IVX] + sm*cp;
    flxi[IVY] = sl*fl[IVY] + sr*fr[IVY];
    flxi[IVZ] = sl*fl[IVZ] + sr*fr[IVZ];
    flxi[IEN] = sl*fl[IEN] + sr*fr[IEN] + sm*cp*am;

    flx(IDN,i) = flxi[IDN];
    flx(ivx,i) = flxi[IVX];
    flx(ivy,i) = flxi[IVY];
    flx(ivz,i) = flxi[IVZ];
    flx(IEN,i) = flxi[IEN];
  }

  return;
}
