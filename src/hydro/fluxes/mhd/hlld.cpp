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
//! \file hlld.cpp
//  \brief HLLD Riemann solver for adiabatic MHD.
//
// REFERENCES:
// - T. Miyoshi & K. Kusano, "A multi-state HLL approximate Riemann solver for ideal
//   MHD", JCP, 208, 315 (2005)
//======================================================================================

// C++ headers
#include <algorithm>  // max(), min()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../hydro.hpp"
#include "../../../eos/eos.hpp"

// this class header
#include "../hydro_fluxes.hpp"

// container to store (density, momentum, total energy, tranverse magnetic field)
// minimizes changes required to adopt athena4.2 version of this solver 
typedef struct Cons1D {
  Real d,mx,my,mz,e,by,bz;
} Cons1D;

#define SMALL_NUMBER 1.0e-8

void HydroFluxes::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real flxi[(NWAVE)];             // temporary variable to store flux
  Real wli[(NWAVE)],wri[(NWAVE)]; // L/R states, primitive variables (input)
  Cons1D ul,ur;                   // L/R states, conserved variables (computed)
  Real spd[5];                    // signal speeds, left to right 
  Cons1D ulst,uldst,urdst,urst;   // Conserved variable for all states 
  Cons1D fl,fr;                   // Fluxes for left & right states

  Real gm1 = pmy_hydro->pmy_block->peos->GetGamma() - 1.0;

#pragma simd
  for (int i=il; i<=iu; ++i){

//--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    wli[IEN]=wl(IEN,i);
    wli[IBY]=wl(IBY,i);
    wli[IBZ]=wl(IBZ,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    wri[IEN]=wr(IEN,i);
    wri[IBY]=wr(IBY,i);
    wri[IBZ]=wr(IBZ,i);

    Real bxi = bx(k,j,i);

    // Compute L/R states for selected conserved variables
    Real bxsq = bxi*bxi;
    Real pbl = 0.5*(bxsq + SQR(wli[IBY]) + SQR(wli[IBZ]));  // magnetic pressure (l/r)
    Real pbr = 0.5*(bxsq + SQR(wri[IBY]) + SQR(wri[IBZ]));
    Real kel = 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
    Real ker = 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));

    ul.d  = wli[IDN];
    ul.mx = wli[IVX]*ul.d;
    ul.my = wli[IVY]*ul.d;
    ul.mz = wli[IVZ]*ul.d;
    ul.e  = wli[IEN]/gm1 + kel + pbl;
    ul.by = wli[IBY];
    ul.bz = wli[IBZ];

    ur.d  = wri[IDN];
    ur.mx = wri[IVX]*ur.d;
    ur.my = wri[IVY]*ur.d;
    ur.mz = wri[IVZ]*ur.d;
    ur.e  = wri[IEN]/gm1 + ker + pbr;
    ur.by = wri[IBY];
    ur.bz = wri[IBZ];

//--- Step 2.  Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)

    Real cfl = pmy_hydro->pmy_block->peos->FastMagnetosonicSpeed(wli,bxi);
    Real cfr = pmy_hydro->pmy_block->peos->FastMagnetosonicSpeed(wri,bxi);

    spd[0] = std::min( wli[IVX]-cfl, wri[IVX]-cfr );
    spd[4] = std::max( wli[IVX]+cfl, wri[IVX]+cfr );
/*
    Real cfmax = std::max(cfl,cfr);
    if(wli[IVX] <= wri[IVX]) {
      spd[0] = wli[IVX] - cfmax;
      spd[4] = wri[IVX] + cfmax;
    } else {
      spd[0] = wri[IVX] - cfmax;
      spd[4] = wli[IVX] + cfmax;
    }
*/

//--- Step 3.  Compute L/R fluxes

    Real ptl = wli[IEN] + pbl; // total pressures L,R
    Real ptr = wri[IEN] + pbr;

    fl.d  = ul.mx;
    fl.mx = ul.mx*wli[IVX] + ptl - bxsq;
    fl.my = ul.my*wli[IVX] - bxi*ul.by;
    fl.mz = ul.mz*wli[IVX] - bxi*ul.bz;
    fl.e  = wli[IVX]*(ul.e + ptl - bxsq) - bxi*(wli[IVY]*ul.by + wli[IVZ]*ul.bz);
    fl.by = ul.by*wli[IVX] - bxi*wli[IVY];
    fl.bz = ul.bz*wli[IVX] - bxi*wli[IVZ];

    fr.d  = ur.mx;
    fr.mx = ur.mx*wri[IVX] + ptr - bxsq;
    fr.my = ur.my*wri[IVX] - bxi*ur.by;
    fr.mz = ur.mz*wri[IVX] - bxi*ur.bz;
    fr.e  = wri[IVX]*(ur.e + ptr - bxsq) - bxi*(wri[IVY]*ur.by + wri[IVZ]*ur.bz);
    fr.by = ur.by*wri[IVX] - bxi*wri[IVY];
    fr.bz = ur.bz*wri[IVX] - bxi*wri[IVZ];

//--- Step 4.  Compute middle and Alfven wave speeds

    Real sdl = spd[0] - wli[IVX];  // S_i-u_i (i=L or R)
    Real sdr = spd[4] - wri[IVX];

    // S_M: eqn (38) of Miyoshi & Kusano 
    spd[2] = (sdr*ur.mx - sdl*ul.mx - ptr + ptl)/(sdr*ur.d - sdl*ul.d);

    Real sdml   = spd[0] - spd[2];  // S_i-S_M (i=L or R)
    Real sdmr   = spd[4] - spd[2];
    // eqn (43) of Miyoshi & Kusano
    ulst.d = ul.d * sdl/sdml;
    urst.d = ur.d * sdr/sdmr;
    Real sqrtdl = sqrt(ulst.d);
    Real sqrtdr = sqrt(urst.d);

    // eqn (51) of Miyoshi & Kusano
    spd[1] = spd[2] - fabs(bxi)/sqrtdl;
    spd[3] = spd[2] + fabs(bxi)/sqrtdr;

//--- Step 5.  Compute intermediate states

    Real ptst = ptl + ul.d*sdl*(sdl-sdml);  // total pressure (star state)
 
  // ul* - eqn (39) of M&K
    ulst.mx = ulst.d * spd[2];
    if (fabs(ul.d*sdl*sdml-bxsq) < (SMALL_NUMBER)*ptst) {
      // Degenerate case
      ulst.my = ulst.d * wli[IVY];
      ulst.mz = ulst.d * wli[IVZ];

      ulst.by = ul.by;
      ulst.bz = ul.bz;
    } else {
      // eqns (44) and (46) of M&K
      Real tmp = bxi*(sdl - sdml)/(ul.d*sdl*sdml - bxsq);
      ulst.my = ulst.d * (wli[IVY] - ul.by*tmp);
      ulst.mz = ulst.d * (wli[IVZ] - ul.bz*tmp);

      // eqns (45) and (47) of M&K
      tmp = (ul.d*SQR(sdl) - bxsq)/(ul.d*sdl*sdml - bxsq);
      ulst.by = ul.by * tmp;
      ulst.bz = ul.bz * tmp;
    }
    Real vbstl=(ulst.mx*bxi+ulst.my*ulst.by+ulst.mz*ulst.bz)/ulst.d; // v_i* dot B_i*
    // eqn (48) of M&K
    ulst.e = (sdl*ul.e - ptl*wli[IVX] + ptst*spd[2] +
              bxi*(wli[IVX]*bxi + wli[IVY]*ul.by + wli[IVZ]*ul.bz - vbstl))/sdml;

  // ur* - eqn (39) of M&K
    urst.mx = urst.d * spd[2];
    if (fabs(ur.d*sdr*sdmr - bxsq) < (SMALL_NUMBER)*ptst) {
      // Degenerate case
      urst.my = urst.d * wri[IVY];
      urst.mz = urst.d * wri[IVZ];

      urst.by = ur.by;
      urst.bz = ur.bz;
    } else {
      // eqns (44) and (46) of M&K 
      Real tmp = bxi*(sdr - sdmr)/(ur.d*sdr*sdmr - bxsq);
      urst.my = urst.d * (wri[IVY] - ur.by*tmp);
      urst.mz = urst.d * (wri[IVZ] - ur.bz*tmp);

      // eqns (45) and (47) of M&K
      tmp = (ur.d*SQR(sdr) - bxsq)/(ur.d*sdr*sdmr - bxsq);
      urst.by = ur.by * tmp;
      urst.bz = ur.bz * tmp;
    }
    Real vbstr=(urst.mx*bxi+urst.my*urst.by+urst.mz*urst.bz)/urst.d; // v_i* dot B_i*
    // eqn (48) of M&K
    urst.e = (sdr*ur.e - ptr*wri[IVX] + ptst*spd[2] +
              bxi*(wri[IVX]*bxi + wri[IVY]*ur.by + wri[IVZ]*ur.bz - vbstr))/sdmr;

  // ul** and ur** - if Bx is near zero, same as *-states
    if(0.5*bxsq < (SMALL_NUMBER)*ptst) {
      uldst = ulst;
      urdst = urst;
    } else {
      Real invsumd = 1.0/(sqrtdl + sqrtdr);
      Real bxsig = (bxi > 0.0 ? 1.0 : -1.0);

      uldst.d = ulst.d;
      urdst.d = urst.d;

      uldst.mx = ulst.mx;
      urdst.mx = urst.mx;

      // eqn (59) of M&K
      Real tmp = invsumd*(sqrtdl*(ulst.my/ulst.d) + sqrtdr*(urst.my/urst.d) +
                 bxsig*(urst.by - ulst.by));
      uldst.my = uldst.d * tmp;
      urdst.my = urdst.d * tmp;

      // eqn (60) of M&K
      tmp = invsumd*(sqrtdl*(ulst.mz/ulst.d) + sqrtdr*(urst.mz/urst.d) +
            bxsig*(urst.bz - ulst.bz));
      uldst.mz = uldst.d * tmp;
      urdst.mz = urdst.d * tmp;

      // eqn (61) of M&K
      tmp = invsumd*(sqrtdl*urst.by + sqrtdr*ulst.by +
                     bxsig*sqrtdl*sqrtdr*((urst.my/urst.d) - (ulst.my/ulst.d)));
      uldst.by = urdst.by = tmp;

      // eqn (62) of M&K
      tmp = invsumd*(sqrtdl*urst.bz + sqrtdr*ulst.bz +
                     bxsig*sqrtdl*sqrtdr*((urst.mz/urst.d) - (ulst.mz/ulst.d)));
      uldst.bz = urdst.bz = tmp;

      // eqn (63) of M&K
      tmp = spd[2]*bxi + (uldst.my*uldst.by + uldst.mz*uldst.bz)/uldst.d;
      uldst.e = ulst.e - sqrtdl*bxsig*(vbstl - tmp);
      urdst.e = urst.e + sqrtdr*bxsig*(vbstr - tmp);
    }

//--- Step 6.  Compute flux

    if(spd[0] >= 0.0){
      // return Fl if flow is supersonic
      flxi[IDN] = fl.d;
      flxi[IVX] = fl.mx;
      flxi[IVY] = fl.my;
      flxi[IVZ] = fl.mz;
      flxi[IEN] = fl.e;
      flxi[IBY] = fl.by;
      flxi[IBZ] = fl.bz;
    } else if(spd[4] <= 0.0){
      // return Fr if flow is supersonic
      flxi[IDN] = fr.d;
      flxi[IVX] = fr.mx;
      flxi[IVY] = fr.my;
      flxi[IVZ] = fr.mz;
      flxi[IEN] = fr.e;
      flxi[IBY] = fr.by;
      flxi[IBZ] = fr.bz;
    } else if(spd[1] >= 0.0) {
      // return Fl*
      flxi[IDN] = fl.d  + spd[0]*(ulst.d  - ul.d);
      flxi[IVX] = fl.mx + spd[0]*(ulst.mx - ul.mx);
      flxi[IVY] = fl.my + spd[0]*(ulst.my - ul.my);
      flxi[IVZ] = fl.mz + spd[0]*(ulst.mz - ul.mz);
      flxi[IEN] = fl.e  + spd[0]*(ulst.e  - ul.e);
      flxi[IBY] = fl.by + spd[0]*(ulst.by - ul.by);
      flxi[IBZ] = fl.bz + spd[0]*(ulst.bz - ul.bz);
    } else if(spd[2] >= 0.0) {
      // return Fl**
      Real tmp = spd[1] - spd[0];
      flxi[IDN] = fl.d  - spd[0]*ul.d  - tmp*ulst.d  + spd[1]*uldst.d;
      flxi[IVX] = fl.mx - spd[0]*ul.mx - tmp*ulst.mx + spd[1]*uldst.mx;
      flxi[IVY] = fl.my - spd[0]*ul.my - tmp*ulst.my + spd[1]*uldst.my;
      flxi[IVZ] = fl.mz - spd[0]*ul.mz - tmp*ulst.mz + spd[1]*uldst.mz;
      flxi[IEN] = fl.e  - spd[0]*ul.e  - tmp*ulst.e  + spd[1]*uldst.e;
      flxi[IBY] = fl.by - spd[0]*ul.by - tmp*ulst.by + spd[1]*uldst.by;
      flxi[IBZ] = fl.bz - spd[0]*ul.bz - tmp*ulst.bz + spd[1]*uldst.bz;
    } else if(spd[3] > 0.0) {
      // return Fr**
      Real tmp = spd[3] - spd[4];
      flxi[IDN] = fr.d  - spd[4]*ur.d  - tmp*urst.d  + spd[3]*urdst.d;
      flxi[IVX] = fr.mx - spd[4]*ur.mx - tmp*urst.mx + spd[3]*urdst.mx;
      flxi[IVY] = fr.my - spd[4]*ur.my - tmp*urst.my + spd[3]*urdst.my;
      flxi[IVZ] = fr.mz - spd[4]*ur.mz - tmp*urst.mz + spd[3]*urdst.mz;
      flxi[IEN] = fr.e  - spd[4]*ur.e  - tmp*urst.e  + spd[3]*urdst.e;
      flxi[IBY] = fr.by - spd[4]*ur.by - tmp*urst.by + spd[3]*urdst.by;
      flxi[IBZ] = fr.bz - spd[4]*ur.bz - tmp*urst.bz + spd[3]*urdst.bz;
    } else {
      // return Fr* 
      flxi[IDN] = fr.d  + spd[4]*(urst.d  - ur.d);
      flxi[IVX] = fr.mx + spd[4]*(urst.mx - ur.mx);
      flxi[IVY] = fr.my + spd[4]*(urst.my - ur.my);
      flxi[IVZ] = fr.mz + spd[4]*(urst.mz - ur.mz);
      flxi[IEN] = fr.e  + spd[4]*(urst.e  - ur.e);
      flxi[IBY] = fr.by + spd[4]*(urst.by - ur.by);
      flxi[IBZ] = fr.bz + spd[4]*(urst.bz - ur.bz);
    }

    flx(IDN,i) = flxi[IDN];
    flx(ivx,i) = flxi[IVX];
    flx(ivy,i) = flxi[IVY];
    flx(ivz,i) = flxi[IVZ];
    flx(IEN,i) = flxi[IEN];
    flx(IBY,i) = flxi[IBY];
    flx(IBZ,i) = flxi[IBZ];

  }      

  return;
}
