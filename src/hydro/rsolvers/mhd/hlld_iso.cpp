//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hlld_iso.cpp
//! \brief HLLD Riemann solver for isothermal MHD.
//!
//! REFERENCES:
//! - A. Mignone, "A simple and accurate Riemann solver for isothermal MHD", JPC, 225,
//!   1427 (2007)

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"
#include "../../hydro.hpp"

// container to store (density, momentum, tranverse magnetic field)
// minimizes changes required to adopt athena4.2 version of this solver
struct Cons1D {
  Real d, mx, my, mz, by, bz;
};

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//! \brief The HLLD Riemann solver for isothermal MHD

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx, const AthenaArray<Real> &bx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr,
                          AthenaArray<Real> &flx,
                          AthenaArray<Real> &ey, AthenaArray<Real> &ez,
                          AthenaArray<Real> &wct, const AthenaArray<Real> &dxw) {
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real flxi[(NWAVE)];             // temporary variable to store flux
  Real wli[(NWAVE)],wri[(NWAVE)]; // L/R states, primitive variables (input)
  Real spd[5];                    // signal speeds, left to right
  constexpr Real SMALL_NUMBER = 1.0e-4;

  Real dfloor = pmy_block->peos->GetDensityFloor();
  Real cs = (pmy_block->peos->GetIsoSoundSpeed());
  Real dt = pmy_block->pmy_mesh->dt;

#pragma omp simd private(flxi,wli,wri,spd)
  for (int i=il; i<=iu; ++i) {
    Cons1D ul,ur;                   // L/R states, conserved variables (computed)
    Cons1D ulst,urst,ucst;          // Conserved variable for all states
    Cons1D fl,fr;                   // Fluxes for left & right states

    //--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    wli[IBY]=wl(IBY,i);
    wli[IBZ]=wl(IBZ,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    wri[IBY]=wr(IBY,i);
    wri[IBZ]=wr(IBZ,i);

    Real bxi = bx(k,j,i);

    // Compute L/R states for selected conserved variables
    ul.d  = wli[IDN];
    ul.mx = wli[IVX]*ul.d;
    ul.my = wli[IVY]*ul.d;
    ul.mz = wli[IVZ]*ul.d;
    ul.by = wli[IBY];
    ul.bz = wli[IBZ];

    ur.d  = wri[IDN];
    ur.mx = wri[IVX]*ur.d;
    ur.my = wri[IVY]*ur.d;
    ur.mz = wri[IVZ]*ur.d;
    ur.by = wri[IBY];
    ur.bz = wri[IBZ];

    //--- Step 2.  Compute L & R wave speeds according to Miyoshi & Kusano, eqn. (67)

    Real cfl = pmy_block->peos->FastMagnetosonicSpeed(wli,bxi);
    Real cfr = pmy_block->peos->FastMagnetosonicSpeed(wri,bxi);

    spd[0] = std::min( wli[IVX]-cfl, wri[IVX]-cfr );
    spd[4] = std::max( wli[IVX]+cfl, wri[IVX]+cfr );

    //--- Step 3.  Compute L/R fluxes

    // total pressures L,R
    Real bxsq = bxi*bxi;
    Real ptl = SQR(cs)*wli[IDN] + 0.5*(bxsq + SQR(wli[IBY]) + SQR(wli[IBZ]));
    Real ptr = SQR(cs)*wri[IDN] + 0.5*(bxsq + SQR(wri[IBY]) + SQR(wri[IBZ]));

    fl.d  = ul.mx;
    fl.mx = ul.mx*wli[IVX] + ptl - bxsq;
    fl.my = ul.my*wli[IVX] - bxi*ul.by;
    fl.mz = ul.mz*wli[IVX] - bxi*ul.bz;
    fl.by = ul.by*wli[IVX] - bxi*wli[IVY];
    fl.bz = ul.bz*wli[IVX] - bxi*wli[IVZ];

    fr.d  = ur.mx;
    fr.mx = ur.mx*wri[IVX] + ptr - bxsq;
    fr.my = ur.my*wri[IVX] - bxi*ur.by;
    fr.mz = ur.mz*wri[IVX] - bxi*ur.bz;
    fr.by = ur.by*wri[IVX] - bxi*wri[IVY];
    fr.bz = ur.bz*wri[IVX] - bxi*wri[IVZ];

    //--- Step 4.  Compute hll averages and Alfven wave speed

    // inverse of difference between right and left signal speeds
    Real idspd = 1.0/(spd[4]-spd[0]);

    // rho component of U^{hll} from Mignone eqn. (15); uses F_L and F_R from eqn. (6)
    Real dhll = (spd[4]*ur.d - spd[0]*ul.d - fr.d + fl.d)*idspd;
    dhll = std::max(dhll, dfloor);
    Real sqrtdhll = std::sqrt(dhll);

    // rho and mx components of F^{hll} from Mignone eqn. (17)
    Real fdhll  = (spd[4]*fl.d  - spd[0]*fr.d  + spd[4]*spd[0]*(ur.d -ul.d ))*idspd;
    Real fmxhll = (spd[4]*fl.mx - spd[0]*fr.mx + spd[4]*spd[0]*(ur.mx-ul.mx))*idspd;

    // ustar from paragraph between eqns. (23) and (24)
    Real ustar = fdhll/dhll;

    // mx component of U^{hll} from Mignone eqn. (15); paragraph referenced
    // above states that mxhll should NOT be used to compute ustar
    Real mxhll = (spd[4]*ur.mx - spd[0]*ul.mx - fr.mx + fl.mx)*idspd;

    // S*_L and S*_R from Mignone eqn. (29)
    spd[1] = ustar - std::abs(bxi)/sqrtdhll;
    spd[3] = ustar + std::abs(bxi)/sqrtdhll;

    //--- Step 5. Compute intermediate states

    // Ul* - eqn. (20) of Mignone
    ulst.d  = dhll;
    ulst.mx = mxhll; // eqn. (24) of Mignone

    Real tmp = (spd[0]-spd[1])*(spd[0]-spd[3]);
    if (std::abs(spd[0]-spd[1]) < (SMALL_NUMBER)*cs) {
      // degenerate case described below eqn. (39)
      ulst.my = ul.my;
      ulst.mz = ul.mz;
      ulst.by = ul.by;
      ulst.bz = ul.bz;
    } else {
      Real mfact = bxi*(ustar-wli[IVX])/tmp;
      Real bfact = (ul.d*SQR(spd[0]-wli[IVX]) - bxsq)/(dhll*tmp);

      ulst.my = dhll*wli[IVY] - ul.by*mfact; // eqn. (30) of Mignone
      ulst.mz = dhll*wli[IVZ] - ul.bz*mfact; // eqn. (31) of Mignone
      ulst.by = ul.by*bfact; // eqn. (32) of Mignone
      ulst.bz = ul.bz*bfact; // eqn. (33) of Mignone
    }

    // Ur* - eqn. (20) of Mignone */
    urst.d  = dhll;
    urst.mx = mxhll; // eqn. (24) of Mignone

    tmp = (spd[4]-spd[1])*(spd[4]-spd[3]);
    if (std::abs(spd[4]-spd[3]) < (SMALL_NUMBER)*cs) {
      // degenerate case described below eqn. (39)
      urst.my = ur.my;
      urst.mz = ur.mz;
      urst.by = ur.by;
      urst.bz = ur.bz;
    } else {
      Real mfact = bxi*(ustar-wri[IVX])/tmp;
      Real bfact = (ur.d*SQR(spd[4]-wri[IVX]) - bxsq)/(dhll*tmp);

      urst.my = dhll*wri[IVY] - ur.by*mfact; // eqn. (30) of Mignone
      urst.mz = dhll*wri[IVZ] - ur.bz*mfact; // eqn. (31) of Mignone
      urst.by = ur.by*bfact; // eqn. (32) of Mignone
      urst.bz = ur.bz*bfact; // eqn. (33) of Mignone
    }

    // Uc*
    Real x = sqrtdhll*(bxi > 0.0 ? 1.0 : -1.0); // from below eqn. (37) of Mignone
    ucst.d  = dhll;  // eqn. (20) of Mignone
    ucst.mx = mxhll; // eqn. (24) of Mignone
    ucst.my = 0.5*(ulst.my + urst.my + (urst.by-ulst.by)*x); // eqn. (34) of Mignone
    ucst.mz = 0.5*(ulst.mz + urst.mz + (urst.bz-ulst.bz)*x); // eqn. (35) of Mignone
    ucst.by = 0.5*(ulst.by + urst.by + (urst.my-ulst.my)/x); // eqn. (36) of Mignone
    ucst.bz = 0.5*(ulst.bz + urst.bz + (urst.mz-ulst.mz)/x); // eqn. (37) of Mignone

    //--- Step 6.  Compute flux

    if (spd[0] >= 0.0) {
      // return Fl if flow is supersonic, eqn. (38a) of Mignone
      flxi[IDN] = fl.d;
      flxi[IVX] = fl.mx;
      flxi[IVY] = fl.my;
      flxi[IVZ] = fl.mz;
      flxi[IBY] = fl.by;
      flxi[IBZ] = fl.bz;
    } else if (spd[4] <= 0.0) {
      // return Fr if flow is supersonic, eqn. (38e) of Mignone
      flxi[IDN] = fr.d;
      flxi[IVX] = fr.mx;
      flxi[IVY] = fr.my;
      flxi[IVZ] = fr.mz;
      flxi[IBY] = fr.by;
      flxi[IBZ] = fr.bz;
    } else if (spd[1] >= 0.0) {
      // return (Fl+Sl*(Ulst-Ul)), eqn. (38b) of Mignone
      flxi[IDN] = fl.d  + spd[0]*(ulst.d  - ul.d);
      flxi[IVX] = fl.mx + spd[0]*(ulst.mx - ul.mx);
      flxi[IVY] = fl.my + spd[0]*(ulst.my - ul.my);
      flxi[IVZ] = fl.mz + spd[0]*(ulst.mz - ul.mz);
      flxi[IBY] = fl.by + spd[0]*(ulst.by - ul.by);
      flxi[IBZ] = fl.bz + spd[0]*(ulst.bz - ul.bz);
    } else if (spd[3] <= 0.0) {
      // return (Fr+Sr*(Urst-Ur)), eqn. (38d) of Mignone
      flxi[IDN] = fr.d  + spd[4]*(urst.d  - ur.d);
      flxi[IVX] = fr.mx + spd[4]*(urst.mx - ur.mx);
      flxi[IVY] = fr.my + spd[4]*(urst.my - ur.my);
      flxi[IVZ] = fr.mz + spd[4]*(urst.mz - ur.mz);
      flxi[IBY] = fr.by + spd[4]*(urst.by - ur.by);
      flxi[IBZ] = fr.bz + spd[4]*(urst.bz - ur.bz);
    } else {
      // return Fcst, eqn. (38c) of Mignone, using eqn. (24)
      flxi[IDN] = dhll*ustar;
      flxi[IVX] = fmxhll;
      flxi[IVY] = ucst.my*ustar - bxi*ucst.by;
      flxi[IVZ] = ucst.mz*ustar - bxi*ucst.bz;
      flxi[IBY] = ucst.by*ustar - bxi*ucst.my/ucst.d;
      flxi[IBZ] = ucst.bz*ustar - bxi*ucst.mz/ucst.d;
    }

    flx(IDN,k,j,i) = flxi[IDN];
    flx(ivx,k,j,i) = flxi[IVX];
    flx(ivy,k,j,i) = flxi[IVY];
    flx(ivz,k,j,i) = flxi[IVZ];
    ey(k,j,i) = -flxi[IBY];
    ez(k,j,i) =  flxi[IBZ];

    wct(k,j,i)=GetWeightForCT(flxi[IDN], wli[IDN], wri[IDN], dxw(i), dt);
  }

  return;
}
