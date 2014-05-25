//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"

//======================================================================================
/*! \file hlle.cpp
 *  \brief HLLE Riemann solver for hydrodynamics
 *
 *  Computes 1D fluxes using the Harten-Lax-van Leer (HLL) Riemann solver.  This flux is
 *  very diffusive, especially for contacts, and so it is not recommended for use in
 *  applications.  However, as shown by Einfeldt et al.(1991), it is positively
 *  conservative (cannot return negative densities or pressure), so it is a useful
 *  option when other approximate solvers fail and/or when extra dissipation is needed.
 *
 * REFERENCES:
 * - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
 *   Springer-Verlag, Berlin, (1999) chpt. 10.
 * - Einfeldt et al., "On Godunov-type methods near low densities", JCP, 92, 273 (1991)
 * - A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and Godunov-type
 *   schemes for hyperbolic conservation laws", SIAM Review 25, 35-61 (1983).
 *====================================================================================*/

/*
void RiemannSolver::HLLE(const int il, const int iu,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  Real cfl,cfr,bp,bm,tmp;
  Real evp, evm,,al,ar;
  Real Fl[NVAR],Fr[NVAR];

  Real Gamma = pmy_fluid_->GetGamma();

#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& d_l=wl(IDN,i);
    Real& vx_l=wl(IVX,i);
    Real& vy_l=wl(IVY,i);
    Real& vz_l=wl(IVZ,i);
    Real& e_l=wl(IEN,i);

    Real& d_r=wr(IDN,i);
    Real& vx_r=wr(IVX,i);
    Real& vy_r=wr(IVY,i);
    Real& vz_r=wr(IVZ,i);
    Real& e_r=wr(IEN,i);

    Real& d_Flx=flx(IDN,i);
    Real& vx_Flx=flx(IVX,i);
    Real& vy_Flx=flx(IVY,i);
    Real& vz_Flx=flx(IVZ,i);
    Real& e_Flx=flx(IEN,i);

// Compute Roe-averaged velocities

    Real sqrtdl = sqrt(d_l);
    Real sqrtdr = sqrt(d_r);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    Real v1roe = (sqrtdl*vx_l + sqrtdr*vx_r)*isdlpdr;
    Real v2roe = (sqrtdl*vy_l + sqrtdr*vy_r)*isdlpdr;
    Real v3roe = (sqrtdl*vz_l + sqrtdr*vz_r)*isdlpdr;

// Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
// rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl

    Real Ul_E = e_l + 0.5*d_l* (vx_l*vx_l + vy_l*vy_l + vz_l*vz_l);
    Real Ur_E = e_r + 0.5*d_r* (vx_r*vx_r + vy_r*vy_r + vz_r*vz_r);
    Real Ul_Mx = d_l*vx_l;
    Real Ur_Mx = d_r*vx_r;
    Real hroe = ((Ul_E + e_l)/sqrtdl + (Ur_E + e_r)/sqrtdr)*isdlpdr;

// Compute Roe-averaged wave speeds

    vsq = v1roe*v1roe + v2roe*v2roe + v3roe*v3roe;
    asq = (Gamma-1.0)*std::max((hroe - 0.5*vsq), TINY_NUMBER);
    a = sqrt(asq);
    evp = v1roe + a;
    evm = v1roe - a;

// Compute the max/min wave speeds based on L/R values and Roe averages

    cfl = sqrt((Gamma*e_l/d_l));
    cfr = sqrt((Gamma*e_r/d_r));

    ar = std::max(evp,(vx_r + cfr));
    al = std::min(evm,(vx_l - cfl));

    bp = ar > 0.0 ? ar : 0.0;
    bm = al < 0.0 ? al : 0.0;

// Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}

    Fl[IDN] = Ul.Mx - bm*Ul.d;
    Fr[IDN] = Ur.Mx - bp*Ur.d;

    Fl[IVX] = Ul.Mx*(Wl.Vx - bm) + Wl.P;
    Fr[IVX] = Ur.Mx*(Wr.Vx - bp) + Wr.P;

    Fl[IVY] = Ul.My*(Wl.Vx - bm);
    Fr[IVY] = Ur.My*(Wr.Vx - bp);

    Fl[IVZ] = Ul.Mz*(Wl.Vx - bm);
    Fr[IVZ] = Ur.Mz*(Wr.Vx - bp);

    Fl[IEN] = Ul.E*(Wl.Vx - bm) + Wl.P*Wl.Vx;
    Fr[IEN] = Ur.E*(Wr.Vx - bp) + Wr.P*Wr.Vx;

// Compute the HLLE flux at interface.

    tmp = 0.5*(bp + bm)/(bp - bm);

    d_Flx  = 0.5*(Fl[IDN]+Fr[IDN]) + (Fl[IDN]-Fr[IDN])*tmp;
    vx_Flx = 0.5*(Fl[IVX]+Fr[IVX]) + (Fl[IVX]-Fr[IVX])*tmp;
    vy_Flx = 0.5*(Fl[IVY]+Fr[IVY]) + (Fl[IVY]-Fr[IVY])*tmp;
    vz_Flx = 0.5*(Fl[IVZ]+Fr[IVZ]) + (Fl[IVZ]-Fr[IVZ])*tmp;
    e_Flx  = 0.5*(Fl[IEN]+Fr[IEN]) + (Fl[IEN]-Fr[IEN])*tmp;
  }

  return;
}
*/
