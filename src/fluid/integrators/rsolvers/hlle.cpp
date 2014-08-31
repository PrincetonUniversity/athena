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
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "../integrators/integrators.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../athena.hpp"         // enums, macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../fluid.hpp"          // Fluid

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

void FluidIntegrator::RiemannSolver(const int k, const int j,
  const int il, const int iu, const int ivx, const int ivy, const int ivz,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  Real cfl,cfr,bp,bm;
  Real evp,evm,al,ar;
  Real fl[NVAR],fr[NVAR];

  Real gamma = pmy_fluid->GetGamma();

#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& d_l=wl(IDN,i);
    Real& vx_l=wl(ivx,i);
    Real& vy_l=wl(ivy,i);
    Real& vz_l=wl(ivz,i);
    Real& p_l=wl(IEN,i);

    Real& d_r=wr(IDN,i);
    Real& vx_r=wr(ivx,i);
    Real& vy_r=wr(ivy,i);
    Real& vz_r=wr(ivz,i);
    Real& p_r=wr(IEN,i);

// Compute Roe-averaged velocities

    Real sqrtdl = sqrt(d_l);
    Real sqrtdr = sqrt(d_r);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    Real v1roe = (sqrtdl*vx_l + sqrtdr*vx_r)*isdlpdr;
    Real v2roe = (sqrtdl*vy_l + sqrtdr*vy_r)*isdlpdr;
    Real v3roe = (sqrtdl*vz_l + sqrtdr*vz_r)*isdlpdr;

// Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
// rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl

    Real ul_e = p_l/(gamma - 1.0) + 0.5*d_l* (vx_l*vx_l + vy_l*vy_l + vz_l*vz_l);
    Real ur_e = p_r/(gamma - 1.0) + 0.5*d_r* (vx_r*vx_r + vy_r*vy_r + vz_r*vz_r);
    Real ul_mx = d_l*vx_l;
    Real ur_mx = d_r*vx_r;
    Real hroe = ((ul_e + p_l)/sqrtdl + (ur_e + p_r)/sqrtdr)*isdlpdr;

// Compute Roe-averaged wave speeds

    Real vsq = v1roe*v1roe + v2roe*v2roe + v3roe*v3roe;
    Real asq = (gamma-1.0)*std::max((hroe - 0.5*vsq), TINY_NUMBER);
    Real a = sqrt(asq);
    evp = v1roe + a;
    evm = v1roe - a;

// Compute the max/min wave speeds based on L/R values and Roe averages

    cfl = sqrt((gamma*p_l/d_l));
    cfr = sqrt((gamma*p_r/d_r));

    ar = std::max(evp,(vx_r + cfr));
    al = std::min(evm,(vx_l - cfl));

    bp = ar > 0.0 ? ar : 0.0;
    bm = al < 0.0 ? al : 0.0;

// Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}

    fl[IDN] = ul_mx - bm*d_l;
    fr[IDN] = ur_mx - bp*d_r;

    fl[ivx] = ul_mx*(vx_l - bm) + p_l;
    fr[ivx] = ur_mx*(vx_r - bp) + p_r;

    fl[ivy] = d_l*vy_l*(vx_l - bm);
    fr[ivy] = d_r*vy_r*(vx_r - bp);

    fl[ivz] = d_l*vz_l*(vx_l - bm);
    fr[ivz] = d_r*vz_r*(vx_r - bp);

    fl[IEN] = ul_e*(vx_l - bm) + p_l*vx_l;
    fr[IEN] = ur_e*(vx_r - bp) + p_r*vx_r;

// Compute the HLLE flux at interface.

    Real& d_flx=flx(IDN,i);
    Real& vx_flx=flx(ivx,i);
    Real& vy_flx=flx(ivy,i);
    Real& vz_flx=flx(ivz,i);
    Real& e_flx=flx(IEN,i);

    Real tmp = 0.5*(bp + bm)/(bp - bm);

    d_flx  = 0.5*(fl[IDN]+fr[IDN]) + (fl[IDN]-fr[IDN])*tmp;
    vx_flx = 0.5*(fl[ivx]+fr[ivx]) + (fl[ivx]-fr[ivx])*tmp;
    vy_flx = 0.5*(fl[ivy]+fr[ivy]) + (fl[ivy]-fr[ivy])*tmp;
    vz_flx = 0.5*(fl[ivz]+fr[ivz]) + (fl[ivz]-fr[ivz])*tmp;
    e_flx  = 0.5*(fl[IEN]+fr[IEN]) + (fl[IEN]-fr[IEN])*tmp;
  }

  return;
}
