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
/*! \file hllc.cpp
 *  \brief HLLC Riemann solver for hydrodynamics
 *====================================================================================*/

void Fluid::HLLC(const int il, const int iu,
       AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  Real sqrtdl,sqrtdr,isdlpdr,v1roe,v2roe,v3roe,hroe;
  Real cfl,cfr,bp,bm;
  Real al,ar; /* Min and Max wave speeds */
  Real am,cp; /* Contact wave speed and pressure */
  Real tl,tr,ml,mr,sl,sm,sr;
  Real Ul_E,Ur_E,Ul_Mx,Ur_Mx,vsq,asq,a,evp,evm;
  Real Fl[NVAR],Fr[NVAR];
  Real qa,qb;

//  Real Gamma = Fluid::GetGamma();
  Real Gamma = 1.6666667;

#pragma simd
  for (int i=il; i<=iu; ++i){

/* Compute Roe-averaged data from left- and right-states */

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
      

    sqrtdl = sqrt(d_l);
    sqrtdr = sqrt(d_r);
    isdlpdr = 1.0/(sqrtdl + sqrtdr);

    v1roe = (sqrtdl*vx_l + sqrtdr*vx_r)*isdlpdr;
    v2roe = (sqrtdl*vy_l + sqrtdr*vy_r)*isdlpdr;
    v3roe = (sqrtdl*vz_l + sqrtdr*vz_r)*isdlpdr;

/* Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */

    Ul_E = e_l + 0.5*d_l*
          (vx_l*vx_l + vy_l*vy_l + vz_l*vz_l);
    Ur_E = e_r + 0.5*d_r*
          (vx_r*vx_r + vy_r*vy_r + vz_r*vz_r);
    Ul_Mx = d_l*vx_l;
    Ur_Mx = d_r*vx_r;
    
    hroe = ((Ul_E + e_l)/sqrtdl + (Ur_E + e_r)/sqrtdr)*isdlpdr;
    vsq = v1roe*v1roe + v2roe*v2roe + v3roe*v3roe;
    asq = (Gamma-1.0)*std::max((hroe - 0.5*vsq), TINY_NUMBER);
    a = sqrt(asq);
    evp = v1roe + a;
    evm = v1roe - a;

/* Compute the max and min wave speeds */

    cfl = sqrt((Gamma*e_l/d_l));
    cfr = sqrt((Gamma*e_r/d_r));

    ar = std::max(evp,(vx_r + cfr));
    al = std::min(evm,(vx_l - cfl));

    bp = ar > 0.0 ? ar : 0.0;
    bm = al < 0.0 ? al : 0.0;

/* Compute the contact wave speed and Pressure */

    tl = e_l + (vx_l - al)*Ul_Mx;
    tr = e_r + (vx_r - ar)*Ur_Mx;

    ml =   Ul_Mx - d_l*al;
    mr = -(Ur_Mx - d_r*ar);

/* Determine the contact wave speed... */
    am = (tl - tr)/(ml + mr);
/* ...and the pressure at the contact surface */
    cp = (ml*tr + mr*tl)/(ml + mr);
    cp = cp > 0.0 ? cp : 0.0;

/* Compute L/R fluxes along the line bm, bp */

    Fl[IDN]  = Ul_Mx - bm*d_l;
    Fr[IDN]  = Ur_Mx - bp*d_r;

    Fl[IVX] = Ul_Mx*(vx_l - bm);
    Fr[IVX] = Ur_Mx*(vx_r - bp);

    Fl[IVY] = d_l*vy_l*(vx_l - bm);
    Fr[IVY] = d_r*vy_r*(vx_r - bp);

    Fl[IVZ] = d_l*vz_l*(vx_l - bm);
    Fr[IVZ] = d_r*vz_r*(vx_r - bp);

    Fl[IVX] += e_l;
    Fr[IVX] += e_r;

    Fl[IEN] = Ul_E*(vx_l - bm) + e_l*vx_l;
    Fr[IEN] = Ur_E*(vx_r - bp) + e_r*vx_r;

/* Compute flux weights or scales */

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

/* Compute the HLLC flux at interface */

    d_Flx = sl*Fl[IDN] + sr*Fr[IDN];
    vx_Flx = sl*Fl[IVX] + sr*Fr[IVX];
    vy_Flx = sl*Fl[IVY] + sr*Fr[IVY];
    vz_Flx = sl*Fl[IVZ] + sr*Fr[IVZ];
    e_Flx = sl*Fl[IEN] + sr*Fr[IEN];

/* Add the weighted contribution of the flux along the contact */

    vx_Flx += sm*cp;
    e_Flx += sm*cp*am;
  }
  return;
}
