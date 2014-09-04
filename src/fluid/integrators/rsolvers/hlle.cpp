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
#include "../integrators.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../../athena.hpp"         // enums, macros, Real
#include "../../../athena_arrays.hpp"  // AthenaArray
#include "../../fluid.hpp"             // Fluid
#include "../../eos/eos.hpp"           // GetGamma

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
  AthenaArray<Real> *pwl, AthenaArray<Real> *pwr, AthenaArray<Real> *pflx)
{
  Real fl[NVAR],fr[NVAR],wli[NVAR],wri[NVAR],flxi[NVAR];

#pragma simd
  for (int i=il; i<=iu; ++i){
    wli[IDN]=(*pwl)(IDN,i);
    wli[IVX]=(*pwl)(ivx,i);
    wli[IVY]=(*pwl)(ivy,i);
    wli[IVZ]=(*pwl)(ivz,i);
    if (NON_BAROTROPIC_EOS) wli[IEN]=(*pwl)(IEN,i);

    wri[IDN]=(*pwr)(IDN,i);
    wri[IVX]=(*pwr)(ivx,i);
    wri[IVY]=(*pwr)(ivy,i);
    wri[IVZ]=(*pwr)(ivz,i);
    if (NON_BAROTROPIC_EOS) wri[IEN]=(*pwr)(IEN,i);

// Compute the max/min wave speeds based on L/R values

    Real al = wli[IVX] - pmy_fluid->pf_eos->SoundSpeed(wli);
    Real ar = wri[IVX] + pmy_fluid->pf_eos->SoundSpeed(wri);

    Real bp = ar > 0.0 ? ar : 0.0;
    Real bm = al < 0.0 ? al : 0.0;

// Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}

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
      fl[IEN] = wli[IEN]/(pmy_fluid->pf_eos->GetGamma() - 1.0) + 0.5*wli[IDN]*
        (wli[IVX]*wli[IVX] + wli[IVY]*wli[IVY] + wli[IVZ]*wli[IVZ]);
      fr[IEN] = wri[IEN]/(pmy_fluid->pf_eos->GetGamma() - 1.0) + 0.5*wri[IDN]*
        (wri[IVX]*wri[IVX] + wri[IVY]*wri[IVY] + wri[IVZ]*wri[IVZ]);
      fl[IEN] *= (wli[IVX] - bm);
      fr[IEN] *= (wri[IVX] - bp);
      fl[IEN] += wli[IEN]*wli[IVX];
      fr[IEN] += wri[IEN]*wri[IVX];
    } else {
      Real iso_cs = pmy_fluid->pf_eos->SoundSpeed(wli);
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

// Compute the HLLE flux at interface.

    Real tmp = 0.5*(bp + bm)/(bp - bm);

    flxi[IDN] = 0.5*(fl[IDN]+fr[IDN]) + (fl[IDN]-fr[IDN])*tmp;
    flxi[IVX] = 0.5*(fl[IVX]+fr[IVX]) + (fl[IVX]-fr[IVX])*tmp;
    flxi[IVY] = 0.5*(fl[IVY]+fr[IVY]) + (fl[IVY]-fr[IVY])*tmp;
    flxi[IVZ] = 0.5*(fl[IVZ]+fr[IVZ]) + (fl[IVZ]-fr[IVZ])*tmp;
    if (NON_BAROTROPIC_EOS) flxi[IEN] = 0.5*(fl[IEN]+fr[IEN]) + (fl[IEN]-fr[IEN])*tmp;

    (*pflx)(IDN,i) = flxi[IDN];
    (*pflx)(ivx,i) = flxi[IVX];
    (*pflx)(ivy,i) = flxi[IVY];
    (*pflx)(ivz,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) (*pflx)(IEN,i) = flxi[IEN];
  }

  return;
}
