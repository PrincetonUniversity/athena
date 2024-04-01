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
//! \file absorption_scattering.cpp
//  \brief implementation of absorption-scattering source terms
//======================================================================================


// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../eos/eos.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../utils/utils.hpp"
#include "../../radiation.hpp"

// this class header
#include "../rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::AbsorptionScattering()
//  \brief

// wmu_cm is the weight in the co-moving frame
// wmu_cm=wmu * 1/(1-vdotn/Crat)^2 / Lorz^2
// tran_coef is (1-vdotn/Crat)*Lorz
// rho is gas density
// tgas is gas temperature
// This function updates normal absorption plus scattering opacity together

Real RadIntegrator::AbsorptionScattering(
    AthenaArray<Real> &wmu_cm, AthenaArray<Real> &tran_coef, Real *sigma_a,
    Real *sigma_p, Real *sigma_pe, Real *sigma_s, Real dt, Real lorz,
    Real rho, Real &tgas, AthenaArray<Real> &implicit_coef, AthenaArray<Real> &ir_cm) {
  const Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->crat;
  Real redfactor=pmy_rad->reduced_c/pmy_rad->crat;
  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->peos->GetGamma();

  // Temporary array
  //AthenaArray<Real> &vncsigma = vncsigma_;
  AthenaArray<Real> &vncsigma2 = vncsigma2_;

  int badcell=0;


  Real coef[2];
  for (int ci=0; ci<2; ++ci)
    coef[ci] = 0.0;

  Real tgasnew = tgas;

  for (int ifr=0; ifr<nfreq; ++ifr) {
    Real suma1=0.0, suma2=0.0, suma3=0.0;
    Real jr_cm=0.0;

    Real dtcsigmar = ct * sigma_a[ifr]; //Rosseland mean absorption
    Real dtcsigmas = ct * sigma_s[ifr]; // scattering
    Real dtcsigmap = ct * sigma_p[ifr];// planck mean absorption
    Real dtcsigmae = ct * sigma_pe[ifr]; // Energy (planck) mean absorption

    Real rdtcsigmar = redfactor * dtcsigmar;
    Real rdtcsigmas = redfactor * dtcsigmas;
    Real rdtcsigmap = redfactor * dtcsigmap;
    Real rdtcsigmae = redfactor * dtcsigmae;


    Real *ircm = &(ir_cm(nang*ifr));
    Real *vn2 = &(vncsigma2(0));
    Real *tcoef = &(tran_coef(0));
    Real *wmu = &(wmu_cm(0));
    Real *imcoef = &(implicit_coef(0));
    for (int n=0; n<nang; n++) {
      Real vn = 1.0/(imcoef[n] + (rdtcsigmar + rdtcsigmas) * tcoef[n]);
      vn2[n] = tcoef[n] * vn;
      Real ir_weight = ircm[n] * wmu[n];

      suma1 += (wmu[n] * vn2[n]);
      suma2 += (ir_weight * vn);
    }
    suma3 = suma1 * (rdtcsigmas + rdtcsigmar - rdtcsigmae);
    suma1 *= (rdtcsigmap);

    // Now solve the equation
    // rho dT/gamma-1=-\gamma Prat c(sigma T^4 - sigma(a1 T^4 + a2)/(1-a3))


    // No need to do this if already in thermal equilibrium
    coef[1] = prat * (dtcsigmap - dtcsigmae * suma1/(1.0-suma3))
              * (gamma - 1.0)/rho;
    coef[0] = -tgas - dtcsigmae * prat * suma2 * (gamma - 1.0)/(rho*(1.0-suma3));

    if (std::abs(coef[1]) > TINY_NUMBER) {
      int flag = FouthPolyRoot(coef[1], coef[0], tgasnew);
      if (flag == -1 || tgasnew != tgasnew) {
        badcell = 1;
        tgasnew = tgas;
      }
    } else {
      tgasnew = -coef[0];
    }
    // even if tr=told, there can be change for intensity, making them isotropic
    if (!badcell) {
      Real emission = tgasnew * tgasnew * tgasnew * tgasnew;

      // get update jr_cm
      jr_cm = (suma1 * emission + suma2)/(1.0-suma3);

      // Update the co-moving frame specific intensity
      Real *irn = &(ir_cm(nang*ifr));
      Real *vn2 = &(vncsigma2(0));
      Real *imcoef = &(implicit_coef(nang*ifr));
      Real *tcoef = &(tran_coef(0));
      for (int n=0; n<nang; n++) {
        irn[n] +=  ((rdtcsigmas + rdtcsigmar - rdtcsigmae) * jr_cm + rdtcsigmap * emission
                    - ((imcoef[n]-1.0)/tcoef[n] + rdtcsigmas
                       + rdtcsigmar) * irn[n]) * vn2[n];
        irn[n] = std::max(irn[n],static_cast<Real>(TINY_NUMBER));
      }
    }
  }

  // Update gas temperature
  // Do not update gas temperature
  //  tgas = tgasnew;
  return tgasnew;
}
