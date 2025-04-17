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
//! \file multigroup_abs_sca.cpp
//  \brief  Add multi-group absorption scattering
//======================================================================================


// C headers

// C++ headers

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
//! \fn RadIntegrator::MultiGroupAbsScat()
//  \brief

// wmu_cm is the weight in the co-moving frame
// wmu_cm=wmu * 1/(1-vdotn/Crat)^2 / Lorz^2
// tran_coef is (1-vdotn/Crat)*Lorz
// rho is gas density
// tgas is gas temperature
// This function updates normal absorption plus scattering opacity together

// ir_cm here is the frequency shifted array
// emission in each frequency bin is assumed to b T^4 emission_spec[ifr]

Real RadIntegrator::MultiGroupAbsScat(
    AthenaArray<Real> &wmu_cm, AthenaArray<Real> &tran_coef, Real *sigma_a, Real *sigma_p,
    Real *sigma_pe, Real *sigma_s, Real dt, Real lorz, Real rho, Real &tgas,
    AthenaArray<Real> &implicit_coef, AthenaArray<Real> &ir_cm) {
  const Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->crat;
  Real redfactor=pmy_rad->reduced_c/pmy_rad->crat;
  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->peos->GetGamma();

  // Temporary array
  AthenaArray<Real> &vncsigma2 = vncsigma2_;

  bool badcell=false;
  Real coef[2];

  Real *suma1 = &(sum_nu1_(0));
  Real *suma2 = &(sum_nu2_(0));
  Real *suma3 = &(sum_nu3_(0));

  Real tgasnew = tgas;
  Real tgas_last = tgas;

  int count_iteration = 0;

  // first, calculate the coefficients that are independnet of emission_spec
  coef[0] = 0.0;
  for (int ifr=0; ifr<nfreq; ++ifr) {
    suma1[ifr]=0.0;
    suma2[ifr]=0.0;
    suma3[ifr]=0.0;

    Real dtcsigmar = ct * sigma_a[ifr];
    Real dtcsigmas = ct * sigma_s[ifr];
    Real dtcsigmae = ct * sigma_pe[ifr];
    Real dtcsigmap = ct * sigma_p[ifr];

    Real rdtcsigmar = redfactor * dtcsigmar;
    Real rdtcsigmas = redfactor * dtcsigmas;
    Real rdtcsigmae = redfactor * dtcsigmae;
    Real rdtcsigmap = redfactor * dtcsigmap;

    Real *ircm = &(ir_cm(ifr*nang));
    Real *vn2 = &(vncsigma2(ifr*nang));
    Real *tcoef = &(tran_coef(0));
    Real *wmu = &(wmu_cm(0));
    Real *imcoef = &(implicit_coef(ifr*nang));
    for (int n=0; n<nang; n++) {
      Real vn = 1.0/(imcoef[n] + (rdtcsigmar + rdtcsigmas) * tcoef[n]);
      vn2[n] = tcoef[n] * vn;
      Real ir_weight = ircm[n] * wmu[n];
      suma1[ifr] += (wmu[n] * vn2[n]);
      suma2[ifr] += (ir_weight * vn);
    }
    suma3[ifr] = suma1[ifr] * (rdtcsigmas + rdtcsigmar - rdtcsigmae);
    suma1[ifr] *= (rdtcsigmap);

    coef[0] += - dtcsigmae * prat * suma2[ifr]
               * (gamma - 1.0)/(rho*(1.0-suma3[ifr]));
  }
  coef[0] += -tgas;

  Real tgas_diff = 1.0;

  while ((count_iteration <= iteration_tgas_) && (tgas_diff > tgas_error_)) {
    coef[1] = 0.0;

    // update emission spectrum
    pmy_rad->UserEmissionSpec(pmy_rad,tgasnew);

    for (int ifr=0; ifr<nfreq; ++ifr) {
      // Real dtcsigmar = ct * sigma_a[ifr];
      Real dtcsigmae = ct * sigma_pe[ifr];
      //      Real dtcsigmas = ct * sigma_s[ifr];
      // Real rdtcsigmar = redfactor * dtcsigmar;
      // Real rdtcsigmae = redfactor * dtcsigmae;
      //      Real rdtcsigmas = redfactor * dtcsigmas;
      Real dtcsigmap = ct * sigma_p[ifr];
      // Real rdtcsigmap = redfactor * dtcsigmap;

      // No need to do this if already in thermal equilibrium
      coef[1] += prat * (dtcsigmap - dtcsigmae * suma1[ifr]/(1.0-suma3[ifr]))
                 * pmy_rad->emission_spec(ifr) * (gamma - 1.0)/rho;
    }

    // The polynomial is
    // coef[1] * x^4 + x + coef[0] == 0

    if (std::abs(coef[1]) > TINY_NUMBER) {
      int flag = FouthPolyRoot(coef[1], coef[0], tgasnew);
      if (flag == -1 || tgasnew != tgasnew) {
        badcell = true;
        tgasnew = tgas;
      }
    } else {
      tgasnew = -coef[0];
    }

    if (badcell) break;
    tgas_diff = std::abs(tgas_last - tgasnew)/tgas_last;
    count_iteration++;
    tgas_last = tgasnew;
  }
  // even if tr=told, there can be change for intensity, making them isotropic
  if (!badcell) {
    Real emission = tgasnew * tgasnew * tgasnew * tgasnew;
    for (int ifr=0; ifr<nfreq; ++ifr) {
      Real dtcsigmar = ct * sigma_a[ifr];
      Real dtcsigmae = ct * sigma_pe[ifr];
      Real dtcsigmas = ct * sigma_s[ifr];
      Real dtcsigmap = ct * sigma_p[ifr];

      Real rdtcsigmar = redfactor * dtcsigmar;
      Real rdtcsigmae = redfactor * dtcsigmae;
      Real rdtcsigmas = redfactor * dtcsigmas;
      Real rdtcsigmap = redfactor * dtcsigmap;

      Real emi_nu = emission * pmy_rad->emission_spec(ifr);

      // get update jr_cm
      Real jr_cm = (suma1[ifr] * emi_nu + suma2[ifr])/(1.0-suma3[ifr]);

      // Update the co-moving frame specific intensity
      Real *irn = &(ir_cm(ifr*nang));
      Real *vn2 = &(vncsigma2(ifr*nang));
      Real *imcoef = &(implicit_coef(nang*ifr));
      Real *tcoef = &(tran_coef(0));
      for (int n=0; n<nang; n++) {
        irn[n] +=  ((rdtcsigmas + rdtcsigmar - rdtcsigmae) * jr_cm + rdtcsigmap * emi_nu
                    - ((imcoef[n]-1.0)/tcoef[n] + rdtcsigmas
                       + rdtcsigmar) * irn[n]) * vn2[n];
        irn[n] = std::max(irn[n], static_cast<Real>(TINY_NUMBER));
      }
    }
  }
  // Update gas temperature
  // Do not update gas temperature
  //  tgas = tgasnew;
  return tgasnew;
}
