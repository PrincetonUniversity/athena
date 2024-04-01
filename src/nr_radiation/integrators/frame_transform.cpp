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
//! \file frame_transform.cpp
//  \brief  Perform frame transformation
//======================================================================================

// C headers

// C++ headers
#include <sstream>  // msg
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
#include "../radiation.hpp"

// this class header
#include "./rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::LabToCom(const Real vx, const Real vy, const Real vz,
//                          AthenaArray<Real> &ir, AthenaArray<Real> &ir_cm)
//  \brief Transform specific intensity from lab frame to co-moving frame
// with flow velocity vx, vy, vz


// this is only used for one frequency group
void RadIntegrator::LabToCom(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          Real *ir_lab, AthenaArray<Real> &ir_cm) {
  //Real& prat = pmy_rad->prat;
  Real invcrat = 1.0/pmy_rad->crat;
  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;

  // square of Lorentz factor
  Real lorzsq = 1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat);

  for (int ifr=0; ifr<nfreq; ++ifr) {
    for (int n=0; n<nang; n++) {
      Real vnc = vx * mux[n] + vy * muy[n] + vz * muz[n];
      vnc = 1.0 - vnc * invcrat;
      Real coef = vnc * vnc * vnc * vnc * lorzsq * lorzsq;
      ir_cm(ifr*nang+n) = ir_lab[ifr*nang+n] * coef;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ComToLab(const Real vx, const Real vy, const Real vz,
//                          AthenaArray<Real> &ir, AthenaArray<Real> &ir_cm)
//  \brief Transform specific intensity from co-moving frame to lab frame
// with flow velocity vx, vy, vz

void RadIntegrator::ComToLab(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          AthenaArray<Real> &ir_cm, Real *ir_lab) {
  Real invcrat = 1.0/pmy_rad->crat;
  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;

  // square of Lorentz factor
  Real lorzsq = 1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat);

  for (int ifr=0; ifr<nfreq; ++ifr) {
    for (int n=0; n<nang; n++) {
      Real vnc = vx * mux[n] + vy * muy[n] + vz * muz[n];

      vnc = 1.0 - vnc * invcrat;
      Real coef = vnc * vnc * vnc * vnc * lorzsq * lorzsq;

      ir_lab[ifr*nang+n] = ir_cm(ifr*nang+n) / coef;
    }
  }


  return;
}

// transform co-moving frame intensity to the lab frame for a given spectrum in
// co-moving frame
void RadIntegrator::ComToLabMultiGroup(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          AthenaArray<Real> &ir_cm, Real *ir_lab) {
  Real invcrat = 1.0/pmy_rad->crat;
  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;
  AthenaArray<Real> split_ratio;
  AthenaArray<int> map_start, map_end;

  // square of Lorentz factor
  Real lorz = std::sqrt(1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat));
  // first, get the lorentz transformation factor
  // now calculate the actual transformation factor

  for (int n=0; n<nang; ++n) {
    Real vnc = vx * mux[n] + vy * muy[n] + vz * muz[n];
    vnc = 1.0 - vnc * invcrat;
    Real cm_nu = vnc * lorz;

    for (int ifr=0; ifr<nfreq; ++ifr)
      ir_ori_(ifr) = ir_cm(ifr*nang+n);

    split_ratio.InitWithShallowSlice(split_ratio_,3,n,1);
    map_start.InitWithShallowSlice(map_bin_start_,2,n,1);
    map_end.InitWithShallowSlice(map_bin_end_,2,n,1);

    MapCmToLabFrequency(cm_nu,split_ratio,map_start,map_end,ir_ori_,ir_done_);

    for (int ifr=0; ifr<nfreq; ++ifr) {
      ir_lab[ifr*nang+n] = ir_done_(ifr)/(cm_nu*cm_nu*cm_nu*cm_nu);
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ComAngle(const Real vx, const Real vy, const Real vz,
//                  Real mux, Real muy, Real muz, Real mux0, Real muy0, Real muz0)
//  \brief Transform angles from lab frame to co-moving frame


void RadIntegrator::ComAngle(const Real vx, const Real vy, const Real vz,
                  Real mux, Real muy, Real muz, Real *mux0, Real *muy0, Real *muz0) {
  Real invcrat = 1.0/pmy_rad->crat;
  // square of Lorentz factor
  Real lorz = 1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat);
  lorz = std::sqrt(lorz);
  Real vdotn = vx * mux + vy * muy + vz * muz;

  Real angcoef = lorz * invcrat * (1.0 - lorz * vdotn * invcrat/(1.0 + lorz));
  Real incoef = 1.0/(lorz*(1.0-vdotn * invcrat));

  (*mux0) = (mux - vx * angcoef) * incoef;
  (*muy0) = (muy - vy * angcoef) * incoef;
  (*muz0) = (muz - vz * angcoef) * incoef;
  return;
}
