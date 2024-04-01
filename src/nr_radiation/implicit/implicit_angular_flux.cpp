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
//! \file implicit_angular_flux.cpp
//  \brief implementation of angular flux implicitly
//======================================================================================


// C headers

// C++ headers
#include <sstream>    // stringstream

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../reconstruct/reconstruction.hpp"
#include "../radiation.hpp"

// class header
#include "../integrators/rad_integrators.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// calculate advective flux for the implicit scheme
// only work for 1D spherical polar
// The equation to solve is
// \partial /\partial t - c/r \partial(\sin^2 zeta I)/\sin \zeta\partial \zeta
// the output is stored in angflx

// calculate the coefficient for angular flux
void RadIntegrator::ImplicitAngularFluxesCoef(const Real wght) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;
  Coordinates *pco=pmb->pcoord;
  std::stringstream msg;

  const int& nzeta = prad->nzeta;
  //int &npsi = prad->npsi;
  AthenaArray<Real> &area_zeta = zeta_area_, &ang_vol = ang_vol_;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
        if (prad->npsi == 0) {
        // first, related all angle to 2*nzeta
          pco->GetGeometryZeta(prad,k,j,i,g_zeta_);
        // The relation is I_m = coef *I_m+1 + const

        // now, all the other angles
        //  Ir_new - ir_old + coef0 * Ir_l + coef1 ir_ new = 0.0;
          for (int n=0; n<2*nzeta; ++n) {
            Real coef0 = wght * g_zeta_(n) * prad->reduced_c *
                       area_zeta(n)/ang_vol(n);
            Real coef1 = -wght * g_zeta_(n+1) * prad->reduced_c *
                       area_zeta(n+1)/ang_vol(n);
          // the equation is
          // ir_new - ir_old + coef0 * ir_new + coef1 * ir_new1 = 0
            if (n == 2*nzeta-1) {
              imp_ang_coef_r_(k,j,i,n) = 0.0;
              imp_ang_coef_(k,j,i,n) = coef0 + coef1;
            } else {
              imp_ang_coef_r_(k,j,i,n) = coef1;
              imp_ang_coef_(k,j,i,n) = coef0;
            }
          }
        } else {
          // first, starting from the zeta angle 2*nzeta-1
          // zeta area is only 0 at nzeta=0, not zero at nzeta
          pco->GetGeometryZeta(prad,k,j,i,g_zeta_);

          // now go from 2*nzeta-2 to 0
          for (int n=0; n<2*nzeta; ++n) {
            Real zeta_coef0 =  wght * g_zeta_(n) * prad->reduced_c;
            Real zeta_coef1 = -wght * g_zeta_(n+1) * prad->reduced_c;
            ImplicitPsiFluxCoef(k,j,i, n, wght, zeta_coef1, zeta_coef0);
          }
        }
      }
}


void RadIntegrator::ImplicitPsiFluxCoef(int k, int j, int i, int n_zeta, Real wght,
                                                   Real zeta_coef1, Real zeta_coef) {
  NRRadiation *prad=pmy_rad;
  Coordinates *pco=prad->pmy_block->pcoord;
  AthenaArray<Real> &area_psi = psi_area_, &ang_vol = ang_vol_, &zeta_area = zeta_area_;
  const int& npsi = prad->npsi;
  // the equation to solve
  //(1+zeta_coef0+zeta_coef1) I + Div F_psi = 0
  pco->GetGeometryPsi(prad,k,j,i,n_zeta,g_psi_);
  // g_psi_ =sin zeta * cot \theta sin\psi/r
  // g_psi_(0) is always 0

  for (int m=0; m<2*npsi; ++m) { // all take the left state
    int ang_num = n_zeta*2*npsi+m;

    Real coef0 = wght * g_psi_(m) * prad->reduced_c *
            area_psi(n_zeta,m)/ang_vol(ang_num);
    Real coef1 = -wght * g_psi_(m+1) * prad->reduced_c *
            area_psi(n_zeta,m+1)/ang_vol(ang_num);
    Real z_coef1 = zeta_coef1 * zeta_area(m,n_zeta+1)/ang_vol(ang_num);
    Real z_coef = zeta_coef * zeta_area(m,n_zeta)/ang_vol(ang_num);
    Real coef_c = 0.0;
    if ((g_psi_(m) < 0) || (g_psi_(m+1) < 0)) {
      imp_ang_psi_l_(k,j,i,ang_num) = coef0;
      imp_ang_psi_r_(k,j,i,ang_num) = 0.0;
      coef_c = coef1;
    } else if ((g_psi_(m) > 0) || (g_psi_(m+1) > 0)) {
      imp_ang_psi_l_(k,j,i,ang_num) = 0.0;
      imp_ang_psi_r_(k,j,i,ang_num) = coef1;
      coef_c = coef0;
    }
    if (n_zeta == 2*prad->nzeta-1) {
      imp_ang_coef_r_(k,j,i,ang_num) = 0.0;
      imp_ang_coef_(k,j,i,ang_num) = coef_c + z_coef + z_coef1;
    } else {
      imp_ang_coef_r_(k,j,i,ang_num) = z_coef1;
      imp_ang_coef_(k,j,i,ang_num) = coef_c + z_coef;
    }
  }
}


void RadIntegrator::ImplicitAngularFluxes(const int k, const int j, const int i,
                                          AthenaArray<Real> &ir_ini) {
  NRRadiation *prad=pmy_rad;
  std::stringstream msg;

  const int& nzeta = prad->nzeta;
  // const int& npsi = prad->npsi;
  const int& nang = prad->nang;
  const int& nfreq = prad->nfreq;

  // AthenaArray<Real> &area_zeta = zeta_area_, &ang_vol = ang_vol_;

  for (int ifr=0; ifr<nfreq; ++ifr) {
    if (prad->npsi == 0) {
      // the angle 2*nzeta-1 does not contribute to ang_flx_
      Real *p_angflx = &(ang_flx_(k,j,i,ifr*nang));
      Real *coef_r = &(imp_ang_coef_r_(k,j,i,0));
      Real *p_ir = &(ir_ini(k,j,i,ifr*nang));
      for (int n=0; n<2*nzeta-1; ++n) {
        p_angflx[n] = -coef_r[n] * p_ir[n+1];
      }
    } else {
      // now go from 2*nzeta-2 to 0
      for (int n=0; n<2*nzeta; ++n) {
        ImplicitPsiFlux(k,j,i, ifr, n, ir_ini);
      }
    }
  }
}


void RadIntegrator::ImplicitPsiFlux(const int k, const int j, const int i,
                      int ifr, int n_zeta, AthenaArray<Real> &ir_ini) {
  NRRadiation *prad=pmy_rad;
  const int& npsi = prad->npsi;
  // m=0
  int ang_num = n_zeta*2*npsi;
  const int& nang = prad->nang;

  Real *psi_l = &(imp_ang_psi_l_(k,j,i,ang_num));
  Real *psi_r = &(imp_ang_psi_r_(k,j,i,ang_num));
  Real *p_ir = &(ir_ini(k,j,i,ifr*nang+ang_num));
  Real *p_angflx = &(ang_flx_(k,j,i,ifr*nang+ang_num));

  if (n_zeta == 2*prad->nzeta-1) {
    // m=0
    p_angflx[0] = -(psi_l[0] * p_ir[2*npsi-1] + psi_r[0] * p_ir[1]);
    for (int m=1; m<2*npsi-1; ++m) { // all take the left state
      p_angflx[m] = -(psi_l[m] * p_ir[m-1] + psi_r[m] * p_ir[m+1]);
    }
    //m=2*npsi-1
    p_angflx[2*npsi-1] = -(psi_l[2*npsi-1] * p_ir[2*npsi-2] + psi_r[2*npsi-1] * p_ir[0]);
  } else {
    Real *p_ir_zetar = &(ir_ini(k,j,i,ifr*nang+(n_zeta+1)*2*npsi));
    Real *zeta_r = &(imp_ang_coef_r_(k,j,i,ang_num));

    p_angflx[0] = -(psi_l[0] * p_ir[2*npsi-1] + psi_r[0] * p_ir[1])
                  -zeta_r[0] * p_ir_zetar[0];
    for (int m=1; m<2*npsi-1; ++m) { // all take the left state
      p_angflx[m] = -(psi_l[m] * p_ir[m-1] + psi_r[m] * p_ir[m+1])
                    -zeta_r[m] * p_ir_zetar[m];
    }
    //m=2*npsi-1
    p_angflx[2*npsi-1] = -(psi_l[2*npsi-1] * p_ir[2*npsi-2] + psi_r[2*npsi-1] * p_ir[0])
                         -zeta_r[2*npsi-1] * p_ir_zetar[2*npsi-1];
  }
}
