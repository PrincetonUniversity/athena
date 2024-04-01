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
//! \file implicit_first_order.cpp
//  \brief implementation of implicit transport with first order reconstruction
//======================================================================================

// C headers

// C++ headers

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


// the flux is
// Smax F(U_L)/(Smax-Smin) - Smin F(U_R) /(Smax-Smin)
// + Smax Smin * (U_R - U_L)/(Smax - Smin)
// Gauss-Seidel iteration, use upper triangle from last iteration
// lower triangle is from next iteration

void RadIntegrator::FirstOrderFluxDivergenceCoef(const Real wght) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;
  Coordinates *pco= pmb->pcoord;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int nang = prad->nang;
  int nfreq = prad->nfreq;

  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pco->CenterWidth1(k,j,is-1,ie+1,dxw1_);
      for (int i=is; i<=ie+1; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          // use the signal speed in each frequency group
          Real sigmal = prad->sigma_a(k,j,i-1,ifr) + prad->sigma_s(k,j,i-1,ifr);
          Real sigmar = prad->sigma_a(k,j,i,ifr) + prad->sigma_s(k,j,i,ifr);
          Real taul = dxw1_(i-1) * sigmal;
          Real taur = dxw1_(i) * sigmar;

          Real f_l = 1.0;
          Real f_r = 1.0;
          taul *= taufact(k,j,i-1);
          taur *= taufact(k,j,i);
          // the choice of taufactor = sum of L and R cell optical depths is made after
          // many trial experiments, particularly to minimize artifacts when there is a
          // sharp opacity jump. other choices such as arithmetic and harmonic mean may
          // work well in some applications
          GetTaufactor(taul+taur,f_r,1);
          GetTaufactor(taul+taur,f_l,-1);

          Real *s1n = &(sfac1_x_(i,ifr*nang));
          Real *s2n = &(sfac2_x_(i,ifr*nang));
          Real *velxn = &(velx_(k,j,i,ifr*nang));
          Real adv = adv_vel(0,k,j,i);
          SignalSpeed(adv, f_l, f_r, velxn, s1n, s2n);
        }
      }
      // calculate x1-flux divergence
      pco->Face1Area(k,j,is,ie+1,x1area);
      for (int i=is; i<=ie; ++i) {
        Real areal = x1area(i);
        Real arear = x1area(i+1);
        Real *smax_l = &(sfac1_x_(i,0));
        Real *smax_r = &(sfac1_x_(i+1,0));
        Real *smin_l = &(sfac2_x_(i,0));
        Real *smin_r = &(sfac2_x_(i+1,0));
        Real *vel_ln = &(velx_(k,j,i,0));
        Real *vel_rn = &(velx_(k,j,i+1,0));

        Real *coef1n = &(const_coef_(k,j,i,0));
        Real *coef1ln= &(const_coef1_l_(k,j,i,0));
        Real *coef1rn= &(const_coef1_r_(k,j,i,0));
        Real *exp_coefn=&(exp_coef_(k,j,i,0));

        Real advl = 0.0;
        Real advr = 0.0;
        if (adv_flag_ > 0) {
          advl = adv_vel(0,k,j,i);
          advr = adv_vel(0,k,j,i+1);
        }
        Real *sm_diff1 = &(sm_diff1_(0));
        Real *sm_diff2 = &(sm_diff2_(0));
        for (int n=0; n<prad->n_fre_ang; ++n) {
          Real diff = smax_l[n] - smin_l[n];
          if (std::abs(diff) < TINY_NUMBER)
            sm_diff1[n] = 0.0;
          else
            sm_diff1[n] = 1.0/diff;
        }
        for (int n=0; n<prad->n_fre_ang; ++n) {
          Real diff = smax_r[n] - smin_r[n];
          if (std::abs(diff) < TINY_NUMBER)
            sm_diff2[n] = 0.0;
          else
            sm_diff2[n] = 1.0/diff;
        }

        SplitVelocity(vel_ln, vel_rn,advl, advr,smax_l,smin_l,smax_r,smin_r);

        for (int n=0; n<prad->n_fre_ang; ++n) {
          Real vl = vel_ln[n] - advl;
          coef1n[n] = -areal * vel_im_l_(n) * sm_diff1[n];
          coef1ln[n] = -areal * smax_l[n] * (vl - smin_l[n]) * sm_diff1[n];
          exp_coefn[n] = -areal * vel_ex_l_(n) * sm_diff1[n];

          Real vr = vel_rn[n] - advr;
          coef1n[n] += arear * vel_im_r_(n) * sm_diff2[n];
          coef1rn[n] = arear * smin_r[n] * (smax_r[n] - vr) * sm_diff2[n];
          exp_coefn[n] += arear * vel_ex_r_(n) * sm_diff2[n];
        }
      }
    }
  }
  // calculate x2-flux
  if (pmb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k) {
      // first, calculate speed
      for (int j=js; j<=je+1; ++j) {
        pco->CenterWidth2(k,j-1,is,ie,dxw1_);
        pco->CenterWidth2(k,j,is,ie,dxw2_);
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real sigmal = prad->sigma_a(k,j-1,i,ifr) + prad->sigma_s(k,j-1,i,ifr);
            Real sigmar = prad->sigma_a(k,j,i,ifr) + prad->sigma_s(k,j,i,ifr);
            Real taul = dxw1_(i) * sigmal;
            Real taur = dxw2_(i) * sigmar;

            Real f_l = 1.0;
            Real f_r = 1.0;
            taul *= taufact(k,j-1,i);
            taur *= taufact(k,j,i);
            GetTaufactor(taul+taur,f_r,1);
            GetTaufactor(taul+taur,f_l,-1);

            Real *s1n = &(sfac1_y_(j,i,ifr*nang));
            Real *s2n = &(sfac2_y_(j,i,ifr*nang));
            Real *velyn = &(vely_(k,j,i,ifr*nang));
            Real adv = adv_vel(1,k,j,i);
            SignalSpeed(adv, f_l, f_r, velyn, s1n, s2n);
          }
        }
      }
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);

        for (int i=is; i<=ie; ++i) {
          Real *vel_ln = &(vely_(k,j,i,0));
          Real *vel_rn = &(vely_(k,j+1,i,0));
          Real *coef2n = &(const_coef_(k,j,i,0));
          Real *coef2ln= &(const_coef2_l_(k,j,i,0));
          Real *coef2rn= &(const_coef2_r_(k,j,i,0));
          Real *exp_coefn=&(exp_coef_(k,j,i,0));

          Real *smax_l = &(sfac1_y_(j,i,0));
          Real *smax_r = &(sfac1_y_(j+1,i,0));
          Real *smin_l = &(sfac2_y_(j,i,0));
          Real *smin_r = &(sfac2_y_(j+1,i,0));

          Real areal = x2area(i);
          Real arear = x2area_p1(i);

          Real advl = 0.0;
          Real advr = 0.0;
          if (adv_flag_ > 0) {
            advl = adv_vel(1,k,j,i);
            advr = adv_vel(1,k,j+1,i);
          }
          Real *sm_diff1 = &(sm_diff1_(0));
          Real *sm_diff2 = &(sm_diff2_(0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real diff = smax_l[n] - smin_l[n];
            if (std::abs(diff) < TINY_NUMBER)
              sm_diff1[n] = 0.0;
            else
              sm_diff1[n] = 1.0/diff;
          }
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real diff = smax_r[n] - smin_r[n];
            if (std::abs(diff) < TINY_NUMBER)
              sm_diff2[n] = 0.0;
            else
              sm_diff2[n] = 1.0/diff;
          }

          SplitVelocity(vel_ln, vel_rn,advl, advr,smax_l,smin_l,smax_r,smin_r);

          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real vl = vel_ln[n] - advl;
            coef2n[n] += -areal * vel_im_l_(n) * sm_diff1[n];
            coef2ln[n] = -areal * smax_l[n] * (vl - smin_l[n]) * sm_diff1[n];
            exp_coefn[n] += -areal * vel_ex_l_(n) * sm_diff1[n];

            Real vr = vel_rn[n] - advr;
            coef2n[n] += arear * vel_im_r_(n) * sm_diff2[n];
            coef2rn[n] = arear * smin_r[n] * (smax_r[n] - vr) * sm_diff2[n];
            exp_coefn[n] += arear * vel_ex_r_(n) * sm_diff2[n];
          }
        }
      }
    }
  }

  // calculate x3-flux divergence
  if (pmb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pco->CenterWidth3(k-1,j,is,ie,dxw1_);
        pco->CenterWidth3(k,j,is,ie,dxw2_);
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real sigmal = prad->sigma_a(k-1,j,i,ifr) + prad->sigma_s(k-1,j,i,ifr);
            Real sigmar = prad->sigma_a(k,j,i,ifr) + prad->sigma_s(k,j,i,ifr);
            Real taul = dxw1_(i) * sigmal;
            Real taur = dxw2_(i) * sigmar;

            Real f_l = 1.0;
            Real f_r = 1.0;
            taul *= taufact(k-1,j,i);
            taur *= taufact(k,j,i);
            GetTaufactor(taul+taur,f_r,1);
            GetTaufactor(taul+taur,f_l,-1);

            Real *s1n = &(sfac1_z_(k,j,i,ifr*nang));
            Real *s2n = &(sfac2_z_(k,j,i,ifr*nang));
            Real *velzn = &(velz_(k,j,i,ifr*nang));
            Real adv = adv_vel(2,k,j,i);
            SignalSpeed(adv, f_l, f_r, velzn, s1n, s2n);
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int i=is; i<=ie; ++i) {
          Real *smax_l = &(sfac1_z_(k,j,i,0));
          Real *smax_r = &(sfac1_z_(k+1,j,i,0));
          Real *smin_l = &(sfac2_z_(k,j,i,0));
          Real *smin_r = &(sfac2_z_(k+1,j,i,0));
          Real *vel_ln = &(velz_(k,j,i,0));
          Real *vel_rn = &(velz_(k+1,j,i,0));
          Real *coef3n = &(const_coef_(k,j,i,0));
          Real *coef3ln= &(const_coef3_l_(k,j,i,0));
          Real *coef3rn= &(const_coef3_r_(k,j,i,0));
          Real *exp_coefn=&(exp_coef_(k,j,i,0));
          Real areal = x3area(i);
          Real arear = x3area_p1(i);
          Real advl = 0.0;
          Real advr = 0.0;
          if (adv_flag_ > 0) {
            advl = adv_vel(2,k,j,i);
            advr = adv_vel(2,k+1,j,i);
          }
          Real *sm_diff1 = &(sm_diff1_(0));
          Real *sm_diff2 = &(sm_diff2_(0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real diff = smax_l[n] - smin_l[n];
            if (std::abs(diff) < TINY_NUMBER)
              sm_diff1[n] = 0.0;
            else
              sm_diff1[n] = 1.0/diff;
          }
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real diff = smax_r[n] - smin_r[n];
            if (std::abs(diff) < TINY_NUMBER)
              sm_diff2[n] = 0.0;
            else
              sm_diff2[n] = 1.0/diff;
          }

          SplitVelocity(vel_ln, vel_rn,advl, advr,smax_l,smin_l,smax_r,smin_r);

          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real vl = vel_ln[n] - advl;
            coef3n[n] += -areal * vel_im_l_(n) * sm_diff1[n];
            coef3ln[n] = -areal * smax_l[n] * (vl - smin_l[n]) * sm_diff1[n];
            exp_coefn[n] += -areal * vel_ex_l_(n) * sm_diff1[n];

            Real vr = vel_rn[n] - advr;
            coef3n[n] += arear * vel_im_r_(n) * sm_diff2[n];
            coef3rn[n] = arear * smin_r[n] * (smax_r[n] - vr) * sm_diff2[n];
            exp_coefn[n] += arear * vel_ex_r_(n) * sm_diff2[n];
          }
        }
      }
    }
  }

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for (int i=is; i<=ie; ++i) {
        Real dtvol = wght/vol(i);
        Real *coef1n = &(const_coef_(k,j,i,0));
        Real *exp_coefn=&(exp_coef_(k,j,i,0));
        Real *coef1ln=&(const_coef1_l_(k,j,i,0));
        Real *coef1rn=&(const_coef1_r_(k,j,i,0));
        for (int n=0; n<prad->n_fre_ang; ++n) {
          coef1n[n] *= dtvol;
          exp_coefn[n] *= dtvol;
          coef1ln[n] *= dtvol;
          coef1rn[n] *= dtvol;
        }
        if (pmb->block_size.nx2 > 1) {
          Real *coef2ln = &(const_coef2_l_(k,j,i,0));
          Real *coef2rn = &(const_coef2_r_(k,j,i,0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            coef2ln[n] *= dtvol;
            coef2rn[n] *= dtvol;
          }
        }
        if (pmb->block_size.nx3 > 1) {
          Real *coef3ln = &(const_coef3_l_(k,j,i,0));
          Real *coef3rn = &(const_coef3_r_(k,j,i,0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            coef3ln[n] *= dtvol;
            coef3rn[n] *= dtvol;
          }
        }
      }
    }
  }
  return;
}

void RadIntegrator::FirstOrderFluxDivergence(const int k, const int j, const int i,
                                                          AthenaArray<Real> &ir_ini) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;

  Real *divn = &(divflx_(k,j,i,0));
  Real *irn = &(ir_ini(k,j,i,0));
  Real *irln = &(ir_ini(k,j,i-1,0));
  Real *irrn = &(ir_ini(k,j,i+1,0));
  Real *exp_coefn = &(exp_coef_(k,j,i,0));
  Real *coef1ln = &(const_coef1_l_(k,j,i,0));
  Real *coef1rn = &(const_coef1_r_(k,j,i,0));

  for (int n=0; n<prad->n_fre_ang; ++n) {
    divn[n] = -(coef1ln[n] * irln[n] + coef1rn[n] * irrn[n]
              + exp_coefn[n] * irn[n]);
  }
  if (adv_flag_ > 0) {
    Real *advn = &(adv_flx_(k,j,i,0));
    for (int n=0; n<prad->n_fre_ang; ++n) {
      divn[n] += advn[n];
    }
  }
  if (pmb->block_size.nx2 > 1) {
    Real *irjln = &(ir_ini(k,j-1,i,0));
    Real *irjrn = &(ir_ini(k,j+1,i,0));
    Real *coef2ln = &(const_coef2_l_(k,j,i,0));
    Real *coef2rn = &(const_coef2_r_(k,j,i,0));
    for (int n=0; n<prad->n_fre_ang; ++n) {
      Real temp = -(coef2ln[n] * irjln[n] + coef2rn[n] * irjrn[n]);
      divn[n] += temp;
    }
  }

  if (pmb->block_size.nx3 > 1) {
    Real *irkln = &(ir_ini(k-1,j,i,0));
    Real *irkrn = &(ir_ini(k+1,j,i,0));
    Real *coef3ln = &(const_coef3_l_(k,j,i,0));
    Real *coef3rn = &(const_coef3_r_(k,j,i,0));
    for (int n=0; n<prad->n_fre_ang; ++n) {
      Real temp = -(coef3ln[n] * irkln[n] + coef3rn[n] * irkrn[n]);
      divn[n] += temp;
    }
  }
  return;
}
