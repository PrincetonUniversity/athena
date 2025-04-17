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
//! \file rad_transport.cpp
//  \brief implementation of radiation integrators
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
#include "rad_integrators.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// calculate the transport flux as
// (v- vel)I + vel*I
void RadIntegrator::CalculateFluxes(AthenaArray<Real> &w,
                                    AthenaArray<Real> &ir, const int order) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;
  Coordinates *pco=pmb->pcoord;

  int nang=prad->nang;
  int nfreq=prad->nfreq;
  // Real invcrat=1.0/pmy_rad->crat;

  // int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2,
  // ncells3 = pmb->ncells3;
  // Real tau_fact;

  AthenaArray<Real> &x1flux=prad->flux[X1DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  // int il, iu, jl, ju, kl, ku;
  // jl = js, ju=je, kl=ks, ku=ke;

  // if (ncells2 > 1) {
  //   if (ncells3 == 1) {
  //     jl=js-1, ju=je+1, kl=ks, ku=ke;
  //   } else {
  //     jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  //   }
  // }

  //-----------------------------------------------------------------------------
  // i-direction
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
          GetTaufactor(taul+taur,f_r,1);
          GetTaufactor(taul+taur,f_l,-1);

          Real *s1n = &(sfac1_x_(i,ifr*nang));
          Real *s2n = &(sfac2_x_(i,ifr*nang));
          Real *velxn = &(velx_(k,j,i,ifr*nang));
          Real adv = adv_vel(0,k,j,i);
          SignalSpeed(adv, f_l, f_r, velxn, s1n, s2n);
        }
      }
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
      }

      // calculate flux with HLLE type flux
      //(smax F(u_L) - Smin F(U_R))/(smax-smin)+smax*smin(u(R)-u_L)/(smax-smin)
      for (int i=is; i<=ie+1; ++i) {
        Real *vel = &(velx_(k,j,i,0));
        Real *smax = &(sfac1_x_(i,0));
        Real *smin = &(sfac2_x_(i,0));
        Real *irln = &(il_(i,0));
        Real *irrn = &(ir_(i,0));
        Real adv = 0.0;
        if (adv_flag_ > 0) {
          adv = adv_vel(0,k,j,i);
        }
        Real *sm_diff = &(sm_diff1_(0));
        for (int n=0; n<prad->n_fre_ang; ++n) {
          Real diff = smax[n] - smin[n];
          if (std::abs(diff) < TINY_NUMBER)
            sm_diff[n] = 0.0;
          else
            sm_diff[n] = 1.0/diff;
        }
        for (int n=0; n<prad->n_fre_ang; ++n) {
          Real vl = vel[n] - adv;
          x1flux(k,j,i,n) = smax[n] * (vl - smin[n]) * irln[n] * sm_diff[n]
                            + smin[n] * (smax[n] - vl) * irrn[n] * sm_diff[n];
        }
        if (adv_flag_ > 0) {
          const Real &advv = adv_vel(0,k,j,i);
          if (advv >0) {
            for (int n=0; n<prad->n_fre_ang; ++n) {
              x1flux(k,j,i,n) += advv * irln[n];
            }
          } else {
            for (int n=0; n<prad->n_fre_ang; ++n) {
              x1flux(k,j,i,n) += advv * irrn[n];
            }
          }
        }
      }
    }
  }


  // j-direction
  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux=prad->flux[X2DIR];

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

      // reconstruction
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, is, ie, ir, -1, il_, ir_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, is, ie, ir, -1, il_, ir_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie, ir, -1, il_, ir_);
      }


      for (int j=js; j<=je+1; ++j) {
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, is, ie, ir, -1, ilb_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, is, ie, ir, -1, ilb_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie, ir, -1, ilb_, ir_);
        }

        for (int i=is; i<=ie; ++i) {
          Real *vel = &(vely_(k,j,i,0));
          Real *smax = &(sfac1_y_(j,i,0));
          Real *smin = &(sfac2_y_(j,i,0));
          Real *irln = &(il_(i,0));
          Real *irrn = &(ir_(i,0));
          Real adv = 0.0;
          if (adv_flag_ > 0) {
            adv = adv_vel(1,k,j,i);
          }
          Real *sm_diff = &(sm_diff1_(0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real diff = smax[n] - smin[n];
            if (std::abs(diff) < TINY_NUMBER)
              sm_diff[n] = 0.0;
            else
              sm_diff[n] = 1.0/diff;
          }

          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real vl = vel[n] - adv;
            x2flux(k,j,i,n) = smax[n] * (vl - smin[n]) * irln[n] * sm_diff[n]
                              + smin[n] * (smax[n] - vl) * irrn[n] * sm_diff[n];
            // x2flux(k,j,i,n) = (smax[n] * vl * irln[n] - smin[n] * vl * irrn[n]
            // + smax[n] * smin[n] * (irrn[n] - irln[n]))/(smax[n] - smin[n]);
          }
          if (adv_flag_ > 0) {
            const Real &advv = adv_vel(1,k,j,i);
            if (advv >0) {
              for (int n=0; n<prad->n_fre_ang; ++n) {
                x2flux(k,j,i,n) += advv * irln[n];
              }
            } else {
              for (int n=0; n<prad->n_fre_ang; ++n) {
                x2flux(k,j,i,n) += advv * irrn[n];
              }
            }
          }
        }
        il_.SwapAthenaArray(ilb_);
      }
    }
  }

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux=prad->flux[X3DIR];
    // First, calculate the transport velocity

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

    for (int j=js; j<=je; ++j) { // this loop ordering is intentional
      // reconstruct the first row

      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, is, ie, ir, -1, il_, ir_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, is, ie, ir, -1, il_, ir_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, is, ie, ir, -1, il_, ir_);
      }


      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k

        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, is, ie, ir, -1, ilb_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, is, ie, ir, -1, ilb_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, is, ie, ir, -1, ilb_, ir_);
        }

        for (int i=is; i<=ie; ++i) {
          Real *vel = &(velz_(k,j,i,0));
          Real *smax = &(sfac1_z_(k,j,i,0));
          Real *smin = &(sfac2_z_(k,j,i,0));
          Real *irln = &(il_(i,0));
          Real *irrn = &(ir_(i,0));
          Real adv = 0.0;
          if (adv_flag_ > 0) {
            adv = adv_vel(2,k,j,i);
          }
          Real *sm_diff = &(sm_diff1_(0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real diff = smax[n] - smin[n];
            if (std::abs(diff) < TINY_NUMBER)
              sm_diff[n] = 0.0;
            else
              sm_diff[n] = 1.0/diff;
          }
          for (int n=0; n<prad->n_fre_ang; ++n) {
            Real vl = vel[n] - adv;
            x3flux(k,j,i,n) = smax[n] * (vl - smin[n]) * irln[n] * sm_diff[n]
                              + smin[n] * (smax[n] - vl) * irrn[n] * sm_diff[n];
            // x3flux(k,j,i,n) = (smax[n] * vl * irln[n] - smin[n] * vl * irrn[n]
            // + smax[n] * smin[n] * (irrn[n] - irln[n]))/(smax[n] - smin[n]);
          }// end n
          if (adv_flag_ > 0) {
            const Real &advv = adv_vel(2,k,j,i);
            if (advv >0) {
              for (int n=0; n<prad->n_fre_ang; ++n) {
                x3flux(k,j,i,n) += advv * irln[n];
              }
            } else {
              for (int n=0; n<prad->n_fre_ang; ++n) {
                x3flux(k,j,i,n) += advv * irrn[n];
              }
            }
          }
        }
        // swap the arrays for the next step
        il_.SwapAthenaArray(ilb_);
      }
    }
  }

  // calculate flux along angular direction
  if ((prad->angle_flag == 1) && (prad->nzeta > 0) && (imp_ang_flx_ == 0)) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          // going through all frequency groups
          for (int ifr=0; ifr<nfreq; ++ifr) {
            const int& nzeta = prad->nzeta;
            const int& npsi = prad->npsi;
            int psi_limit=2*npsi;
            if (npsi == 0) psi_limit=1;

            for (int m=0; m<psi_limit; ++m) {
              // TODO(KGF): forced to comment out 4x pragmas here and 2x in
              // ../get_moments.cpp to stop link-time segfault with icpx 2023.0.0
              // and other problems on 2022.2.0
//#pragma omp simd
              for (int n=0; n<nzeta*2; ++n) {
                int ang_num = ifr*nang + n*psi_limit+m;
                q_zeta_(n+NGHOST) = ir(k,j,i,ang_num);
              }
              // Add ghost zones
//#pragma omp simd
              for (int n=1; n<=NGHOST; ++n) {
                q_zeta_(NGHOST-n) = q_zeta_(NGHOST+n-1);
              }
//#pragma omp simd
              for (int n=1; n<=NGHOST; ++n) {
                q_zeta_(2*nzeta+NGHOST+n-1) = q_zeta_(2*nzeta+NGHOST-n);
              }
              int zl = NGHOST-1;
              int zu = 2*nzeta+NGHOST;
              if (order == 1) {
                pmb->precon->DonorCellZeta(prad,zl,zu,q_zeta_,
                                           ql_zeta_,qr_zeta_);
              } else {
                pmb->precon->PiecewiseLinearZeta(prad,zl,zu,q_zeta_,
                                                 ql_zeta_,qr_zeta_);
              }

              // zeta flux
              pco->GetGeometryZeta(prad,k,j,i,g_zeta_);
              for (int n=0; n<nzeta*2+1; ++n) {
                int ang_num = n*psi_limit+m;
                Real g_coef = g_zeta_(n);
                if (g_coef > 0)
                  zeta_flux_(k,j,i,ifr,ang_num) =  -prad->reduced_c * qr_zeta_(n+NGHOST)
                                                   * g_zeta_(n);
                else if (g_coef < 0)
                  zeta_flux_(k,j,i,ifr,ang_num) =  -prad->reduced_c * ql_zeta_(n+NGHOST)
                                                   * g_zeta_(n);
              }
            }
          }
        }
      }
    }
  }

  // Now calculate psi flux
  if ((prad->angle_flag == 1) && (prad->npsi > 0) && (imp_ang_flx_ == 0)) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            const int& nzeta = prad->nzeta;
            const int& npsi = prad->npsi;
            int zeta_limit=2*nzeta;
            if (nzeta == 0) zeta_limit=1;

            for (int n=0; n<zeta_limit; ++n) {
              Real *qpsi = &(q_psi_(NGHOST));
              Real *irm = &(ir(k,j,i,ifr*nang+n*2*npsi));
              for (int m=0; m<npsi*2; ++m) {
                qpsi[m] = irm[m];
              }
              // Add ghost zones
              // phi is periodic
              // the first one qpsi[NGHOST]
              // The last one is qpsi[NGHOST+2*npsi-1]
              for (int m=1; m<=NGHOST; ++m) {
                q_psi_(NGHOST-m) = q_psi_(2*npsi+NGHOST-m);
              }
              for (int m=1; m<=NGHOST; ++m) {
                q_psi_(2*npsi+NGHOST+m-1) = q_psi_(NGHOST+m-1);
              }
              int pl = NGHOST-1;
              int pu = 2*npsi+NGHOST;
              if (order == 1) {
                pmb->precon->DonorCellPsi(prad,pl,pu,q_psi_,
                                          ql_psi_,qr_psi_);
              } else {
                pmb->precon->PiecewiseLinearPsi(prad,pl,pu,q_psi_,
                                                ql_psi_,qr_psi_);
              }

              // psi flux
              if (nzeta > 0)
                pco->GetGeometryPsi(prad,k,j,i,n,g_psi_);
              else
                pco->GetGeometryPsi(prad,k,j,i,g_psi_);
              for (int m=0; m<npsi*2+1; ++m) {
                int ang_num = n*2*npsi+m;
                Real g_coef = g_psi_(m);
                if (g_coef > 0)
                  psi_flux_(k,j,i,ifr,ang_num) =  -prad->reduced_c * qr_psi_(m+NGHOST)
                                                  * g_coef;
                else if (g_coef < 0)
                  psi_flux_(k,j,i,ifr,ang_num) =  -prad->reduced_c * ql_psi_(m+NGHOST)
                                                  * g_coef;
              }
            }
          }
        }
      }
    }
  }
}


// calculate advective flux for the implicit scheme
void RadIntegrator::CalculateFluxes(AthenaArray<Real> &ir, const int order) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;

  // int nang=prad->nang;
  // int nfreq=prad->nfreq;

  // int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2, ncells3 = pmb->ncells3;

  AthenaArray<Real> &x1flux=prad->flux[X1DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if (ncells2 > 1)  {
    if (ncells3 == 1) {
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    } else {
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  //--------------------------------------------------------------------------------------
  // i-direction

  // First, do reconstruction of the specific intensities
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
      }

      // calculate flux with velocity times the interface state
      for (int i=is; i<=ie+1; ++i) {
        const Real &vel = adv_vel(0,k,j,i);
        if (vel >0) {
          for (int n=0; n<prad->n_fre_ang; ++n) {
            x1flux(k,j,i,n) = vel * il_(i,n);
          }
        } else {
          for (int n=0; n<prad->n_fre_ang; ++n) {
            x1flux(k,j,i,n) = vel * ir_(i,n);
          }
        }
      }
    }
  }

  // j-direction
  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux=prad->flux[X2DIR];
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (ncells3 ==  1) // 2D
      kl = ks, ku = ke;
    else // 3D
      kl = ks-1, ku = ke+1;

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, ir, -1, il_, ir_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, ir, -1, il_, ir_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, ir, -1, il_, ir_);
      }

      for (int j=js; j<=je+1; ++j) {
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, ir, -1, ilb_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, ir, -1, ilb_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, ir, -1, ilb_, ir_);
        }

        // calculate flux with velocity times the interface state
        for (int i=il; i<=iu; ++i) {
          const Real &vel = adv_vel(1,k,j,i);
          if (vel > 0) {
            for (int n=0; n<prad->n_fre_ang; ++n) {
              x2flux(k,j,i,n) = vel * il_(i,n);
            }
          } else {
            for (int n=0; n<prad->n_fre_ang; ++n) {
              x2flux(k,j,i,n) = vel * ir_(i,n);
            }
          }
        }
        //swap the array for the next step
        il_.SwapAthenaArray(ilb_);
      }
    }
  }

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux=prad->flux[X3DIR];

    il =is-1, iu=ie+1, jl=js-1, ju=je+1;

    // First, calculate the transport velocity
    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, ir, -1, il_, ir_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, ir, -1, il_, ir_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, ir, -1, il_, ir_);
      }

      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, ir, -1, ilb_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, ir, -1, ilb_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, ir, -1, ilb_, ir_);
        }

        // calculate flux with velocity times the interface state
        for (int i=il; i<=iu; ++i) {
          const Real &vel = adv_vel(2,k,j,i);
          if (vel > 0) {
            for (int n=0; n<prad->n_fre_ang; ++n) {
              x3flux(k,j,i,n) = vel * il_(i,n);
            }
          } else {
            for (int n=0; n<prad->n_fre_ang; ++n) {
              x3flux(k,j,i,n) = vel * ir_(i,n);
            }
          }
        }

        // swap the arrays for the next step
        il_.SwapAthenaArray(ilb_);
      }
    }
  }
}


// add flux divergence due to advection for the implicit scheme
void RadIntegrator::FluxDivergence(const Real wght) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;
  Coordinates *pco= pmb->pcoord;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_;

  AthenaArray<Real> &x1flux=prad->flux[X1DIR];
  AthenaArray<Real> &x2flux=prad->flux[X2DIR];
  AthenaArray<Real> &x3flux=prad->flux[X3DIR];

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pco->CellVolume(k,j,is,ie,vol);
      pco->Face1Area(k,j,is,ie+1,x1area);
      for (int i=is; i<=ie; ++i) {
        Real *divn = &(adv_flx_(k,j,i,0));
        Real *flxr = &(x1flux(k,j,i+1,0));
        Real *flxl = &(x1flux(k,j,i,0));
        Real areal = x1area(i);
        Real arear = x1area(i+1);
        for (int n=0; n<prad->n_fre_ang; ++n) {
          divn[n] = -wght * (arear * flxr[n] - areal * flxl[n])/vol(i);
        }
      }
    }
  }
  if (pmb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pco->CellVolume(k,j,is,ie,vol);
        pco->Face2Area(k,j  ,is,ie,x2area   );
        pco->Face2Area(k,j+1,is,ie,x2area_p1);
        for (int i=is; i<=ie; ++i) {
          Real *divn = &(adv_flx_(k,j,i,0));
          Real *flxr = &(x2flux(k,j+1,i,0));
          Real *flxl = &(x2flux(k,j,i,0));
          Real areal = x2area(i);
          Real arear = x2area_p1(i);
          for (int n=0; n<prad->n_fre_ang; ++n) {
            divn[n] += -wght * (arear * flxr[n] - areal * flxl[n])/vol(i);
          }
        }
      }
    }
  }

  if (pmb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pco->CellVolume(k,j,is,ie,vol);
        pco->Face3Area(k,j  ,is,ie,x3area   );
        pco->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int i=is; i<=ie; ++i) {
          Real *divn = &(adv_flx_(k,j,i,0));
          Real *flxr = &(x3flux(k+1,j,i,0));
          Real *flxl = &(x3flux(k,j,i,0));
          Real areal = x3area(i);
          Real arear = x3area_p1(i);
          for (int n=0; n<prad->n_fre_ang; ++n) {
            divn[n] += -wght * (arear * flxr[n] - areal * flxl[n])/vol(i);
          }
        }
      }
    }
  }
}


void RadIntegrator::FluxDivergence(const Real wght, AthenaArray<Real> &ir_in,
                                   AthenaArray<Real> &ir_out) {
  NRRadiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;
  int nfreq=prad->nfreq;
  int nang=prad->nang;

  AthenaArray<Real> &x1flux=prad->flux[X1DIR];
  AthenaArray<Real> &x2flux=prad->flux[X2DIR];
  AthenaArray<Real> &x3flux=prad->flux[X3DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  AthenaArray<Real> &area_zeta = zeta_area_, &area_psi = psi_area_,
                      &ang_vol = ang_vol_, &dflx_ang = dflx_ang_;
  const int& nzeta = prad->nzeta;
  const int& npsi = prad->npsi;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for (int i=is; i<=ie; ++i) {
        Real *flxr = &(x1flux(k,j,i+1,0));
        Real *flxl = &(x1flux(k,j,i,0));
        Real *flxn = &(dflx(i,0));
        for (int n=0; n<prad->n_fre_ang; ++n) {
          flxn[n] = (x1area(i+1) *flxr[n] - x1area(i)*flxl[n]);
        }
      }

      // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for (int i=is; i<=ie; ++i) {
          Real *flxr = &(x2flux(k,j+1,i,0));
          Real *flxl = &(x2flux(k,j,i,0));
          Real *flxn = &(dflx(i,0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            flxn[n] += (x2area_p1(i)*flxr[n] - x2area(i)*flxl[n]);
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int i=is; i<=ie; ++i) {
          Real *flxr = &(x3flux(k+1,j,i,0));
          Real *flxl = &(x3flux(k,j,i,0));
          Real *flxn = &(dflx(i,0));
          for (int n=0; n<prad->n_fre_ang; ++n) {
            flxn[n] += (x3area_p1(i)*flxr[n] - x3area(i)*flxl[n]);
          }
        }
      }
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for (int i=is; i<=ie; ++i) {
        Real *irin = &(ir_in(k,j,i,0));
        Real *iro = &(ir_out(k,j,i,0));
        Real *flxn = &(dflx(i,0));
        for (int n=0; n<prad->n_fre_ang; ++n) {
          iro[n] = std::max(irin[n]-wght*flxn[n]/vol(i), static_cast<Real>(TINY_NUMBER));
        }
      }

      // add angular flux
      if ((prad->angle_flag == 1) && (imp_ang_flx_ == 0)) {
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for (int n=0; n<prad->nang; ++n)
              dflx_ang(n) = 0.0;
            if (nzeta * npsi > 0) {
              for (int m=0; m<2*npsi; ++m) {
//#pragma omp simd
                for (int n=0; n<2*nzeta; ++n) {
                  int ang_num = n*2*npsi + m;
                  int ang_num1 = (n+1)*2*npsi+m;
                  dflx_ang(ang_num) += (area_zeta(m,n+1) * zeta_flux_(k,j,i,ifr,ang_num1)
                                        - area_zeta(m,n) * zeta_flux_(k,j,i,ifr,ang_num));
                }
              }
              // now psi flux
              for (int n=0; n<2*nzeta; ++n) {
                Real *flxn = &(dflx_ang(n*2*npsi));
                Real *areapsi = &(area_psi(n,0));
                Real *psiflx = &(psi_flux_(k,j,i,ifr,n*2*npsi));
                for (int m=0; m<2*npsi; ++m) {
                  flxn[m] += (areapsi[m+1] * psiflx[m+1]
                              - areapsi[m] * psiflx[m]);
                }
              }
            } else if (nzeta >0) {
              Real *flxn = &(dflx_ang(0));
              Real *areazeta = &(area_zeta(0));
              Real *zetaflx = &(zeta_flux_(k,j,i,ifr,0));
              for (int n=0; n<2*nzeta; ++n) {
                flxn[n] += (areazeta[n+1] * zetaflx[n+1]
                            - areazeta[n] * zetaflx[n]);
              }
            } else if (npsi > 0) {
              Real *flxn = &(dflx_ang(0));
              Real *areapsi = &(area_psi(0));
              Real *psiflx = &(psi_flux_(k,j,i,ifr,0));
              for (int m=0; m<2*npsi; ++m) {
                flxn[m] += (areapsi[m+1] * psiflx[m+1]
                            - areapsi[m] * psiflx[m]);
              }
            }
            // apply the flux divergence back
            // We need ir_out, not ir_in, as ir_out is already partially updated
            Real *iro = &(ir_out(k,j,i,ifr*nang));
            Real *flxn = &(dflx_ang(0));
            Real *angv = &(ang_vol(0));
            for (int n=0; n<prad->nang; ++n) {
              iro[n] = std::max(iro[n]-wght*flxn[n]/angv[n],
                                static_cast<Real>(TINY_NUMBER));
            }
          }
        }
      }
    }
  }
  if (prad->angle_flag == 1 && (imp_ang_flx_ == 1)) {
    ImplicitAngularFluxesCoef(wght);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          ImplicitAngularFluxes(k,j,i,ir_out);

          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real *p_angflx  = &(ang_flx_(k,j,i,ifr*nang));
            Real *iro = &(ir_out(k,j,i,ifr*nang));
            Real *ang_coef = &(imp_ang_coef_(k,j,i,0));
            for (int n=0; n<nang; ++n) {
              iro[n] += p_angflx[n];
              iro[n] /= (1.0 + ang_coef[n]);
            }
          }
        }
      }
    }
  }
}
