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
//! \file rad_integrators.cpp
//  \brief implementation of radiation integrators
//======================================================================================

// C headers

// C++ headers
#include <algorithm>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../implicit/radiation_implicit.hpp"
#include "../radiation.hpp"
#include "rad_integrators.hpp"

static constexpr Real tau_asymptotic_lim=1.0e-3;

RadIntegrator::RadIntegrator(NRRadiation *prad, ParameterInput *pin) {
  pmy_rad = prad;

  MeshBlock *pmb = prad->pmy_block;
  Coordinates *pco=pmb->pcoord;

  int nang=prad->nang;
  int nfreq=prad->nfreq;
  rad_xorder = pin->GetOrAddInteger("time","rad_xorder",2);
  if (rad_xorder == 3) {
    if (NGHOST < 3) {
      std::stringstream msg;
      msg << "### FATAL ERROR in radiation reconstruction constructor" << std::endl
          << "rad_xorder=" << rad_xorder <<
          " (PPM) reconstruction selected, but nghost=" << NGHOST << std::endl
          << "Reconfigure with --nghost=3  " <<std::endl;
      ATHENA_ERROR(msg);
    }
  }


  // factor to separate the diffusion and advection part
  Real taucell = pin->GetOrAddReal("radiation","taucell",5);
  tau_flag_ = pin->GetOrAddInteger("radiation","tau_scheme",1);
  compton_flag_=pin->GetOrAddInteger("radiation","Compton",0);
  compton_t_=pin->GetOrAddInteger("radiation","Compton_t",0);
  split_compton_=pin->GetOrAddInteger("radiation","Split_compton",0);
  if ((!IM_RADIATION_ENABLED) || (compton_flag_ == 0))
    split_compton_ = 0;

  adv_flag_=pin->GetOrAddInteger("radiation","Advection",1);
  flux_correct_flag_ = pin->GetOrAddInteger("radiation","CorrectFlux",0);
  imp_ang_flx_ = pin->GetOrAddInteger("radiation","implicit_ang_flux",0);

  doppler_flag_ = pin->GetOrAddInteger("radiation","doppler_flag",1);

  // multi-group iteration
  // iterations < iterative_tgas_
  iteration_tgas_  = pin->GetOrAddInteger("radiation","iterative_tgas",5);
  iteration_compton_  = pin->GetOrAddInteger("radiation","iterative_compton",5);
  tgas_error_ =   pin->GetOrAddReal("radiation","gas_error",1.e-6);
  compton_error_ =   pin->GetOrAddReal("radiation","gas_error",1.e-6);
  // maximum number of bins each frequency bin will map to, default is nfreq/2
  nmax_map_ = pin->GetOrAddInteger("radiation","max_map_bin",0);

  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2,
      ncells3 = pmb->ncells3;



  x1face_area_.NewAthenaArray(ncells1+1);
  if (ncells2 > 1) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if (ncells3 > 1) {
    x3face_area_.NewAthenaArray(ncells1);
    x3face_area_p1_.NewAthenaArray(ncells1);
  }
  cell_volume_.NewAthenaArray(ncells1);


  cwidth2_.NewAthenaArray(ncells1);
  cwidth3_.NewAthenaArray(ncells1);

  dflx_.NewAthenaArray(ncells1,prad->n_fre_ang);

  // arrays for spatial recontruction
  il_.NewAthenaArray(ncells1,prad->n_fre_ang);
  ilb_.NewAthenaArray(ncells1,prad->n_fre_ang);

  ir_.NewAthenaArray(ncells1,prad->n_fre_ang);

  sfac1_x_.NewAthenaArray(ncells1,prad->n_fre_ang);
  sfac2_x_.NewAthenaArray(ncells1,prad->n_fre_ang);
  if (ncells2 > 1) {
    sfac1_y_.NewAthenaArray(ncells2,ncells1,prad->n_fre_ang);
    sfac2_y_.NewAthenaArray(ncells2,ncells1,prad->n_fre_ang);
  }
  if (ncells3 > 1) {
    sfac1_z_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    sfac2_z_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  }



  sm_diff1_.NewAthenaArray(prad->n_fre_ang);
  sm_diff2_.NewAthenaArray(prad->n_fre_ang);
  vel_ex_l_.NewAthenaArray(prad->n_fre_ang);
  vel_ex_r_.NewAthenaArray(prad->n_fre_ang);
  vel_im_l_.NewAthenaArray(prad->n_fre_ang);
  vel_im_r_.NewAthenaArray(prad->n_fre_ang);

  adv_vel.NewAthenaArray(3,ncells3,ncells2,ncells1);


  if (IM_RADIATION_ENABLED) {
    limiter_.NewAthenaArray(ncells1,prad->n_fre_ang);

    if (ncells2 > 1) {
      limiterj_.NewAthenaArray(ncells2,ncells1,prad->n_fre_ang);
    }
    if (ncells3 > 1) {
      limiterk_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    }
    dql_.NewAthenaArray(prad->n_fre_ang);
    dqr_.NewAthenaArray(prad->n_fre_ang);

    const_coef_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    exp_coef_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    const_coef1_l_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    const_coef1_r_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    if (ncells2 > 1) {
      const_coef2_l_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
      const_coef2_r_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    }
    if (ncells3 > 1) {
      const_coef3_l_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
      const_coef3_r_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    }
    divflx_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);

    left_coef1_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    left_coef2_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    left_coef3_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    adv_flx_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  }

  implicit_coef_.NewAthenaArray(prad->n_fre_ang);

  tgas_.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas_new_.NewAthenaArray(ncells3,ncells2,ncells1);
  vel_source_.NewAthenaArray(ncells3,ncells2,ncells1,3);
  taufact.NewAthenaArray(ncells3,ncells2,ncells1);


  rad_source.NewAthenaArray(4,ncells3,ncells2,ncells1);
  delta_source.NewAthenaArray(4,nfreq);


  vel_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  velx_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  vely_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  velz_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);

  vncsigma_.NewAthenaArray(nang);
  vncsigma2_.NewAthenaArray(prad->n_fre_ang);
  wmu_cm_.NewAthenaArray(nang);
  tran_coef_.NewAthenaArray(nang);
  cm_to_lab_.NewAthenaArray(nang);
  ir_cm_.NewAthenaArray(prad->n_fre_ang);
  dxw1_.NewAthenaArray(ncells1);
  dxw2_.NewAthenaArray(ncells1);

  //----------------------------------------------------
  // array for multi-group
  if (nfreq > 1) {
    // put the angle index first when map in frequency space

    ir_buff_.NewAthenaArray(prad->n_fre_ang);

    if (nmax_map_ == 0)  nmax_map_ = nfreq/2+1;
    split_ratio_.NewAthenaArray(nang,nfreq,nmax_map_);
    // this needs to be done for each angle, all frequency groups
    fre_map_matrix_.NewAthenaArray(nfreq,nmax_map_);
    delta_nu_n_.NewAthenaArray(nfreq);
    map_bin_start_.NewAthenaArray(nang,nfreq);
    map_bin_end_.NewAthenaArray(nang,nfreq);
    map_count_.NewAthenaArray(nfreq);
    nu_shift_.NewAthenaArray(nfreq);
    ir_face_.NewAthenaArray(nfreq);
    ir_ori_.NewAthenaArray(nfreq);
    ir_done_.NewAthenaArray(nfreq);

    com_b_face_coef_.NewAthenaArray(nfreq);
    com_d_face_coef_.NewAthenaArray(nfreq);
    com_b_coef_l_.NewAthenaArray(nfreq);
    com_b_coef_r_.NewAthenaArray(nfreq);
    com_d_coef_l_.NewAthenaArray(nfreq);
    com_d_coef_r_.NewAthenaArray(nfreq);

    nf_rhs_.NewAthenaArray(nfreq);
    nf_n0_.NewAthenaArray(nfreq);
    new_j_nu_.NewAthenaArray(nfreq);
  }

  sum_nu3_.NewAthenaArray(nfreq);
  sum_nu2_.NewAthenaArray(nfreq);
  sum_nu1_.NewAthenaArray(nfreq);
  eq_sol_.NewAthenaArray(nfreq);

  //----------------------------------------------------
  // array for multi-group

  // set the default taufact
  for (int k=0; k<ncells3; ++k)
    for (int j=0; j<ncells2; ++j)
      for (int i=0; i<ncells1; ++i) {
        taufact(k,j,i) = taucell;
      }

  // The angular grid coefficients
  ang_flx_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  imp_ang_coef_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  imp_ang_coef_r_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  imp_ang_psi_r_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  imp_ang_psi_l_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);

  if (prad->angle_flag == 1) {
    const int& nzeta = prad->nzeta;
    const int& npsi = prad->npsi;
    if (nzeta > 0) {
      g_zeta_.NewAthenaArray(2*nzeta+1);
      q_zeta_.NewAthenaArray(2*nzeta+2*NGHOST);
      ql_zeta_.NewAthenaArray(2*nzeta+2*NGHOST);
      qr_zeta_.NewAthenaArray(2*nzeta+2*NGHOST);


      if (npsi > 0) {
        zeta_flux_.NewAthenaArray(ncells3,ncells2,ncells1,nfreq,(2*nzeta+1)*2*npsi);
        zeta_area_.NewAthenaArray(2*npsi,2*nzeta+1);
      } else {
        zeta_flux_.NewAthenaArray(ncells3,ncells2,ncells1,nfreq,2*nzeta+1);
        zeta_area_.NewAthenaArray(2*nzeta+1);
      }
      pco->ZetaArea(prad, zeta_area_);
    }

    if (npsi > 0) {
      g_psi_.NewAthenaArray(2*npsi+1);
      q_psi_.NewAthenaArray(2*npsi+2*NGHOST);
      ql_psi_.NewAthenaArray(2*npsi+2*NGHOST);
      qr_psi_.NewAthenaArray(2*npsi+2*NGHOST);
      if (nzeta > 0) {
        psi_flux_.NewAthenaArray(ncells3,ncells2,ncells1,nfreq,2*nzeta*(2*npsi+1));
        psi_area_.NewAthenaArray(2*nzeta,2*npsi+1);
      } else {
        psi_flux_.NewAthenaArray(ncells3,ncells2,ncells1,nfreq,2*npsi+1);
        psi_area_.NewAthenaArray(2*npsi+1);
      }
      pco->PsiArea(prad, psi_area_);
    }

    dflx_ang_.NewAthenaArray(nang);
    ang_vol_.NewAthenaArray(nang);
    pco->AngularVol(prad, ang_vol_);
  }

  // calculate the advection velocity at the cell faces
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if (ncells2 > 1) {
    if (ncells3 == 1) {
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    } else {
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }

  // calculate velx_
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // get the velocity at the interface
      for (int i=is-1; i<=ie+1; ++i) {
        Real dxl = pco->x1f(i)-pco->x1v(i-1);
        Real dxr = pco->x1v(i) - pco->x1f(i);
        Real factl = dxr/(dxl+dxr);
        Real factr = dxl/(dxl+dxr);
        for (int ifr=0; ifr<nfreq; ++ifr) {
          Real *cosx = &(prad->mu(0,k,j,i-1,0));
          Real *cosx1 = &(prad->mu(0,k,j,i,0));
          Real *veln = &(velx_(k,j,i,ifr*nang));
          for (int n=0; n<nang; ++n) {
            // linear intepolation between x1v(i-1), x1f(i), x1v(i)
            veln[n] = prad->reduced_c *
                      (factl * cosx[n] + factr * cosx1[n]);
          }
        }
      }
    }
  }

  // calculate vely_
  if (ncells2 > 1) {
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (ncells3 >  1) // 2D
      kl = ks-1, ku = ke+1;

    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
        // get the velocity at the interface
        for (int i=il; i<=iu; ++i) {
          Real dxl = pco->x2f(j)-pco->x2v(j-1);
          Real dxr = pco->x2v(j) - pco->x2f(j);
          Real factl = dxr/(dxl+dxr);
          Real factr = dxl/(dxl+dxr);
          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real *cosy = &(prad->mu(1,k,j-1,i,0));
            Real *cosy1 = &(prad->mu(1,k,j,i,0));
            Real *veln = &(vely_(k,j,i,ifr*nang));
            for (int n=0; n<nang; ++n) {
              // linear intepolation between x2v(j-1), x2f(j), x2v(j)
              veln[n] = prad->reduced_c *
                        (factl * cosy[n] + factr * cosy1[n]);
            }
          }
        }
      }
    }
  }

  // calculate vely_
  if (ncells3 > 1) {
    il =is-1, iu=ie+1, jl=js-1, ju=je+1;

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // get the velocity at the interface
        for (int i=il; i<=iu; ++i) {
          Real dxl = pco->x3f(k) - pco->x3v(k-1);
          Real dxr = pco->x3v(k) - pco->x3f(k);
          Real factl = dxr/(dxl+dxr);
          Real factr = dxl/(dxl+dxr);
          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real *cosz = &(prad->mu(2,k-1,j,i,0));
            Real *cosz1 = &(prad->mu(2,k,j,i,0));
            Real *veln = &(velz_(k,j,i,ifr*nang));
            for (int n=0; n<nang; ++n) {
              // linear intepolation between x2v(j-1), x2f(j), x2v(j)
              veln[n] = prad->reduced_c *
                        (factl * cosz[n] + factr * cosz1[n]);
            }
          }
        }
      }
    }
  }
}

void RadIntegrator::GetTgasVel(MeshBlock *pmb, const Real dt,
                               AthenaArray<Real> &u, AthenaArray<Real> &w,
                               AthenaArray<Real> &bcc, AthenaArray<Real> &ir) {
  Real gm1 = pmb->peos->GetGamma() - 1.0;

  Real rho_floor = pmb->peos->GetDensityFloor();

  NRRadiation *prad=pmb->pnrrad;
  Coordinates *pco=pmb->pcoord;

  const Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;

  const int &nang =prad->nang;
  const int &nfreq=prad->nfreq;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=0; k<pmb->ncells3; ++k) {
    for (int j=0; j<pmb->ncells2; ++j) {
      for (int i=0; i<pmb->ncells1; ++i) {
        // for implicit update, using the quantities from the partially
        // updated u, not from w
        Real rho = u(IDN,k,j,i);
        rho = std::max(rho,rho_floor);
        Real vx = u(IM1,k,j,i)/rho;
        Real vy = u(IM2,k,j,i)/rho;
        Real vz = u(IM3,k,j,i)/rho;
        Real pb = 0.0;
        if (MAGNETIC_FIELDS_ENABLED)
          pb = 0.5*(SQR(bcc(IB1,k,j,i))+SQR(bcc(IB2,k,j,i))
                    +SQR(bcc(IB3,k,j,i)));

        Real vsq = vx * vx + vy * vy + vz * vz;
        Real tgas = u(IEN,k,j,i) - pb - 0.5*rho*vsq;
        tgas = gm1*tgas/rho;
        tgas = std::max(tgas,pmb->pnrrad->t_floor_(k,j,i));
        tgas = std::min(tgas,pmb->pnrrad->t_ceiling_(k,j,i));
        tgas_(k,j,i) = tgas;
        // Do not use the velocity directly in strongly radiation pressure
        // dominated regime
        // use the predicted velocity based on moment equatio

        // calculate radiation energy density
        Real er = 0.0;
        for (int ifr=0; ifr<nfreq; ++ifr) {
          Real *irn = &(ir(k,j,i,ifr*nang));
          Real *weight = &(prad->wmu(0));
          Real er_freq = 0.0;
#pragma omp simd reduction(+:er_freq)
          for (int n=0; n<nang; ++n) {
            er_freq += weight[n] * irn[n];
          }
          er += er_freq;
        }

        // now the velocity term,
        // using velocity from current stage
        vx = w(IVX,k,j,i);
        vy = w(IVY,k,j,i);
        vz = w(IVZ,k,j,i);
        vsq = vx * vx + vy * vy + vz * vz;
        if (prat * er * invcrat * invcrat > rho) {
          PredictVel(ir,k,j,i, 0.5 * dt, rho, &vx, &vy, &vz);
          vsq = vx * vx + vy * vy + vz * vz;
        }

        Real ratio = std::sqrt(vsq) * invcrat;
        // Limit the velocity to be smaller than the speed of light
        if (ratio > prad->vmax) {
          Real factor = prad->vmax/ratio;
          vx *= factor;
          vy *= factor;
          vz *= factor;
        }
        vel_source_(k,j,i,0) = vx;
        vel_source_(k,j,i,1) = vy;
        vel_source_(k,j,i,2) = vz;
      }
    }
  }
  // Now get interface velocity
  // vx
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pco->CenterWidth1(k,j,is-1,ie+1,dxw1_);
      for (int i=is; i<=ie+1; ++i) {
        Real tau = 0.0;
        for (int ifr=0; ifr<nfreq; ++ifr) {
          Real sigmal = prad->sigma_a(k,j,i-1,ifr) + prad->sigma_s(k,j,i-1,ifr);
          sigmal *= taufact(k,j,i-1);
          Real sigmar = prad->sigma_a(k,j,i,ifr) + prad->sigma_s(k,j,i,ifr);
          sigmar *= taufact(k,j,i);
          tau += (dxw1_(i-1) * sigmal + dxw1_(i) * sigmar);
        }

        Real factor = 0.0;
        GetTaufactorAdv(tau,factor);
        Real vl = vel_source_(k,j,i-1,0);
        Real vr = vel_source_(k,j,i,0);
        adv_vel(0,k,j,i) = factor*(vl + (pco->x1f(i) - pco->x1v(i-1)) *
                                   (vr - vl)/(pco->x1v(i) - pco->x1v(i-1)));
      }
    }
  }
  if (je > js) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        Real ratio = (pco->x2f(j) - pco->x2v(j-1))/
                     (pco->x2v(j) - pco->x2v(j-1));
        pco->CenterWidth2(k,j-1,is,ie,dxw1_);
        pco->CenterWidth2(k,j,is,ie,dxw2_);
        for (int i=is; i<=ie; ++i) {
          Real tau = 0.0;
          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real sigmal = prad->sigma_a(k,j-1,i,ifr) + prad->sigma_s(k,j-1,i,ifr);
            sigmal *= taufact(k,j-1,i);
            Real sigmar = prad->sigma_a(k,j,i,ifr) + prad->sigma_s(k,j,i,ifr);
            sigmar *= taufact(k,j,i);
            tau += (dxw1_(i) * sigmal + dxw2_(i) * sigmar);
          }

          Real factor = 0.0;
          GetTaufactorAdv(tau,factor);

          Real vl = vel_source_(k,j-1,i,1);
          Real vr = vel_source_(k,j,i,1);
          adv_vel(1,k,j,i) = factor*(vl +  ratio * (vr - vl));
        }
      }
    }
  }

  if (ke > ks) {
    for (int k=ks; k<=ke+1; ++k) {
      Real ratio = (pco->x3f(k) - pco->x3v(k-1))/
                   (pco->x3v(k) - pco->x3v(k-1));
      for (int j=js; j<=je; ++j) {
        pco->CenterWidth3(k-1,j,is,ie,dxw1_);
        pco->CenterWidth3(k,j,is,ie,dxw2_);
        for (int i=is; i<=ie; ++i) {
          Real tau = 0.0;
          for (int ifr=0; ifr<nfreq; ++ifr) {
            Real sigmal = prad->sigma_a(k-1,j,i,ifr) + prad->sigma_s(k-1,j,i,ifr);
            sigmal *= taufact(k-1,j,i);
            Real sigmar = prad->sigma_a(k,j,i,ifr) + prad->sigma_s(k,j,i,ifr);
            sigmar *= taufact(k,j,i);
            tau += (dxw1_(i) * sigmal + dxw2_(i) * sigmar);
          }

          Real factor = 0.0;
          GetTaufactorAdv(tau,factor);

          Real vl = vel_source_(k-1,j,i,2);
          Real vr = vel_source_(k,j,i,2);
          adv_vel(2,k,j,i) = factor * (vl +  ratio * (vr - vl));
        }
      }
    }
  }
}


// f_l is the factor in the opposite direction of vel
// f_r is the factor in the same direction of vel
void RadIntegrator::SignalSpeed(const Real adv, const Real f_l,
                                const Real f_r, Real *vel, Real *smax, Real *smin) {
  for (int n=0; n<pmy_rad->nang; ++n) {
    if (vel[n] > 0.0) {
      smax[n] = f_r * vel[n];
      smin[n] = -f_l * vel[n];
    } else {
      smax[n] = -f_l * vel[n];
      smin[n] = f_r * vel[n];
    }
  }
  if (adv_flag_ == 0) {
    for (int n=0; n<pmy_rad->nang; ++n) {
      smax[n] += std::abs(adv);
      smin[n] -= std::abs(adv);
    }
  }
}

void RadIntegrator::SplitVelocity(Real *vel_l, Real *vel_r, const Real advl,
                                  const Real advr, Real *smax_l, Real *smin_l,
                                  Real *smax_r, Real *smin_r) {
  int tot_ang = pmy_rad->n_fre_ang;
  int iteration = pmy_rad->pmy_block->pmy_mesh->pimrad->ite_scheme;
  // the left side
  if (iteration == 0) {
    for (int n=0; n<tot_ang; ++n) {
      Real vl = vel_l[n] - advl;
      vel_ex_l_(n) = 0.0;
      vel_im_l_(n) = smin_l[n] * (smax_l[n] - vl);
    }
    for (int n=0; n<tot_ang; ++n) {
      Real vr = vel_r[n] - advr;
      vel_ex_r_(n) = 0.0;
      vel_im_r_(n) = smax_r[n] * (vr - smin_r[n]);
    }
  } else if (iteration == 1) {
    for (int n=0; n<tot_ang; ++n) {
      Real vl = vel_l[n] - advl;
      if (vl > 0.0) {
        vel_ex_l_(n) = -smin_l[n] * vl;
        vel_im_l_(n) = smin_l[n] * smax_l[n];
      } else {
        vel_ex_l_(n) = 0.0;
        vel_im_l_(n) = smin_l[n] * (smax_l[n] - vl);
      }
    }// end n

    for (int n=0; n<tot_ang; ++n) {
      Real vr = vel_r[n] - advr;
      if (vr > 0.0) {
        vel_ex_r_(n) = 0.0;
        vel_im_r_(n) = smax_r[n] * (vr - smin_r[n]);
      } else {
        vel_ex_r_(n) = smax_r[n] * vr;
        vel_im_r_(n) = -smax_r[n] * smin_r[n];
      }
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [SplitVelocity]"
        << std::endl << "ite_scheme '" << iteration << "' not allowed!";
    ATHENA_ERROR(msg);
  }
}

void RadIntegrator::GetTaufactor(const Real tau, Real &factor1, int dir) {
  std::stringstream msg;

  if (dir > 0) {
    if (tau_flag_ == 1) {
      Real tausq = tau * tau;
      factor1 = 1.0 - 0.5 * tausq;
      if (tausq > tau_asymptotic_lim) {
        factor1  = (1.0 - std::exp(-tausq))/tausq;
      }
      factor1 = std::sqrt(factor1);
    } else if (tau_flag_ == 2) {
      Real tausq = tau;
      if (tausq > 1)
        factor1 = 1.0/tausq;
      else
        factor1 = 1.0;

    } else {
      msg << "### FATAL ERROR in function [GetTaufactor]"
          << std::endl << "tau_flag_ '" << tau_flag_ << "' not allowed!";
      ATHENA_ERROR(msg);
    }
  } else {
    if (tau_flag_ == 1) {
      Real tausq = tau * tau;
      factor1 = tausq;
      if (tausq > tau_asymptotic_lim) {
        factor1  = (1.0 - std::exp(-tausq*tausq))/tausq;
      }
      factor1 = std::sqrt(factor1);
    } else if (tau_flag_ == 2) {
      Real tausq = tau;
      if (tausq > 1)
        factor1 = 1.0/tausq;
      else
        factor1 = tausq;
    } else {
      msg << "### FATAL ERROR in function [GetTaufactor]"
          << std::endl << "tau_flag_ '" << tau_flag_ << "' not allowed!";
      ATHENA_ERROR(msg);
    }
  }
}


// The tau factor for advection velocity
void RadIntegrator::GetTaufactorAdv(const Real tau, Real &factor) {
  Real tausq = tau * tau;
  factor = tausq - 0.5 * tausq * tausq;
  if (tausq > tau_asymptotic_lim)
    factor = (1.0 - std::exp(-tausq));
  return;
}


void RadIntegrator::PredictVel(AthenaArray<Real> &ir, int k, int j, int i,
                               Real dt, Real rho, Real *vx, Real *vy, Real *vz) {
  NRRadiation *prad = pmy_rad;

  const Real &prat = prad->prat;
  Real invcrat = 1.0/prad->crat;
  Real ct = dt * prad->reduced_c;
  const int& nang =prad->nang;
  const int& nfreq=prad->nfreq;
  // first, calculate the moments
  Real er =0.0, fr1=0.0, fr2=0.0, fr3=0.0,
      pr11=0.0,pr12=0.0,pr13=0.0,pr22=0.0,
      pr23=0.0,pr33=0.0;
  Real *weight = &(prad->wmu(0));
  for (int ifr=0; ifr<nfreq; ++ifr) {
    Real er_f = 0.0,fr1_f=0.0,fr2_f=0.0,fr3_f=0.0,
        pr11_f=0.0,pr12_f=0.0,pr13_f=0.0,pr22_f=0.0,
        pr23_f=0.0,pr33_f=0.0;

    Real *irn = &(ir(k,j,i,ifr*nang));
    Real *cosx = &(prad->mu(0,k,j,i,0));
    Real *cosy = &(prad->mu(1,k,j,i,0));
    Real *cosz = &(prad->mu(2,k,j,i,0));
    for (int n=0; n<nang; ++n) {
      Real irweight = weight[n] * irn[n];
      er_f   += irweight;
      fr1_f  += irweight * cosx[n];
      fr2_f  += irweight * cosy[n];
      fr3_f  += irweight * cosz[n];
      pr11_f += irweight * cosx[n] * cosx[n];
      pr12_f += irweight * cosx[n] * cosy[n];
      pr13_f += irweight * cosx[n] * cosz[n];
      pr22_f += irweight * cosy[n] * cosy[n];
      pr23_f += irweight * cosy[n] * cosz[n];
      pr33_f += irweight * cosz[n] * cosz[n];
    }


    er   += er_f;
    fr1  += fr1_f;
    fr2  += fr2_f;
    fr3  += fr3_f;
    pr11 += pr11_f;
    pr12 += pr12_f;
    pr13 += pr13_f;
    pr22 += pr22_f;
    pr23 += pr23_f;
    pr33 += pr33_f;
  }
  // calculate the frequency integrated opacity
  Real grey_sigma_s = 0.0;
  Real grey_sigma_a = 0.0;
  for (int ifr=0; ifr<nfreq; ++ifr) {
    grey_sigma_s += prad->sigma_s(k,j,i,ifr);
    grey_sigma_a += prad->sigma_a(k,j,i,ifr);
  }
  Real dtcsigma = ct * (grey_sigma_s + grey_sigma_a);

  Real vx0 = (*vx);
  Real vy0 = (*vy);
  Real vz0 = (*vz);

  Real m0x = prat * fr1 * invcrat + rho * vx0;
  Real m0y = prat * fr2 * invcrat + rho * vy0;
  Real m0z = prat * fr3 * invcrat + rho * vz0;


  Real vx11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr11)
              * invcrat * invcrat;
  Real vy11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr22)
              * invcrat * invcrat;
  Real vz11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr33)
              * invcrat * invcrat;
  Real vx12 = dtcsigma * prat * pr12 * invcrat * invcrat;
  Real vx13 = dtcsigma * prat * pr13 * invcrat * invcrat;
  Real vy12 = dtcsigma * prat * pr23 * invcrat * invcrat;
  Real rhs1 = rho * vx0 + dtcsigma * m0x;
  Real rhs2 = rho * vy0 + dtcsigma * m0y;
  Real rhs3 = rho * vz0 + dtcsigma * m0z;

  Real factor = vx11 * vy11 * vz11 - vy11 * vx13 * vx13 + 2.0 * vx12 * vx13 * vy12
                - vx11 * vy12 * vy12 - vx12 * vx12 * vz11;
  factor = 1.0/factor;

  (*vx) = factor*(rhs3*(vx12*vy12 - vx13*vy11) + rhs2*(
      vy12*vx13 - vx12*vz11) + rhs1*(vy11*vz11 - vy12*vy12));

  (*vy) = factor*(rhs3*(vx12*vx13 - vx11*vy12) + rhs2*(
      vx11*vz11 - vx13*vx13) + rhs1*(vx13*vy12 - vx12*vz11));

  (*vz) = factor*(rhs3*(vx11*vy11 - vx12*vx12) + rhs2*(
      vx12*vx13 - vx11*vy12) + rhs1*(vx12*vy12 - vx13*vy11));
}
