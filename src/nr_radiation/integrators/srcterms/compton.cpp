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
//! \file compton.cpp
//  \brief  Add compton source terms
//======================================================================================

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
//! \fn RadIntegrator::Compton()
//  \brief

// wmu_cm is the weight in the co-moving frame
// wmu_cm=wmu * 1/(1-vdotn/Crat)^2 / Lorz^2
// tran_coef is (1-vdotn/Crat)*Lorz
// rho is gas density
// tgas is gas temperature
// This function only update the absorption term in each cell

void RadIntegrator::Compton(AthenaArray<Real> &wmu_cm,
                            AthenaArray<Real> &tran_coef, Real *sigma_s,
                            Real dt, Real lorz, Real rho, Real &tgas,
                            AthenaArray<Real> &ir_cm) {
  const Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->crat;
  Real redfactor=pmy_rad->reduced_c/pmy_rad->crat;

  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->peos->GetGamma();
  Real telectron = 1.0/pmy_rad->telectron;

  // Polynomial coefficients for Compton
  Real coef[2];
  for (int i=0; i<2; ++i)
    coef[i] = 0.0;

  Real tgasnew = tgas;

  for (int ifr=0; ifr<nfreq; ++ifr) {
    Real dtcsigma = ct * sigma_s[ifr];
    Real rdtcsigma = redfactor * dtcsigma;

    // Add Simple Compton scattering using the partically updated jr and ir
    if (rdtcsigma > TINY_NUMBER) {
     // Calculate the sum \int gamma (1-vdotn/c) dw_0 4 dt csigma_s /T_e
      Real suma1 = 0.0, suma2 = 0.0, jr_cm=0.0, source=0.0;
      Real *irn = &(ir_cm(nang*ifr));
      Real *wmun = &(wmu_cm(0));
      Real *tcoef = &(tran_coef(0));
      for (int n=0; n<nang; n++) {
         jr_cm += irn[n] * wmun[n];
         suma1 += tcoef[n] * wmun[n] * 4.0 * rdtcsigma * telectron;
      }
      suma2 = 4.0 * lorz * prat * dtcsigma*(gamma-1.0)*telectron/rho;

      Real tr = std::sqrt(std::sqrt(jr_cm));
      Real trnew;

      if (std::abs(tr - tgas) > 1.e-12) {
        coef[1] = (1.0 + suma2* jr_cm)/(suma1 * jr_cm);
        coef[0] = -(1.0+suma2*jr_cm)/suma1-tgas;
        int flag = FouthPolyRoot(coef[1], coef[0], trnew);
        if (flag == -1) {
          trnew = tr;
          tgasnew = tgas;
          source = 0.0;
        } else {
          Real jrnew = trnew * trnew * trnew * trnew;
          tgasnew = (jrnew - jr_cm)/(suma1*jr_cm) + trnew;
          source = rdtcsigma * 4.0 * jr_cm * telectron * (tgasnew - trnew);
        }
      }

      // Update the co-moving frame specific intensity
      for (int n=0; n<nang; n++) {
        irn[n] += source * tcoef[n];
      }
    }
  }

  // Update gas temperature
  tgas = tgasnew;
  return;
}


void RadIntegrator::MultiGroupCompton(
    AthenaArray<Real> &wmu_cm, AthenaArray<Real> &tran_coef, Real dt, Real lorz, Real rho,
    Real &tgas_ini, Real &tgas, AthenaArray<Real> &ir_cm) {
  int bd_cell_flag = 0;

  const Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->crat;
  Real redfactor=pmy_rad->reduced_c/pmy_rad->crat;

  const int& nang=pmy_rad->nang;
  const int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->peos->GetGamma();
  Real one_telectron = 1.0/pmy_rad->telectron;
  //correction factor for high order terms of energy change per scattering
  // f_theta is 1 when kT/m_ec^2 -->0
  Real f_theta=ComptCorrection(tgas);

  //Compton scattering coefficient always use the frequency
  // independent kappa_es
  // rho * kappa_es=sigma_es

  //--------------------------------------------------------------------
  // Estimate gas temperature based on integrated Kompaneets Equation

  Real compt_coef = f_theta*lorz*ct*redfactor*rho*pmy_rad->kappa_es
                    *one_telectron;
  Real coef_a1=4*compt_coef;
  Real coef_a2=(gamma -1.0)*prat*coef_a1/(redfactor*rho);

  // now calculate integral of J
  Real *j_nu = &(sum_nu1_(0));
  Real *wmun = &(wmu_cm(0));
  Real *nu_grid = &(pmy_rad->nu_grid(0));
  Real *nu_cen = &(pmy_rad->nu_cen(0));
  Real *delta_nu = &(pmy_rad->delta_nu(0));

  //--------------------------------------------------------------------
  // do not include the correction factor when estimate T
  Real tgas_new = tgas;

  for (int ifr=0; ifr<nfreq; ++ifr) {
    j_nu[ifr] = 0.0;
    Real *irn = &(ir_cm(nang*ifr));
    for (int n=0; n<nang; ++n) {
      j_nu[ifr] += wmun[n] * irn[n];
    }
  }
  Real sum_jnu = 0.0;
  for (int ifr=0; ifr<nfreq; ++ifr) {
    sum_jnu += j_nu[ifr];
  }
  if (compton_t_ > 0) {
    // calculate integral \int nu J_nu d\nu
    Real sum_nu_jnu = 0.0;
    for (int ifr=0; ifr<nfreq-1; ++ifr)
      sum_nu_jnu += nu_cen[ifr] * j_nu[ifr];

  //The last frequency bin is special
  // determine the effective blackbody temperature
  //from frequency nu_grid[0], nu_grid[nfreq-1]
  // This is nu/Tr for the last bin

//  Real nu_tr = pmy_rad->EffectiveBlackBody(j_nu[nfreq-1], nu_grid[nfreq-1]);
//  Real eff_tr = nu_grid[nfreq-1]/nu_tr;

    Real nu_jnu_last_bin = 0.0;
    Real jnu_sq_last_bin = 0.0;


    pmy_rad->ConvertBBJWien(j_nu[nfreq-1], nu_grid[nfreq-1], tgas,
                         nu_jnu_last_bin, jnu_sq_last_bin);

  // convert from j_nu to \int nu J_nu d\nu
//  Real nu_jnu_last_bin = pmy_rad->BBJtoJnu(j_nu[nfreq-1],nu_grid[nfreq-1]);
//  Real jnu_sq_last_bin = pmy_rad->BBJToJONuSq(j_nu[nfreq-1],nu_grid[nfreq-1]);

    sum_nu_jnu += nu_jnu_last_bin;

  // calculate the integral \int (J_nu/nu)^2 d\nu
    Real sum_jonu_sq = 0.0;
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      sum_jonu_sq += j_nu[ifr] * j_nu[ifr]/(nu_cen[ifr] * nu_cen[ifr]
                     * delta_nu[ifr]);
    }
  // the last bin
    sum_jonu_sq += jnu_sq_last_bin;

  // construct the coefficients
    Real coef_r4 = sum_jnu + coef_a2 * sum_jnu * sum_jnu;
    Real coef_r5 = 0.25*coef_a1*sum_nu_jnu+coef_a1*PI_FOUR_POWER*sum_jonu_sq/60.0;
    Real coef_rhs = sum_jnu + coef_a1 * sum_jnu * tgas
                  + coef_a2 * sum_jnu * sum_jnu;

    Real r_ini = 1.0;
    Real r_four_sol = coef_rhs/(coef_r4 + coef_r5 * r_ini);
    r_four_sol = std::pow(r_four_sol,0.25);
    int count=0;
    while((count < 10) && (std::abs(r_ini-r_four_sol)/r_four_sol > 1.e-6)) {
      r_ini = r_four_sol;
      r_four_sol = coef_rhs/(coef_r4 + coef_r5 * r_ini);
      r_four_sol = std::pow(r_four_sol,0.25);
      count++;
    }

     tgas_new = tgas + ((gamma-1.0)*prat/(redfactor*rho)) * sum_jnu *
                 (1.0 - r_four_sol * r_four_sol * r_four_sol * r_four_sol);
  }
  //--------------------------------------------------------------------
  Real n_nusq_last_bin = 0.0;
  Real nf_last = 0.0;
// update
  pmy_rad->ConvertBBJWien2(j_nu[nfreq-1], nu_grid[nfreq-1], tgas_new,
                          n_nusq_last_bin,  nf_last);


  //-------------------------------------------------------------------------
  // now solve the Kompaneets equation

  // frequency centers are: [0] [1] [2]....| [f-2] | [f-1]
  // frequency faces are: [0] [1] [2]... [f-2]   [f-1]
  // do not calculate flux for the last face

  // calculate the photon number density
  // n_nu = (pi^4/15) J_nu/nu^3
  Real *n_nu = &(sum_nu2_(0));
  for (int ifr=0; ifr<nfreq-1; ++ifr)
    n_nu[ifr] = (PI_FOUR_POWER/15.0) * j_nu[ifr]
               /(nu_cen[ifr]*nu_cen[ifr]*nu_cen[ifr]*delta_nu[ifr]);

  // get the conserved total number density n \nu^2 d\nu
  // if we assume blackbody spectrum, \int n \nu^2 d\nu
  // = T^3 \int x^2/(exp(x)-1) dx
//  Real n_nusq_last_bin = pmy_rad->ConvertBBJNNu2(j_nu[nfreq-1],nu_grid[nfreq-1]);

  Real sum_n_nusq = n_nusq_last_bin;
  for (int ifr=0; ifr<nfreq-1; ++ifr)
    sum_n_nusq += n_nu[ifr] * nu_cen[ifr] * nu_cen[ifr] * delta_nu[ifr];

  // determine the quasi-equilibrium solution
  Real eq_sol_c  = QuasiEqSol(tgas_new, sum_n_nusq);

  // now calculate the quasi-equilibrium solution at the frequency center
  // equalibrium sol: n_nu = 1/(c exp(nu/T)-1);

  Real *eq_sol = &(eq_sol_(0));
  for (int ifr=0; ifr<nfreq-1; ++ifr) {
    eq_sol[ifr] = 1.0/(eq_sol_c * std::exp(nu_cen[ifr]/tgas_new) - 1);
  }

  // now determine the delta_coefficient at the frequency interface
  Real *delta_coef = &(sum_nu3_(0));
  // use the weight based on frequency grid
  delta_coef[0] = 0.5;
  // for nearly 0 photon case
  // assuming linear slope
  if (eq_sol_c > 1.e16) {
    for (int ifr=1; ifr<nfreq-1; ++ifr) {
      delta_coef[ifr] = (nu_cen[ifr]-nu_grid[ifr])/(nu_cen[ifr]-nu_cen[ifr-1]);
    }
  } else {
    for (int ifr=1; ifr<nfreq-1; ++ifr) {
      Real dndnu = (eq_sol[ifr]-eq_sol[ifr-1])
                 /(nu_cen[ifr]-nu_cen[ifr-1]);
      Real eq_n_coef = 1-4.0*tgas_new*dndnu;
      Real eq_n_face = 0.5*(std::sqrt(eq_n_coef) - 1.0);
      delta_coef[ifr] = std::max(Real(0.0),(eq_n_face - eq_sol[ifr])/(
          eq_sol[ifr-1]-eq_sol[ifr]));
    }
  }
  delta_coef[nfreq-1] = 0.5;

  //-------------------------------------------------------------------------
  // construct the matrix coefficients
  Real *com_b_face_coef = &(com_b_face_coef_(0));
  Real *com_d_face_coef = &(com_d_face_coef_(0));
  Real *com_b_coef_l = &(com_b_coef_l_(0));
  Real *com_b_coef_r = &(com_b_coef_r_(0));
  Real *com_d_coef_l = &(com_d_coef_l_(0));
  Real *com_d_coef_r = &(com_d_coef_r_(0));

  // the conservative equation to solve is
  // dn_center/dt * \nu_center^2 d\nu
  // = (Flux_rface - Flux_lface)

  com_b_face_coef[0] = 0.0;
  com_d_face_coef[0] = 0.0;
  // the face coefficient: comp_t_coef * nu_f^4 * (T/((d\nu))+(1+n)(1-\delta))
  for (int ifr=1; ifr<nfreq-1; ++ifr) {
    Real a_coef = compt_coef*nu_grid[ifr]*nu_grid[ifr]*nu_grid[ifr]*nu_grid[ifr];
    Real n_face = (1-delta_coef[ifr]) * n_nu[ifr]+delta_coef[ifr]*n_nu[ifr-1];
    Real tdnu_face=tgas_new/(nu_cen[ifr]-nu_cen[ifr-1]);
    com_b_face_coef[ifr] = a_coef*(tdnu_face+(1.0+n_face)*(1.0-delta_coef[ifr]));
    com_d_face_coef[ifr] = a_coef*((1+n_face)*delta_coef[ifr]-tdnu_face);
  }

  // for the first bin, the left boundary nu_grid[0]=0
  // com_a_coef_l[0] is 0
  // split the left and right hand side coefficients
  for (int ifr=0; ifr<nfreq-2; ++ifr) {
    Real bottom = 1.0/(nu_cen[ifr]*nu_cen[ifr]*delta_nu[ifr]);
    com_b_coef_r[ifr]=com_b_face_coef[ifr+1]*bottom;
    com_d_coef_r[ifr]=com_d_face_coef[ifr+1]*bottom;
  }

  for (int ifr=0; ifr<nfreq-1; ++ifr) {
    Real bottom = 1.0/(nu_cen[ifr]*nu_cen[ifr]*delta_nu[ifr]);
    com_b_coef_l[ifr]=com_b_face_coef[ifr]*bottom;
    com_d_coef_l[ifr]=com_d_face_coef[ifr]*bottom;
  }

  // also calculate left hand side coefficients for nfreq-2

  // now we invert the matrix
  // -B_r * n_f+1 + (1 - D_r + B_l) * n_f + D_l * n_f-1 = n_f^old
  // now solve the matrix from the left boundary
  // relate the right variable with the left variable
  // using Gauss elimination

  // each variable at nf is related to nf+1 as
  // n_f =  nf_n0 * n_{f+1} + nf_rhs

  Real *nf_rhs = &(nf_rhs_(0));
  Real *nf_n0 = &(nf_n0_(0));

  nf_rhs[0] = n_nu[0]/(1.0-com_d_coef_r[0]+com_b_coef_l[0]);
  nf_n0[0] = com_b_coef_r[0]/(1.0-com_d_coef_r[0]+com_b_coef_l[0]);

  for (int ifr=1; ifr<nfreq-2; ++ifr) {
    //eliminate n_{f-1}, get a relation between n_f and n_{f+1}
    // n_{f-1} = nf_n0[f-1] * n_{f} + nf_rhs[f-1]
    Real ef = (1.0-com_d_coef_r[ifr]+com_b_coef_l[ifr]);
    Real new_ef = 1.0 + com_d_coef_l[ifr] * nf_n0[ifr-1]/ef;
    nf_rhs[ifr] = n_nu[ifr]/ef-com_d_coef_l[ifr] * nf_rhs[ifr-1]/ef;
    if (std::abs(new_ef) < TINY_NUMBER) {
      bd_cell_flag = 1;
      break;
    }
    nf_rhs[ifr] /= new_ef;
    nf_n0[ifr] = com_b_coef_r[ifr]/(ef*new_ef);
  }

  // for the flux at the last face, we calculate the flux explicitly
  //  Real nf_last=1.0/(\lambda exp(nu_tr)-1.0);
  Real a_coef_last = compt_coef*nu_grid[nfreq-1]*nu_grid[nfreq-1]
                     *nu_grid[nfreq-1]*nu_grid[nfreq-1]
                     /(nu_cen[nfreq-2]*nu_cen[nfreq-2]*delta_nu[nfreq-2]);

  Real flux_last_explicit = 0.0;
  Real flux_last_im_coef = 0.0;

  if (nf_last < n_nu[nfreq-2]) {
    flux_last_explicit=a_coef_last*(nf_last*(1.0+nf_last)
                                    +tgas_new*nf_last*2.0/delta_nu[nfreq-2]);
    flux_last_im_coef=a_coef_last*tgas_new*2.0/delta_nu[nfreq-2];
  }

  // We use the bin nfreq-2 to get the solution
  Real n0_coef = (1.0+com_b_coef_l[nfreq-2]+flux_last_im_coef)
                 +com_d_coef_l[nfreq-2]*nf_n0[nfreq-3];
  Real rhs_coef = n_nu[nfreq-2]+flux_last_explicit
                  -com_d_coef_l[nfreq-2]*nf_rhs[nfreq-3];

  if (std::abs(n0_coef) < TINY_NUMBER) {
    bd_cell_flag = 1;
    return;
  }

  // update new solution at nfreq-2
  n_nu[nfreq-2] = rhs_coef/n0_coef;
  for (int ifr=nfreq-3; ifr>=0; --ifr) {
    n_nu[ifr] = nf_n0[ifr] * n_nu[ifr+1] + nf_rhs[ifr];
  }

  Real flux_last = flux_last_explicit - flux_last_im_coef * n_nu[nfreq-2];
  // update n_nusq in the last bin
  n_nusq_last_bin -= flux_last *
                     (nu_cen[nfreq-2]*nu_cen[nfreq-2]*delta_nu[nfreq-2]);

  n_nusq_last_bin = std::max(n_nusq_last_bin,static_cast<Real>(TINY_NUMBER));

  //-------------------------------------------------------------------------
  // now go from update n_nu to new_j_nu
  if (bd_cell_flag == 0) {
    Real *new_j_nu = &(new_j_nu_(0));
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      new_j_nu[ifr] = 15.0*ONE_PI_FOUR_POWER*n_nu[ifr]*nu_cen[ifr]
                      *nu_cen[ifr]*nu_cen[ifr]*delta_nu[ifr];
    }
    new_j_nu[nfreq-1] = pmy_rad->InverseConvertBBJNNu2Wien(n_nusq_last_bin,
                                                           nu_grid[nfreq-1], tgas_new);
    // now calculate the total j
    Real sum_new_jnu = 0.0;
    sum_new_jnu = new_j_nu[nfreq-1];
    for (int ifr=0; ifr<nfreq-1; ++ifr)
      sum_new_jnu += new_j_nu[ifr];

    // now apply energy conservation

    // Update intensity in each frequency bin
    // rho * tgas/(gamma-1) + (prat/redfactor)* sum_jnu =
    // rho * tgas_new/(gamma-1) + (prat/redfactor) * sum_new_jnu

    // now update tgas_new via energy conservation
    Real tgas_test = (prat/redfactor)*(sum_jnu-sum_new_jnu)
                     *((gamma-1.0)/rho) + tgas_ini;
    if (tgas_test < TINY_NUMBER)
      tgas_test = tgas_new;
    tgas = tgas_test;
    // now update intensity isotropically
    Real *tcoef = &(tran_coef(0));
    for (int ifr=0; ifr<nfreq; ++ifr) {
      Real *irn = &(ir_cm(nang*ifr));
      for (int n=0; n<nang; n++) {
        irn[n] += tcoef[n]* (new_j_nu[ifr] - j_nu[ifr]);
        irn[n] = std::max(irn[n],static_cast<Real>(TINY_NUMBER));
      }
    }
  }
}

Real RadIntegrator::QuasiEqSol(Real &tgas, Real &tot_n) {
  Real eq_sol_c = 1.0;
  Real n_t3 = tot_n/(tgas*tgas*tgas);
  if ((n_t3 < 2.3739) && (n_t3 > 0.6932)) {
    eq_sol_c = 1.948 * std::pow(n_t3, -1.016) + 0.1907;
  } else if (n_t3 <= 0.6932) {
    eq_sol_c = (1.0 + std::sqrt(1.0 + 0.25 * n_t3))/n_t3;
  }
  return eq_sol_c;
}

Real RadIntegrator::ComptCorrection(Real &tgas) {
  // this is kT/m_ec^2
  //Real theta=tgas/pmy_rad->telectron;
  //  Real f_theta=(1.0+3.683*theta+4.0*theta*theta)/(1.0+theta);
  Real f_theta = 1.0;
  return f_theta;
}
