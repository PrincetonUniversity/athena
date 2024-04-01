#ifndef NR_RADIATION_INTEGRATORS_RAD_INTEGRATORS_HPP_
#define NR_RADIATION_INTEGRATORS_RAD_INTEGRATORS_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file rad_integrators.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../../athena_arrays.hpp"
#include "../../task_list/task_list.hpp"
#include "../radiation.hpp"

class MeshBlock;
class ParameterInput;
class NRRadiation;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer
class RadIntegrator {
  friend class NRRadiation;
  friend class IMRadiation;
 public:
  RadIntegrator(NRRadiation *prad, ParameterInput *pin);

  NRRadiation *pmy_rad;

  void FluxDivergence(const Real wght, AthenaArray<Real> &ir_in,
                                       AthenaArray<Real> &ir_out);
  void FluxDivergence(const Real wght);


  void FirstOrderFluxDivergenceCoef(const Real wght);
  void FirstOrderFluxDivergence(const int k, const int j, const int i,
                                 AthenaArray<Real> &ir);


  void ImplicitAngularFluxesCoef(const Real wght);
  void ImplicitAngularFluxes(const int k, const int j, const int i,
                             AthenaArray<Real> &ir);

  void ImplicitPsiFluxCoef(int k, int j, int i, int n_zeta, Real wght,
            Real zeta_coefr, Real zeta_coefl);
  void ImplicitPsiFlux(int k, int j, int i, int ifr, int n_zeta, AthenaArray<Real> &ir);

  void CalculateFluxes(AthenaArray<Real> &w,
                       AthenaArray<Real> &ir, const int order);
  void CalculateFluxes(AthenaArray<Real> &ir, const int order);

  void CalSourceTerms(MeshBlock *pmb, const Real dt, const int k, const int j,
                      const int i, AthenaArray<Real> &u,
                      AthenaArray<Real> &ir_ini, AthenaArray<Real> &ir);

  void AddSourceTerms(MeshBlock *pmb, AthenaArray<Real> &u);

  void GetHydroSourceTerms(MeshBlock *pmb,
                       AthenaArray<Real> &ir_ini, AthenaArray<Real> &ir);


  Real AbsorptionScattering(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real *sigma_a, Real *sigma_p,
          Real *sigma_ae, Real *sigma_s, Real dt, Real lorz, Real rho, Real &tgas,
          AthenaArray<Real> &implicit_coef_, AthenaArray<Real> &ir_cm);

  void GetTgasVel(MeshBlock *pmb, const Real dt,
    AthenaArray<Real> &u, AthenaArray<Real> &w,
    AthenaArray<Real> &bcc, AthenaArray<Real> &ir);


  void Compton(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real *sigma_s,
          Real dt, Real lorz, Real rho, Real &tgas, AthenaArray<Real> &ir_cm);



  void LabToCom(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          Real *ir_lab, AthenaArray<Real> &ir_cm);

  void ComToLab(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          AthenaArray<Real> &ir_cm, Real *ir_lab);

  void ComAngle(const Real vx, const Real vy, const Real vz,
          Real mux, Real muy, Real muz, Real *mux0, Real *muy0, Real *muz0);

  //====================================
  // multi-group functions

  void MapLabToCmFrequency(Real &tran_coef, AthenaArray<Real> &split_ratio,
                   AthenaArray<int> &map_start, AthenaArray<int> &map_end,
                   AthenaArray<Real> &ir_cm, AthenaArray<Real> &ir_shift);
  void MapCmToLabFrequency(Real &tran_coef, AthenaArray<Real> &split_ratio,
                   AthenaArray<int> &map_start, AthenaArray<int> &map_end,
                      AthenaArray<Real> &ir_shift, AthenaArray<Real> &ir_cm);

  void GetCmMCIntensity(AthenaArray<Real> &ir_cm, AthenaArray<Real> &delta_nu_n,
                                                    AthenaArray<Real> &ir_face);

  void ForwardSplitting(Real &tran_coef,
                      AthenaArray<Real> &ir_cm, AthenaArray<Real> &ir_face,
                      AthenaArray<Real> &split_ratio,
                      AthenaArray<int> &map_start,AthenaArray<int> &map_end);

  void BackwardSplitting(Real &tran_coef,
                      AthenaArray<Real> &ir_cm, AthenaArray<Real> &ir_face,
                      AthenaArray<Real> &split_ratio,
                      AthenaArray<int> &map_start,AthenaArray<int> &map_end);


  void MapIrcmFrequency(AthenaArray<Real> &split_ratio,
                          AthenaArray<int> &map_start, AthenaArray<int> &map_end,
                  AthenaArray<Real> &input_array, AthenaArray<Real> &shift_array);

  bool FreMapMatrix(AthenaArray<Real> &split_ratio, Real &tran_coef,
    AthenaArray<int> &map_bin_start, AthenaArray<int> &map_bin_end,
    AthenaArray<int> &map_count, AthenaArray<Real> &map_matrix);


  bool InverseMapFrequency(Real &tran_coef,
         AthenaArray<int> &map_count, AthenaArray<Real> &map_matrix,
         AthenaArray<Real> &input_array, AthenaArray<Real> &shift_array);


  void SplitFrequencyBinLinear(int &l_bd, int &r_bd,
                  Real *nu_lab, Real &nu_l, Real &nu_r, Real &ir_l,
                                      Real &ir_r, Real *split_ratio);



  void ComToLabMultiGroup(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          AthenaArray<Real> &ir_cm, Real *ir_lab);

  Real MultiGroupAbsScat(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real *sigma_a, Real *sigma_p,
          Real *sigma_ae, Real *sigma_s, Real dt, Real lorz, Real rho, Real &tgas,
          AthenaArray<Real> &implicit_coef, AthenaArray<Real> &ir_cm);

  void AddMultiGroupCompt(MeshBlock *pmb, const Real dt,
        AthenaArray<Real> &u, AthenaArray<Real> &ir);

  // multigroup compton scattering function
  void MultiGroupCompton(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real dt, Real lorz, Real rho,
          Real &tgas_ini, Real &tgas, AthenaArray<Real> &ir_cm);

  Real QuasiEqSol(Real &tgas, Real &tot_n);
  Real ComptCorrection(Real &tgas);

  void GetTaufactor(const Real tau, Real &factor1, int dir);
  void GetTaufactorAdv(const Real tau, Real &factor);

  void PredictVel(AthenaArray<Real> &ir, int k, int j, int i, Real dt, Real rho,
                  Real *vx, Real *vy, Real *vz);

  void SignalSpeed(const Real adv, const Real factor1,
                 const Real factor2, Real *vel, Real *smax, Real *smin);

  void SplitVelocity(Real *vel_l, Real *vel_r, const Real advl,
            const Real advr, Real *smax_l, Real *smin_l, Real *smax_r, Real *smin_r);

  int rad_xorder;
  AthenaArray<Real> adv_vel; // the advectioin velocity that we separate
  AthenaArray<Real> taufact;
  AthenaArray<Real> rad_source; // store the radiation source terms
  AthenaArray<Real> delta_source; // calculate delta E, delta Fr per frequency
                                  // then sum the different across fre

 private:
  AthenaArray<Real> vel_, velx_,vely_,velz_;
  AthenaArray<Real> il_, ilb_, ir_;// for recontruction
                          // temporary array to store the flux, velocity
  AthenaArray<Real> vncsigma_, vncsigma2_, wmu_cm_, tran_coef_, ir_cm_;
  AthenaArray<Real> cm_to_lab_;
  AthenaArray<Real> g_zeta_, q_zeta_, ql_zeta_, qr_zeta_, zeta_flux_, zeta_area_;
  AthenaArray<Real> g_psi_, q_psi_, ql_psi_, qr_psi_, psi_flux_, psi_area_;
  AthenaArray<Real> dflx_ang_, ang_vol_;
  AthenaArray<Real> tgas_, vel_source_, tgas_new_; // array to store gas temperature,
                                        // velocity for source term
  // these are temporary arrays for multi-group source terms
  AthenaArray<Real> sum_nu3_, sum_nu2_, sum_nu1_;
  AthenaArray<Real> eq_sol_;
  // The coefficients to solve the Kompaneets equation
  AthenaArray<Real> com_b_face_coef_, com_d_face_coef_;
  AthenaArray<Real> com_b_coef_l_, com_b_coef_r_;
  AthenaArray<Real> com_d_coef_l_, com_d_coef_r_;
  AthenaArray<Real> nf_rhs_, nf_n0_, new_j_nu_;

  int tau_flag_;
  // flag to add simple Compton scattering
  // compton_t_ flag to estimate gas temperature or not for Compton scattering
  int compton_flag_, compton_t_, split_compton_;
  int adv_flag_; // flag used to indicate whether separate
                 // advection flux from diffustion flux or not.
  int imp_ang_flx_;

  int doppler_flag_;

  int flux_correct_flag_; // flag to do second order flux crrection or not.
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_, dflx_, cwidth2_, cwidth3_;
  AthenaArray<Real> adv_flx_;

  AthenaArray<Real> const_coef_, exp_coef_;
  AthenaArray<Real> const_coef1_l_, const_coef1_r_;
  AthenaArray<Real> const_coef2_l_, const_coef2_r_;
  AthenaArray<Real> const_coef3_l_, const_coef3_r_;
  AthenaArray<Real> divflx_, implicit_coef_, ang_flx_, imp_ang_coef_;
  AthenaArray<Real> imp_ang_coef_r_;
  AthenaArray<Real> imp_ang_psi_l_, imp_ang_psi_r_;
  AthenaArray<Real> left_coef1_, left_coef2_, left_coef3_;
  AthenaArray<Real> limiter_, limiterj_, limiterk_, dql_, dqr_;
  AthenaArray<Real> sfac1_x_, sfac2_x_;
  AthenaArray<Real> sfac1_y_, sfac2_y_;
  AthenaArray<Real> sfac1_z_, sfac2_z_;
  AthenaArray<Real> sm_diff1_, sm_diff2_;
  AthenaArray<Real> vel_ex_l_, vel_im_l_, vel_ex_r_, vel_im_r_;
  AthenaArray<Real> dxw1_, dxw2_;

  //----------------------------------------------------
  // array for multi-group
  // This is the coefficient in front I when calculate frequency flux

  // This is the actual flux in frequency space
  // shift amount from the frequency boundary
  AthenaArray<Real> split_ratio_, fre_map_matrix_;
  AthenaArray<Real> delta_nu_n_; // width of frequency bins for each angle
  AthenaArray<Real> ir_buff_;
  AthenaArray<Real> ir_face_, ir_ori_, ir_done_;
  AthenaArray<int> map_bin_start_, map_bin_end_;
  AthenaArray<int> map_count_;
  AthenaArray<Real> nu_shift_;
  int iteration_tgas_, iteration_compton_;
  Real tgas_error_, compton_error_;
  int nmax_map_; //maximum number of frequency bins that each bin will map to
};

#endif // NR_RADIATION_INTEGRATORS_RAD_INTEGRATORS_HPP_
