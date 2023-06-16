#ifndef CR_INTEGRATORS_CR_INTEGRATORS_HPP_
#define CR_INTEGRATORS_CR_INTEGRATORS_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../task_list/task_list.hpp"
#include "../cr.hpp"

class MeshBlock;
class ParameterInput;
class CosmicRay;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer
class CRIntegrator {
  friend class CosmicRay;
 public:
  CRIntegrator(CosmicRay *pcr, ParameterInput *pin);
  ~CRIntegrator();

  CosmicRay *pmy_cr;

  void FluxDivergence(const Real wght, AthenaArray<Real> &cr_out);
  void CalculateFluxes(AthenaArray<Real> &w,
          AthenaArray<Real> &bcc, AthenaArray<Real> &cr, const int order);

  void CRFlux(int fdir, int il, int iu,
              AthenaArray<Real> &w_l, AthenaArray<Real> &w_r,
              AthenaArray<Real> &vdiff_l, AthenaArray<Real> &vdiff_r,
              AthenaArray<Real> &flx);
  void AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc, AthenaArray<Real> &ucr);
  int cr_xorder;

 private:
  AthenaArray<Real> new_sol_;
  AthenaArray<Real> ucr_l_, ucr_r_, ucr_lb_; // for reconstruction
  AthenaArray<Real> vdiff_l_, vdiff_r_;

  AthenaArray<Real> grad_pc_, ec_source_, coord_source_;
  AthenaArray<Real> ucr_vel_; //array for reconstruction

  // temporary array to store the flux
  Real taufact_;
  int vel_flx_flag_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_, dflx_, cwidth2_, cwidth3_;
};

#endif // CR_INTEGRATORS_CR_INTEGRATORS_HPP_
