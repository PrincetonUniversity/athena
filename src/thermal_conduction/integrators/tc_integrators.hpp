#ifndef TCINTEGRATORS_HPP
#define TCINTEGRATORS_HPP
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

class Hydro;
class ParameterInput;
class ThermalConduction;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class TCIntegrator {
  friend class ThermalConduction;
public:
  TCIntegrator(ThermalConduction *ptc, ParameterInput *pin);
  ~TCIntegrator();
  
  ThermalConduction *pmy_tc;

  int tc_xorder;

  void AddSourceTerms(MeshBlock *pmb, const Real dt,  AthenaArray<Real> &u,
                 AthenaArray<Real> &w, AthenaArray<Real> &u_tc);

  void CalculateFluxes(AthenaArray<Real> &w,
            AthenaArray<Real> &bcc, AthenaArray<Real> &tc, const int order);


  void FluxDivergence(const Real wght, AthenaArray<Real> &w, 
                                       AthenaArray<Real> &tc_out);

  void TCFlux(int fdir, int il, int iu, 
      AthenaArray<Real> &w_l, AthenaArray<Real> &w_r,  
      AthenaArray<Real> &vdiff_l, AthenaArray<Real> &vdiff_r,      
                                      AthenaArray<Real> &flx); 




private:
  AthenaArray<Real> utc_l_, utc_r_, utc_lb_; // for reconstruction
  AthenaArray<Real> vdiff_, vdiff_l_, vdiff_r_;

  AthenaArray<Real> tc_esource_, coord_source_;
  AthenaArray<Real> utc_rho_t_; //array for reconstruction

    // temporary array to store the flux
  Real taufact_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_, dflx_, cwidth2_, cwidth3_;


};

#endif // TCINTEGRATORS_HPP
