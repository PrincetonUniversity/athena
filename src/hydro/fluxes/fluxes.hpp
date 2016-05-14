#ifndef FLUXES_HPP
#define FLUXES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file fluxes.hpp
//  \brief defines class HydroFluxes, data and functions for hydro/MHD fluxes
//======================================================================================

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;
class MeshBlock;

//! \class HydroFluxes
//  \brief member functions implement various flux functions for hydro/MHD

class HydroFluxes {
public:
  HydroFluxes(Hydro *phydro, ParameterInput *pin);
  ~HydroFluxes();

  Hydro *pmy_hydro;  // ptr to Hydro containing this HydroIntegrator

  void CalculateFluxes(MeshBlock *pmb, AthenaArray<Real> &u, AthenaArray<Real> &w,
    FaceField &b, AthenaArray<Real> &bcc, const int step);

  void FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &u, AthenaArray<Real> &w,
    FaceField &b, AthenaArray<Real> &bcc, const int step);

  void RiemannSolver(const int k, const int j, const int il, const int iu,
    const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx);

private:
  // scratch space used in integrator
  AthenaArray<Real> wl_, wr_, flx_; 
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> bb_normal_;    // normal magnetic field, for (SR/GR)MHD
  AthenaArray<Real> lambdas_p_l_;  // most positive wavespeeds in left state
  AthenaArray<Real> lambdas_m_l_;  // most negative wavespeeds in left state
  AthenaArray<Real> lambdas_p_r_;  // most positive wavespeeds in right state
  AthenaArray<Real> lambdas_m_r_;  // most negative wavespeeds in right state
  AthenaArray<Real> g_, gi_;       // metric and inverse, for some GR Riemann solvers
  AthenaArray<Real> cons_;         // conserved state, for some GR Riemann solvers
};
#endif // HYDRO_INTEGRATOR_HPP
