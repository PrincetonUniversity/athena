#ifndef FLUID_INTEGRATOR_HPP
#define FLUID_INTEGRATOR_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file fluid_integrator.hpp
//  \brief defines class HydroIntegrator, data and functions to integrate fluid
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Forward declarations
class MeshBlock;
class Hydro;
class ParameterInput;

//! \class HydroIntegrator
//  \brief member functions implement various integration algorithms for the fluid

class HydroIntegrator {
public:
  HydroIntegrator(Hydro *pf, ParameterInput *pin);
  ~HydroIntegrator();

  Hydro *pmy_fluid;  // ptr to Hydro containing this HydroIntegrator

  void OneStep(MeshBlock *pmb, AthenaArray<Real> &u, AthenaArray<Real> &w,
    InterfaceField &b, AthenaArray<Real> &bcc, const int step);
//  void Correct(MeshBlock *pmb);

  void RiemannSolver(const int k, const int j, const int il, const int iu,
    const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx);

  void PiecewiseLinearX1(const int k, const int j,
    const int il, const int iu,
    const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX2(const int k, const int j,
    const int il, const int iu,
    const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX3(const int k, const int j,
    const int il, const int iu,
    const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX1(const int k, const int j,
    const int il, const int iu,
    const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX2(const int k, const int j,
    const int il, const int iu,
    const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX3(const int k, const int j,
    const int il, const int iu,
    const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

private:
  // scratch space used in integrator
  AthenaArray<Real> wl_, wr_, flx_, jflx_j_, kflx_k_; 
  AthenaArray<Real> face_area_, face_area_p1_, cell_volume_;
  AthenaArray<Real> bb_normal_;    // normal magnetic field, for (SR/GR)MHD
  AthenaArray<Real> lambdas_p_l_;  // most positive wavespeeds in left state
  AthenaArray<Real> lambdas_m_l_;  // most negative wavespeeds in left state
  AthenaArray<Real> lambdas_p_r_;  // most positive wavespeeds in right state
  AthenaArray<Real> lambdas_m_r_;  // most negative wavespeeds in right state
  AthenaArray<Real> g_, gi_;       // metric and inverse, for some GR Riemann solvers
  AthenaArray<Real> cons_;         // conserved state, for some GR Riemann solvers
};
#endif
