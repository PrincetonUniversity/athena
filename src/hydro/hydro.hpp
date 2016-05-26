#ifndef HYDRO_HPP
#define HYDRO_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file hydro.hpp
//  \brief definitions for Hydro class
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class MeshBlock;
class ParameterInput;
class HydroFluxes;
class HydroSourceTerms;

//! \class Hydro
//  \brief hydro data and functions

class Hydro {
friend class Field;
public:
  Hydro(MeshBlock *pmb, ParameterInput *pin);
  ~Hydro();
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Hydro

  AthenaArray<Real> u,w;      // conserved and primitive variables
  AthenaArray<Real> u1,w1;    // conserved and primitive variables at intermediate step
  AthenaArray<Real> flux[3];  // conserved and primitive variables
  AthenaArray<Real> g, g_inv; // metric and its inverse
  AthenaArray<Real> ifov;     // internal hydro output variables for analysis

  HydroFluxes *pflux;
  HydroSourceTerms *psrc;

  Real NewBlockTimeStep(MeshBlock *pmb);    // computes new timestep on a MeshBlock
  void CopyOrAverageHydro(AthenaArray<Real> &a, AthenaArray<Real> &b, 
    AthenaArray<Real> &c, Real f);
  void FluxDivergence(MeshBlock *pmb,AthenaArray<Real> &u,
    AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &bcc, const int step);
  void CalculateFluxes(MeshBlock *pmb, AthenaArray<Real> &u, AthenaArray<Real> &w,
    FaceField &b, AthenaArray<Real> &bcc, const int step);
  void RiemannSolver(const int k, const int j, const int il, const int iu,
    const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx);

private:
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep
  // scratch space used to compute fluxes
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
#endif // HYDRO_HPP
