#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro.hpp
//  \brief definitions for Hydro class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
//#include "../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
class HydroSourceTerms;
class HydroDiffusion;
struct IntegratorWeight;

//! \class Hydro
//  \brief hydro data and functions

class Hydro {
friend class Field;
public:
  Hydro(MeshBlock *pmb, ParameterInput *pin);
  ~Hydro();

  // data
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Hydro
  // conserved and primitive variables
  AthenaArray<Real> u,w;      // time-integrator memory register #1
  AthenaArray<Real> u1,w1;    // time-integrator memory register #2
  AthenaArray<Real> u2;       // time-integrator memory register #3
  // (no more than MAX_NREGISTER allowed)

  AthenaArray<Real> flux[3];  // face-averaged flux vector

  HydroSourceTerms *psrc;
  HydroDiffusion *phdif;

  // functions
  Real NewBlockTimeStep(void);    // computes new timestep on a MeshBlock
  void WeightedAveU(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
    AthenaArray<Real> &u_in2, const Real wght[3]);
  void AddFluxDivergenceToAverage(AthenaArray<Real> &w, AthenaArray<Real> &bcc,
    const Real wght, AthenaArray<Real> &u_out);
  void CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
    AthenaArray<Real> &bcc, int order);
  void RiemannSolver(const int kl, const int ku, const int jl, const int ju,
    const int il, const int iu, const int ivx, const AthenaArray<Real> &bx,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
    AthenaArray<Real> &e1, AthenaArray<Real> &e2);

  void AddGravityFlux(void);
  void AddGravityFluxWithGflx(void);
  void CalculateGravityFlux(AthenaArray<Real> &phi_in);
  void CorrectGravityFlux(void);

private:
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep
  // scratch space used to compute fluxes
  AthenaArray<Real> wl_, wr_;
  AthenaArray<Real> dxw_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> dflx_;
  AthenaArray<Real> bb_normal_;    // normal magnetic field, for (SR/GR)MHD
  AthenaArray<Real> lambdas_p_l_;  // most positive wavespeeds in left state
  AthenaArray<Real> lambdas_m_l_;  // most negative wavespeeds in left state
  AthenaArray<Real> lambdas_p_r_;  // most positive wavespeeds in right state
  AthenaArray<Real> lambdas_m_r_;  // most negative wavespeeds in right state
  AthenaArray<Real> g_, gi_;       // metric and inverse, for some GR Riemann solvers
  AthenaArray<Real> cons_;         // conserved state, for some GR Riemann solvers

  // self-gravity
  AthenaArray<Real> gflx[3], gflx_old[3]; // gravity tensor (old Athena style)


  TimeStepFunc_t UserTimeStep_;
};
#endif // HYDRO_HYDRO_HPP_
