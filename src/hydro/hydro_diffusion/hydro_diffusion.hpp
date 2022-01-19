#ifndef HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
#define HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_diffusion.hpp
//! \brief defines class HydroDiffusion
//!
//! Contains data and functions that implement the diffusion processes

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;
class Coordinates;
class HydroDiffusion;

// currently must be free functions for compatibility with user-defined fn via fn pointers
void ConstViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &w,
                    const AthenaArray<Real> &bc,
                    int is, int ie, int js, int je, int ks, int ke);

void ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &w,
                     const AthenaArray<Real> &bc,
                     int is, int ie, int js, int je, int ks, int ke);

//! \class HydroDiffusion
//! \brief data and functions for physical diffusion processes in the hydro

class HydroDiffusion {
 public:
  HydroDiffusion(Hydro *phyd, ParameterInput *pin);

  // data
  bool hydro_diffusion_defined;
  Real nu_iso, nu_aniso; // viscosity coeff
  AthenaArray<Real> visflx[3]; // viscous stress tensor
  AthenaArray<Real> nu; // viscosity array

  Real kappa_iso, kappa_aniso; // thermal conduction coeff
  AthenaArray<Real> cndflx[3]; // thermal stress tensor
  AthenaArray<Real> kappa; // conduction array

  // array indices for hydro diffusion (conduction & viscosity) variants: directionality
  // should not be scoped (C++11) since enumerators are only used as "int" to index arrays
  enum DiffProcess {iso=0, aniso=1};

  // functions
  // iprim is primitives in the inertial system
  void CalcDiffusionFlux(const AthenaArray<Real> &prim, const AthenaArray<Real> &iprim,
                         const AthenaArray<Real> &bcc);
  // TODO(felker): Rename+move out of this class. Confusing w/ Hydro::AddDiffusionFluxes()
  // See note in hydro_diffusion.cpp.
  void AddDiffusionFlux(AthenaArray<Real> *flx_src, AthenaArray<Real> *flx_des);
  void AddDiffusionEnergyFlux(AthenaArray<Real> *flux_src, AthenaArray<Real> *flux_des);
  void ClearFlux(AthenaArray<Real> *flx);
  void SetDiffusivity(const AthenaArray<Real> &w, const AthenaArray<Real> &bc);
  void NewDiffusionDt(Real &dt_vis, Real &dt_cnd);

  // viscosity
  void ViscousFluxIso(const AthenaArray<Real> &p, const AthenaArray<Real> &p_i,
                      AthenaArray<Real> *flx);
  void ViscousFluxAniso(const AthenaArray<Real> &p, const AthenaArray<Real> &p_i,
                        AthenaArray<Real> *flx);

  // thermal conduction
  void ThermalFluxIso(const AthenaArray<Real> &p, AthenaArray<Real> *flx);
  void ThermalFluxAniso(const AthenaArray<Real> &p, AthenaArray<Real> *flx);

 private:
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroDiffusion
  MeshBlock *pmb_;    // ptr to meshblock containing this HydroDiffusion
  Coordinates *pco_;  // ptr to coordinates class
  AthenaArray<Real> div_vel_; // divergence of velocity
  AthenaArray<Real> x1area_, x2area_, x2area_p1_, x3area_, x3area_p1_;
  AthenaArray<Real> vol_;
  AthenaArray<Real> fx_, fy_, fz_;
  AthenaArray<Real> dx1_, dx2_, dx3_;
  AthenaArray<Real> nu_tot_, kappa_tot_;

  // functions pointer to calculate spatial dependent coefficients
  ViscosityCoeffFunc CalcViscCoeff_;
  ConductionCoeffFunc CalcCondCoeff_;

  // auxiliary functions to calculate viscous flux
  void DivVelocity(const AthenaArray<Real> &prim, AthenaArray<Real> &divv);
  void FaceXdx(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceXdy(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceXdz(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdx(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdy(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdz(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdx(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdy(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdz(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &prim, AthenaArray<Real> &len);
};
#endif // HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
