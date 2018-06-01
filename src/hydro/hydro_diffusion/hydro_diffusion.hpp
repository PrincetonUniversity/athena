#ifndef HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
#define HYDRO_HYDRO_DIFFUSION_HYDRO_DIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_diffusion.hpp
//  \brief defines class HydroDiffusion
//  Contains data and functions that implement the diffusion processes

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;
class Coordinates;
class HydroDiffusion;


void ConstViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bc, int is, int ie, int js, int je, int ks, int ke);

void  ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bc, int is, int ie, int js, int je, int ks, int ke);


enum {ISO=0, ANI=1};

//! \class HydroDiffusion
//  \brief data and functions for physical diffusion processes in the hydro

class HydroDiffusion {
public:
  HydroDiffusion(Hydro *phyd, ParameterInput *pin);
  ~HydroDiffusion();

  // data
  bool hydro_diffusion_defined;
  Real nu_iso, nu_aniso; // viscosity coeff
  AthenaArray<Real> visflx[3]; // viscous stress tensor
  AthenaArray<Real> nu; // viscosity array

  Real kappa_iso, kappa_aniso; // thermal conduction coeff
  AthenaArray<Real> cndflx[3]; // thermal stress tensor
  AthenaArray<Real> kappa; // conduction array

  // functions
  void CalcHydroDiffusionFlux(const AthenaArray<Real> &p, const AthenaArray<Real> &c,
                                    AthenaArray<Real> *flx);
  void AddHydroDiffusionFlux(AthenaArray<Real> *flx_src, AthenaArray<Real> *flx_des);
  void AddHydroDiffusionEnergyFlux(AthenaArray<Real> *flux_src,
                                   AthenaArray<Real> *flux_des);
  void ClearHydroFlux(AthenaArray<Real> *flx);
  void SetHydroDiffusivity(AthenaArray<Real> &w, AthenaArray<Real> &bc);
  void NewHydroDiffusionDt(Real &dt_vis, Real &dt_cnd);

  // viscosity
  void ViscousFlux_iso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                             AthenaArray<Real> *flx);
  void ViscousFlux_aniso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                               AthenaArray<Real> *flx);

  // thermal conduction
  void ThermalFlux_iso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                             AthenaArray<Real> *flx);
  void ThermalFlux_aniso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                               AthenaArray<Real> *flx);

private:
  MeshBlock *pmb_;    // ptr to meshblock containing this HydroDiffusion
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroDiffusion
  Coordinates *pco_;  // ptr to coordinates class
  AthenaArray<Real> divv_; // divergence of velocity
  AthenaArray<Real> x1area_,x2area_,x2area_p1_,x3area_,x3area_p1_;
  AthenaArray<Real> vol_;
  AthenaArray<Real> fx_,fy_,fz_;
  AthenaArray<Real> dx1_,dx2_,dx3_;
  AthenaArray<Real> nu_tot_,kappa_tot_;

  // functions pointer to calculate spatial dependent coefficients
  ViscosityCoeff_t CalcViscCoeff_;
  ConductionCoeff_t CalcCondCoeff_;

  // auxiliary functions to calculate viscous flux
  void Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv);
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
