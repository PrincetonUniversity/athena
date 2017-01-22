#ifndef HYDRO_DIFFUSION_HPP
#define HYDRO_DIFFUSION_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file diffusion.hpp
//  \brief defines class HydroDiffusion
//  Contains data and functions that implement the diffusion processes

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;
class Coordinates;
//! \class HydroDiffusion
//  \brief data and functions for physical diffusion processes in the hydro

class HydroDiffusion {
public:
  HydroDiffusion(Hydro *phyd, ParameterInput *pin);
  ~HydroDiffusion();

  // accessors

  // data
  bool hydro_diffusion_defined;
  AthenaArray<Real> diflx[3]; //stress tensor

  // functions
  void CalcHydroDiffusionFlux(const AthenaArray<Real> &p, const AthenaArray<Real> &c, AthenaArray<Real> *flx);
  void AddHydroDiffusionFlux(AthenaArray<Real> *flx);
  void AddEnergyFlux(const AthenaArray<Real> &bc, AthenaArray<Real> *flx);
  // viscosity
  void Viscosity(const AthenaArray<Real> &p,const AthenaArray<Real> &c, AthenaArray<Real> *flx);
  Real NewDtDiff(Real len, int k, int j, int i);
  //Real nuiso1(int n, int k, int j, int i);
  //Real cnuiso2(int k, int j, int i);
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

private:
  MeshBlock *pmb_;    // ptr to meshblock containing this HydroDiffusion
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroDiffusion
  Coordinates *pco_;  // ptr to coordinates class
  Real nuiso_;        // iso viscosity
  AthenaArray<Real> divv_; // divergence of velocity
  AthenaArray<Real> x1area_,x2area_,x2area_p1_,x3area_,x3area_p1_;
  AthenaArray<Real> vol_;
  AthenaArray<Real> fx_,fy_,fz_;
  //AthenaArray<Real> dflx_,flx_;
};
#endif
