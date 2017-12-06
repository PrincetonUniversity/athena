#ifndef FIELD_DIFFUSION_HPP
#define FIELD_DIFFUSION_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field_diffusion.hpp
//  \brief defines class FieldDiffusion
//  Contains data and functions that implement the diffusion processes

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Field;
class ParameterInput;
class Coordinates;
//! \class HydroDiffusion
//  \brief data and functions for physical diffusion processes in the hydro

class FieldDiffusion {
public:
  FieldDiffusion(Field *pfld, ParameterInput *pin);
  ~FieldDiffusion();

  // accessors

  // data
  bool field_diffusion_defined;
  EdgeField emf; //stress tensor

  // functions
  void CalcFieldDiffusionEMF(const FaceField &bi, const AthenaArray<Real> &bc, EdgeField &e);
  void AddFieldDiffusionEMF(EdgeField &e);
  void AddEnergyFlux(const AthenaArray<Real> &bc, AthenaArray<Real> *flx);//add resistive dissipation to energy
  // resistivity
  void Resistivity(const FaceField &bi,const AthenaArray<Real> &bc, EdgeField &e);
  Real NewDtFldDiff(Real len, int k, int j, int i);
  Real Eta_Ohmic(int k, int j, int i);
  void CurlB(const FaceField &bi, const AthenaArray<Real> &bc, EdgeField &crnt);
  //Real cnuiso2(int k, int j, int i);
  //void Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv);
  //void FaceXdx(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceXdy(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceXdz(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceYdx(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceYdy(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceYdz(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceZdx(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceZdy(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  //void FaceZdz(const int k, const int j, const int il, const int iu,
  //  const AthenaArray<Real> &prim, AthenaArray<Real> &len);

private:
  MeshBlock *pmb_;    // ptr to meshblock containing this HydroDiffusion
  Field *pmy_field_;  // ptr to Field containing this FieldDiffusion
  //Hydro *pmy_hydro_;  // ptr to Hydro containing this FieldDiffusion
  Coordinates *pco_;  // ptr to coordinates class
  Real ieta_, etaO_;  // Ohmic diffusion coeff
  EdgeField j_;       // curl of B at cell-center
  //AthenaArray<Real> x1area_,x2area_,x2area_p1_,x3area_,x3area_p1_;
  //AthenaArray<Real> vol_;
  //AthenaArray<Real> fx_,fy_,fz_;
  //AthenaArray<Real> dflx_,flx_;
};
#endif
