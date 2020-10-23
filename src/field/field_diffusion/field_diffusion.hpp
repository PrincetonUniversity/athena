#ifndef FIELD_FIELD_DIFFUSION_FIELD_DIFFUSION_HPP_
#define FIELD_FIELD_DIFFUSION_FIELD_DIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field_diffusion.hpp
//! \brief defines class FieldDiffusion

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Field;
class Hydro;
class ParameterInput;
class Coordinates;

class FieldDiffusion;

//! currently must be free function for compatibility with user-defined fn via fn pointers
void ConstDiffusivity(FieldDiffusion *pfdif, MeshBlock *pmb, const AthenaArray<Real> &w,
                      const AthenaArray<Real> &bmag,
                      const int is, const int ie, const int js, const int je,
                      const int ks, const int ke);

//! \class FieldDiffusion
//! \brief data and functions for physical diffusion processes in the field

class FieldDiffusion {
 public:
  FieldDiffusion(MeshBlock *pmb, ParameterInput *pin);

  // data
  MeshBlock* pmy_block;
  bool field_diffusion_defined;
  Real eta_ohm, eta_hall, eta_ad;
  AthenaArray<Real> etaB; // 4-dim array, covering O/H/A altogether
  EdgeField e_oa, e_h;     // edge-centered electric field from non-ideal MHD
  FaceField pflux;        // face-centered energy (Poynting) flux

  AthenaArray<Real> jfx, jfy, jfz; // interface current density (for HLL Riemann solver)
  AthenaArray<Real> jcc;     // cell-centered current density (for the integrator)

  // array indices for magnetic diffusion types
  // should not be scoped (C++11) since enumerators are only used as "int" to index arrays
  enum DiffProcess {ohmic=0, hall=1, ambipolar=2};
  // TODO(felker) Unlike HydroDiffusion::DiffProcess, not using optional "unscoped enum
  // name" qualifier when referencing the enumerators in other files. Be consistent

  // alternative to unscoped (possibly anonymous) for int constants:
  // static constexpr int n_ohmic = 0;
  // static constexpr int n_hall = 1;
  // static constexpr int n_ambi = 2;

  // functions
  void CalcDiffusionEMF(FaceField &bi, const AthenaArray<Real> &bc, EdgeField &e);
  void AddEMF(const EdgeField &e_src, EdgeField &e_des);
  void ClearEMF(EdgeField &e);
  void CalcCurrent(FaceField &b);
  void AddPoyntingFlux (FaceField &p_src);

  // functions to calculate diffusivities and timesteps
  void SetDiffusivity(const AthenaArray<Real> &w, const AthenaArray<Real> &bcc);
  void NewDiffusionDt(Real &dt_oa, Real &dt_h);

  // non-ideal MHD EMFs
  void OhmicEMF(const FaceField &b, const AthenaArray<Real> &bc, EdgeField &e);
  //void HallEMF(const FaceField &b, const AthenaArray<Real> &bc, EdgeField &e);
  void AmbipolarEMF(const FaceField &b, const AthenaArray<Real> &bc, EdgeField &e);

  // functions for energy flux
  void PoyntingFlux(EdgeField &e, const AthenaArray<Real> &bcc);

 private:
  AthenaArray<Real> bmag_; // B field strength
  EdgeField jedge_;       // curl of B
  //EdgeField eh1_, eh2_, eh3_; // scratch arrays for the Hall integrator

  FieldDiffusionCoeffFunc CalcMagDiffCoeff_; // calculate magnetic diffusivities

  AthenaArray<Real> face_area_, face_area_p1_, edge_length_, edge_length_m1_;
  AthenaArray<Real>  cell_volume_;
  AthenaArray<Real> dx1_, dx2_, dx3_, len_;
  AthenaArray<Real> eta_tot_;
};
#endif // FIELD_FIELD_DIFFUSION_FIELD_DIFFUSION_HPP_
