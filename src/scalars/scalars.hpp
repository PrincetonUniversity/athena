#ifndef SCALARS_SCALARS_HPP_
#define SCALARS_SCALARS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file scalars.hpp
//  \brief definitions for PassiveScalars class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;

//! \class PassiveScalars
//  \brief

// TODO(felker): consider renaming to Scalars
class PassiveScalars {
 public:
  // TODO(felker): pin is currently only used for checking ssprk5_4, otherwise unused.
  // Leaving as ctor parameter in case of run-time "nscalars" option
  PassiveScalars(MeshBlock *pmb, ParameterInput *pin);

  // public data:
  // mass fraction of each species
  AthenaArray<Real> s, s1, s2;  // (no more than MAX_NREGISTER allowed)
  AthenaArray<Real> s_flux[3];  // face-averaged flux vector

  // fourth-order intermediate quantities
  AthenaArray<Real> s_cc;      // cell-centered approximations

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_s_;

  CellCenteredBoundaryVariable sbvar;

  // public functions:
  // KGF: use inheritance for these functions / overall class?
  void AddFluxDivergence(const Real wght, AthenaArray<Real> &s_out);
  void CalculateFluxes(AthenaArray<Real> &s, const int order);

 private:
  MeshBlock* pmy_block;
  // scratch space used to compute fluxes
  AthenaArray<Real> sl_, sr_, slb_;
  AthenaArray<Real> dxw_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> dflx_;

  // fourth-order
  // 4D scratch arrays
  AthenaArray<Real> scr1_nkji_, scr2_nkji_;
  AthenaArray<Real> sl3d_, sr3d_;
  AthenaArray<Real> laplacian_l_fc_, laplacian_r_fc_;

  void ComputeUpwindFlux(const int k, const int j, const int il,
                         const int iu, // CoordinateDirection dir,
                         AthenaArray<Real> &sl, AthenaArray<Real> &sr,
                         AthenaArray<Real> &mass_flx,
                         AthenaArray<Real> &flx_out);
};
#endif // SCALARS_SCALARS_HPP_
