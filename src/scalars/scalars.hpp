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

class PassiveScalars {
 public:
  PassiveScalars(MeshBlock *pmb, ParameterInput *pin);

  // public data:
  // mass fraction of each species
  AthenaArray<Real> s;
  AthenaArray<Real> s_flux[3];  // face-averaged flux vector

  // fourth-order intermediate quantities
  AthenaArray<Real> s_cc;      // cell-centered approximations

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_s_;

  CellCenteredBoundaryVariable sbvar;

  // public functions:
  // KGF: use inheritance for these functions / overall class?
  void AddFluxDivergenceToAverage(AthenaArray<Real> &w, const Real wght);
  void CalculateFluxes(AthenaArray<Real> &w, const int order);

 private:
  MeshBlock* pmy_block;
  // scratch space used to compute fluxes
  AthenaArray<Real> wl_, wr_, wlb_;
  AthenaArray<Real> dxw_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;

  // fourth-order hydro
  // 4D scratch arrays
  AthenaArray<Real> scr1_nkji_;
  AthenaArray<Real> wl3d_, wr3d_;
  AthenaArray<Real> laplacian_l_fc_, laplacian_r_fc_;
};
#endif // SCALARS_SCALARS_HPP_
