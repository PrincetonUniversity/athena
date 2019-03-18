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
//  \brief hydro data and functions

class PassiveScalars {
  // KGF: friend classes? Hydro?
 public:
  PassiveScalars(MeshBlock *pmb, ParameterInput *pin);
  ~PassiveScalars();

  // public data:
  CellCenteredBoundaryVariable *psbval;

  // mass fraction of each species
  AthenaArray<Real> s;
  AthenaArray<Real> s_flux[3];  // face-averaged flux vector

  // fourth-order intermediate quantities
  AthenaArray<Real> s_cc;      // cell-centered approximations

  // public functions:

  // KGF: make the following Hydro member fn a free function:
  // void WeightedAveU(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
  //                   AthenaArray<Real> &u_in2, const Real wght[3]);

  // KGF: use inheritance for these functions / overall class?
  void AddFluxDivergenceToAverage(AthenaArray<Real> &w, const Real wght);
  // AthenaArray<Real> &w, AthenaArray<Real> &bcc,
  //   const Real wght, AthenaArray<Real> &u_out);
  void CalculateFluxes(AthenaArray<Real> &w, const int order);
      // AthenaArray<Real> &w, FaceField &b,
      // AthenaArray<Real> &bcc, const int order);

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
