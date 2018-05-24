#ifndef RECONSTRUCT_RECONSTRUCTION_HPP_
#define RECONSTRUCT_RECONSTRUCTION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.hpp
//  \brief defines class Reconstruction, data and functions for spatial reconstruction

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// Forward declarations
class MeshBlock;
class ParameterInput;

//! \class Reconstruction
//  \brief member functions implement various spatial reconstruction algorithms

class Reconstruction {
public:
  Reconstruction(MeshBlock *pmb, ParameterInput *pin);
  ~Reconstruction();

  // data
  // order and type of reconstruction algorithm
  int xorder;   // order of hydro reconstruction
  bool characteristic_reconstruction;  // TRUE for characteristic recon
  bool uniform_limiter[3]; // TRUE to use the PLM or PPM limiter option w/o coord terms
  AthenaArray<Real> c1i,c2i,c3i,c4i,c5i,c6i;  // coefficients for PPM in x1
  AthenaArray<Real> hplus_ratio_i, hminus_ratio_i; // for curvilinear PPMx1
  AthenaArray<Real> c1j,c2j,c3j,c4j,c5j,c6j;  // coefficients for PPM in x2
  AthenaArray<Real> hplus_ratio_j, hminus_ratio_j; // for curvilinear PPMx2
  AthenaArray<Real> c1k,c2k,c3k,c4k,c5k,c6k;  // coefficients for PPM in x3
  AthenaArray<Real> hplus_ratio_k, hminus_ratio_k; // for curvilinear PPMx3

  // functions
  // linear transformations of vectors between primitive and characteristic variables
  static void LeftEigenmatrixDotVector(MeshBlock *pmb, const int ivx,
    const int il, const int iu, const AthenaArray<Real> &b1, const AthenaArray<Real> &w,
    AthenaArray<Real> &vect);
  static void RightEigenmatrixDotVector(MeshBlock *pmb, const int ivx,
    const int il, const int iu, const AthenaArray<Real> &b1, const AthenaArray<Real> &w,
    AthenaArray<Real> &vect);

  // reconstruction functions of various orders in each dimension
  static void DonorCellX1(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void DonorCellX2(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void DonorCellX3(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PiecewiseLinearX1(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PiecewiseLinearX2(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PiecewiseLinearX3(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PiecewiseParabolicX1(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PiecewiseParabolicX2(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PiecewiseParabolicX3(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

private:
  MeshBlock* pmy_block_;  // ptr to MeshBlock containing this Reconstruction

  // scratch arrays used in PLM and PPM reconstruction functions
  AthenaArray<Real> scr01_i_,scr02_i_,scr03_i_,scr04_i_,scr05_i_;
  AthenaArray<Real> scr06_i_,scr07_i_,scr08_i_,scr09_i_,scr10_i_;
  AthenaArray<Real> scr11_i_,scr12_i_,scr13_i_,scr14_i_;
  AthenaArray<Real> scr1_ni_, scr2_ni_, scr3_ni_, scr4_ni_, scr5_ni_;
  AthenaArray<Real> scr6_ni_, scr7_ni_, scr8_ni_;
};
#endif // RECONSTRUCT_RECONSTRUCTION_HPP_
