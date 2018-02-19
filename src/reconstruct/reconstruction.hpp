#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP
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

  // order and type of reconstruction algorithm
  int xorder;
  bool characteristic_reconstruction;

  // functions to perform linear transformations of vectors between primitive and
  // characteristic variables
  static void LeftEigenmatrixDotVector(MeshBlock *pmb, const int ivx,
    const int il, const int iu, const AthenaArray<Real> &b1, const AthenaArray<Real> &w,
    AthenaArray<Real> &vect);
  static void VectorDotRightEigenmatrix(MeshBlock *pmb, const int ivx,
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

  void PiecewiseLinearX1(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseLinearX2(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseLinearX3(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PPMX1(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PPMX2(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PPMX3(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PPMUniformX1(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PPMUniformX2(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  static void PPMUniformX3(MeshBlock *pmb, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &w,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

private:
  MeshBlock* pmy_block_;  // ptr to MeshBlock containing this Reconstruction

  // 1D scratch arrays used in PLM and PPM reconstruction functions
  AthenaArray<Real> bx_;
  AthenaArray<Real> dph_, dph_p1_;
  AthenaArray<Real> qplus_, qminus_;
  AthenaArray<Real> dqf_plus_, dqf_minus_, d2qf_;
  AthenaArray<Real> d2qc_m1_, d2qc_, d2qc_p1_;

  // 2D scratch arrays used in PLM functions
  AthenaArray<Real> dwl_,dwr_,wc_,dw2_,dwm_;

};
#endif // RECONSTRUCTION_HPP
