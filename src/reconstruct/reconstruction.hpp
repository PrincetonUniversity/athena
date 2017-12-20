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

// Reconstruction function pointer prototype
typedef void (*ReconstructFunc_t) (Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

//! \class Reconstruction
//  \brief member functions implement various spatial reconstruction algorithms

class Reconstruction {
public:
  Reconstruction(MeshBlock *pmb, ParameterInput *pin);
  ~Reconstruction();

  // reconstruction function pointers in each direction
  ReconstructFunc_t ReconstructFuncX1, ReconstructFuncX2, ReconstructFuncX3;

  static void DonorCellX1(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void DonorCellX2(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void DonorCellX3(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PiecewiseLinearX1(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PiecewiseLinearX2(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PiecewiseLinearX3(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PiecewiseLinearUniformX1(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PiecewiseLinearUniformX2(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PiecewiseLinearUniformX3(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PPMX1(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PPMX2(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PPMX3(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PPMUniformX1(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PPMUniformX2(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  static void PPMUniformX3(Coordinates *pco, const int kl, const int ku,
    const int jl, const int ju  , const int il, const int iu, const AthenaArray<Real> &q,
    const int nin, const int nout, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

private:
  MeshBlock* pmy_block_;  // ptr to MeshBlock containing this Reconstruction

  // 1D scratch arrays used in PPM reconstruction functions
  AthenaArray<Real> dph_, dph_p1_;
  AthenaArray<Real> qplus_, qminus_;
  AthenaArray<Real> dqf_plus_, dqf_minus_, d2qf_;
  AthenaArray<Real> d2qc_m1_, d2qc_, d2qc_p1_;
};
#endif // RECONSTRUCTION_HPP
