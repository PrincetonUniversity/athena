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

  void PiecewiseLinearX1(const int kl, const int ku, const int jl, const int ju  ,
    const int il, const int iu, const AthenaArray<Real> &q, const int nin, const int nout,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX2(const int kl, const int ku, const int jl, const int ju  ,
    const int il, const int iu, const AthenaArray<Real> &q, const int nin, const int nout,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX3(const int kl, const int ku, const int jl, const int ju  ,
    const int il, const int iu, const AthenaArray<Real> &q, const int nin, const int nout,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX1(const int kl, const int ku, const int jl, const int ju  ,
    const int il, const int iu, const AthenaArray<Real> &q, const int nin, const int nout,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX2(const int kl, const int ku, const int jl, const int ju  ,
    const int il, const int iu, const AthenaArray<Real> &q, const int nin, const int nout,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX3(const int kl, const int ku, const int jl, const int ju  ,
    const int il, const int iu, const AthenaArray<Real> &q, const int nin, const int nout,
    AthenaArray<Real> &ql, AthenaArray<Real> &qr);

private:
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this Reconstruction
};
#endif // RECONSTRUCTION_HPP
