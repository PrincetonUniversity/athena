#ifndef RECONSTRUCT_RECONSTRUCTION_HPP_
#define RECONSTRUCT_RECONSTRUCTION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.hpp
//! \brief defines class Reconstruction, data and functions for spatial reconstruction

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// Forward declarations
class MeshBlock;
class ParameterInput;

//! \class Reconstruction
//! \brief member functions implement various spatial reconstruction algorithms

class Reconstruction {
 public:
  Reconstruction(MeshBlock *pmb, ParameterInput *pin);

  // data
  // switches for reconstruction method variants:
  int xorder;   // roughly the formal order of accuracy of overall reconstruction method
  bool characteristic_projection; // reconstruct on characteristic or primitive hydro vars
  bool uniform[3], curvilinear[2];
  // (Cartesian reconstruction formulas are used for x3 azimuthal coordinate in both
  // cylindrical and spherical-polar coordinates)

  // related fourth-order solver switches
  const bool correct_ic, correct_err; // used in Mesh::Initialize() and ProblemGenerator()

  // x1-sliced arrays of interpolation coefficients and limiter parameters:
  AthenaArray<Real> c1i, c2i, c3i, c4i, c5i, c6i;  // coefficients for PPM in x1
  AthenaArray<Real> hplus_ratio_i, hminus_ratio_i; // for curvilinear PPMx1
  AthenaArray<Real> c1j, c2j, c3j, c4j, c5j, c6j;  // coefficients for PPM in x2
  AthenaArray<Real> hplus_ratio_j, hminus_ratio_j; // for curvilinear PPMx2
  AthenaArray<Real> c1k, c2k, c3k, c4k, c5k, c6k;  // coefficients for PPM in x3
  AthenaArray<Real> hplus_ratio_k, hminus_ratio_k; // for curvilinear PPMx3

  // functions
  // linear transformations of vectors between primitive and characteristic variables
  void LeftEigenmatrixDotVector(
      const int ivx, const int il, const int iu,
      const AthenaArray<Real> &b1, const AthenaArray<Real> &w, AthenaArray<Real> &vect);
  void RightEigenmatrixDotVector(
      const int ivx, const int il, const int iu,
      const AthenaArray<Real> &b1, const AthenaArray<Real> &w, AthenaArray<Real> &vect);

  // reconstruction functions of various orders in each dimension
  void DonorCellX1(const int k, const int j, const int il, const int iu,
                   const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                   AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void DonorCellX2(const int k, const int j, const int il, const int iu,
                   const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                   AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void DonorCellX3(const int k, const int j, const int il, const int iu,
                   const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                   AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseLinearX1(const int k, const int j, const int il, const int iu,
                         const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                         AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseLinearX2(const int k, const int j, const int il, const int iu,
                         const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                         AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseLinearX3(const int k, const int j, const int il, const int iu,
                         const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                         AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseParabolicX1(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                            AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseParabolicX2(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                            AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void PiecewiseParabolicX3(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
                            AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  // overloads for non-fluid (cell-centered Hydro prim. and magnetic field) reconstruction
  void DonorCellX1(const int k, const int j, const int il, const int iu,
                   const AthenaArray<Real> &q,
                   AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX2(const int k, const int j, const int il, const int iu,
                   const AthenaArray<Real> &q,
                   AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX3(const int k, const int j, const int il, const int iu,
                   const AthenaArray<Real> &q,
                   AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX1(const int k, const int j, const int il, const int iu,
                         const AthenaArray<Real> &q,
                         AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX2(const int k, const int j, const int il, const int iu,
                         const AthenaArray<Real> &q,
                         AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX3(const int k, const int j, const int il, const int iu,
                         const AthenaArray<Real> &q,
                         AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseParabolicX1(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &q,
                            AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseParabolicX2(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &q,
                            AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseParabolicX3(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &q,
                            AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  // overloads for cell centered variables with memory order as [k,j,i,n]
  // Notice that the default order in Athena++ is [n,k,j,i]
  void DonorCellX1(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &q, const int array_order,
                   AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX2(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &q, const int array_order,
                   AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellX3(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &q, const int array_order,
                   AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX1(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &q, const int array_order,
                         AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX2(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &q, const int array_order,
                         AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearX3(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &q, const int array_order,
                         AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseParabolicX1(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &q, const int array_order,
                            AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseParabolicX2(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &q, const int array_order,
                            AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseParabolicX3(const int k, const int j, const int il, const int iu,
                            const AthenaArray<Real> &q, const int array_order,
                            AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellZeta(NRRadiation *prad, const int zs, const int ze,
      AthenaArray<Real> &q,
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void DonorCellPsi(NRRadiation *prad, const int ps, const int pe,
      AthenaArray<Real> &q,
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void PiecewiseLinearZeta(NRRadiation *prad,
      const int zs, const int ze, AthenaArray<Real> &q,
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);


  void PiecewiseLinearPsi(NRRadiation *prad,
      const int zs, const int ze, AthenaArray<Real> &q,
      AthenaArray<Real> &ql, AthenaArray<Real> &qr);

 private:
  MeshBlock* pmy_block_;  // ptr to MeshBlock containing this Reconstruction

  // scratch arrays used in PLM and PPM reconstruction functions
  AthenaArray<Real> scr01_i_, scr02_i_, scr03_i_, scr04_i_, scr05_i_;
  AthenaArray<Real> scr06_i_, scr07_i_, scr08_i_, scr09_i_, scr10_i_;
  AthenaArray<Real> scr11_i_, scr12_i_, scr13_i_, scr14_i_;
  AthenaArray<Real> scr1_ni_, scr2_ni_, scr3_ni_, scr4_ni_, scr5_ni_;
  AthenaArray<Real> scr6_ni_, scr7_ni_, scr8_ni_;

  // scratch arrays for arrays with different ordering
  int nvar_;// the maximum numver of variables for reconstruction
  AthenaArray<Real> scr1_in_, scr2_in_, scr3_in_, scr4_in_, scr5_in_;
  AthenaArray<Real> scr1_in2_, scr2_in2_, scr3_in2_, scr4_in2_;
  AthenaArray<Real> scr1_nn_, scr2_nn_, scr3_nn_, scr4_nn_;
  AthenaArray<Real> scr6_in_, scr7_in_, scr8_in_;
};
#endif // RECONSTRUCT_RECONSTRUCTION_HPP_
