#ifndef FFT_BLOCK_FFT_HPP_
#define FFT_BLOCK_FFT_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file block_fft.hpp
//! \brief defines interface class to Plimpton's fftMPI

// C headers

// C++ headers
#include <complex>
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/meshblock_tree.hpp"

#ifdef FFT
#include <fftw3.h>
#ifdef MPI_PARALLEL
#include <mpi.h>
#include "fftmpi/fft2d.h"
#include "fftmpi/fft3d.h"
#endif // MPI_PARALLEL
#endif


//! \class BlockFFT
//! \brief interface to the Plimpton's fftMPI

class BlockFFT {
 public:
  explicit BlockFFT(MeshBlock *pmb);
  virtual ~BlockFFT();

  // data
  const int ndim;                   // number of dimensions
  const int is, ie, js, je, ks, ke; // meshblock indices
  const int Nx1, Nx2, Nx3;          // mesh size (active zones)
  const int nx1, nx2, nx3;          // meshblock size (active zones)

  // global index for input data layout
  const int in_ilo, in_ihi, in_jlo, in_jhi, in_klo, in_khi;
  int fast_ilo, fast_ihi, fast_jlo, fast_jhi, fast_klo, fast_khi;
  int fast_nx1, fast_nx2, fast_nx3;
  int mid_ilo, mid_ihi, mid_jlo, mid_jhi, mid_klo, mid_khi;
  int mid_nx1, mid_nx2, mid_nx3;
  int slow_ilo, slow_ihi, slow_jlo, slow_jhi, slow_klo, slow_khi;
  int slow_nx1, slow_nx2, slow_nx3;
#ifdef FFT
#ifdef MPI_PARALLEL
  FFTMPI_NS::FFT3d *pf3d;
#endif
#endif

  // functions
  void LoadSource(const AthenaArray<Real> &src);
  void RetrieveResult(AthenaArray<Real> &dst);
  virtual void ExecuteForward();
  virtual void ExecuteBackward();
  virtual void ApplyKernel();

 protected:
  MeshBlock *pmy_block_;
  std::complex<Real> *in_;
};

#endif // FFT_BLOCK_FFT_HPP_
