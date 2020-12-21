#ifndef FFT_ATHENA_FFT_HPP_
#define FFT_ATHENA_FFT_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_fft.hpp
//! \brief defines FFT class which implements parallel FFT using MPI/OpenMP

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
#include "plimpton/fft_2d.h"
#include "plimpton/fft_3d.h"
#endif // MPI_PARALLEL
#endif

#ifdef FFT
#ifdef MPI_PARALLEL // parallel FFT
struct AthenaFFTPlan {
  struct fft_plan_3d *plan3d;
  struct fft_plan_2d *plan2d;
  fftw_plan plan;
  int dir;
  int dim;
};
#else // serial FFT
struct AthenaFFTPlan {
  fftw_plan plan;
  int dir;
  int dim;
};
#endif
#else // no FFT
struct AthenaFFTPlan {
  void *plan;
  int dir;
  int dim;
};
#endif

class Mesh;
class MeshBlock;
class ParameterInput;
class FFTBlock;
class FFTDriver;
class TurbulenceDriver;

class AthenaFFTIndex{
 public:
  explicit AthenaFFTIndex(const AthenaFFTIndex *psrc);
  AthenaFFTIndex(int dim, LogicalLocation loc, RegionSize msize, RegionSize bsize);

  // mesh size
  Real Lx[3];
  int Nx[3];
  // MPI decomposition
  int np[3], ip[3];
  // local size and indices
  int nx[3], is[3], ie[3];

  int iloc[3],ploc[3];

  void SetLocalIndex();

  void SwapAxis(int ref_axis);
  void PermuteAxis(int npermute);
  void SwapProc(int ref_axis);
  void PermuteProc(int npermute);
  void PrintIndex();

  friend class FFTDriver;
  friend class FFTBlock;
 private:
  int dim_;
  int npermute_, swap_;
  template <typename T> void Swap_(T loc[], int ref_axis);
  template <typename T> void Permute_(T loc[], int npermute);
};

//! \class FFTBlock
//! \brief

class FFTBlock {
 public:
  FFTBlock(FFTDriver *pfd, LogicalLocation iloc, int igid,
           RegionSize msize, RegionSize bsize);
  virtual ~FFTBlock();

  void LoadSource(const AthenaArray<Real> &src, bool nu, int ngh,
                  LogicalLocation loc, RegionSize bsize);
  void RetrieveResult(AthenaArray<Real> &dst, bool nu, int ngh,
                      LogicalLocation loc, RegionSize bsize);
  virtual void ApplyKernel(int mode);

  std::int64_t GetIndex(const int i, const int j, const int k);
  std::int64_t GetIndex(const int i, const int j, const int k, AthenaFFTIndex *pidx);

  std::int64_t GetGlobalIndex(const int i, const int j, const int k);

  void DestroyPlan(AthenaFFTPlan *plan);
  void InitializeMPI();
  void Execute(AthenaFFTPlan *plan);
  void Execute(AthenaFFTPlan *plan, std::complex<Real> *data);
  void Execute(AthenaFFTPlan *plan, std::complex<Real> *in_data,
               std::complex<Real> *out_data);

  enum class AthenaFFTDirection {forward=-1, backward=1};

  AthenaFFTPlan *QuickCreatePlan(std::complex<Real> *data,AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(int nfast, std::complex<Real> *data,
                            AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(int nfast, int nslow, std::complex<Real> *data,
                            AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(int nfast, int nmid, int nslow,
                            std::complex<Real> *data,
                            AthenaFFTDirection dir);

  void ExecuteForward() {Execute(fplan_);}
  void ExecuteBackward() {Execute(bplan_);}

  void PrintSource(int in);
  void PrintNormFactor() {std::cout << norm_factor_ << std::endl;}

  void SetNormFactor(Real norm) { norm_factor_=norm;}

  int Nx[3], nx[3], disp[3];
  int kNx[3], knx[3], kdisp[3];
  Real dkx[3], dx1, dx2, dx3;

  friend class TurbulenceDriver;
  friend class FFTDriver;
  friend class Mesh;

 protected:
  FFTDriver *pmy_driver_;
  std::int64_t cnt_, gcnt_;
  int gid_;
  std::complex<Real> *in_, *out_;
  AthenaFFTPlan *fplan_, *bplan_;
  AthenaFFTIndex *f_in_, *f_out_, *b_in_, *b_out_;
  Real norm_factor_;
  int dim_;

  LogicalLocation loc_;
  RegionSize msize_, bsize_;
  AthenaFFTIndex orig_idx_;
#ifdef MPI_PARALLEL
  int decomp_, pdim_;
  int permute0_, permute1_, permute2_;
  bool swap1_, swap2_;
#endif
};

//! \class FFTDriver
//! \brief FFT driver

class FFTDriver {
 public:
  FFTDriver(Mesh *pm, ParameterInput *pin);
  virtual ~FFTDriver();

  int npx1, npx2, npx3, nmb;
  FFTBlock *pmy_fb;

  void QuickCreatePlan();
  void InitializeFFTBlock(bool set_norm);
  // small functions
  int GetNumFFTBlocks() { return nblist_[Globals::my_rank]; }

  friend class FFTBlock;
  friend class Mesh;

 protected:
  std::int64_t gcnt_;
  int nranks_, nblocks_;
  int *ranklist_, *nslist_, *nblist_;
  Mesh *pmy_mesh_;
  RegionSize fft_mesh_size_, fft_block_size_;
  LogicalLocation *fft_loclist_;
#ifdef MPI_PARALLEL
  int decomp_, pdim_;
#endif
  const int dim_;
#ifdef MPI_PARALLEL
  MPI_Comm MPI_COMM_FFT;
#endif
};

#ifdef MPI_PARALLEL
namespace DecompositionNames{
const unsigned int x_decomp = 1<<0;
const unsigned int y_decomp = 1<<1;
const unsigned int z_decomp = 1<<2;
const unsigned int xy_decomp = x_decomp | y_decomp;
const unsigned int yz_decomp = y_decomp | z_decomp;
const unsigned int xz_decomp = x_decomp | z_decomp;
const unsigned int xyz_decomp = x_decomp | y_decomp | z_decomp;
};
#endif

#endif // FFT_ATHENA_FFT_HPP_
