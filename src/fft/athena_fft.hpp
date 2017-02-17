#ifndef ATHENA_FFT_HPP
#define ATHENA_FFT_HPP

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_fft.hpp
//  \brief defines FFT class which implements parallel FFT using MPI/OpenMP

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

enum AthenaFFTDirection { AthenaFFTForward = -1, AthenaFFTBackward = 1 };

#ifdef FFT
#include "fftw3.h"
typedef fftw_complex AthenaFFTComplex;

#ifdef MPI_PARALLEL
#include "mpi.h"
#include "plimpton/fft_3d.h"
#include "plimpton/fft_2d.h"
typedef fftw_complex AthenaFFTComplex;
typedef struct AthenaFFTPlan{
  struct fft_plan_3d *plan3d;
  struct fft_plan_2d *plan2d;
  fftw_plan plan;
  int dir;
  int dim;
} AthenaFFTPlan;
#else
typedef struct AthenaFFTPlan{
  fftw_plan plan;
  int dir;
  int dim;
} AthenaFFTPlan;
#endif
#else // FFT
typedef Real AthenaFFTComplex[2];
typedef struct AthenaFFTPlan{
  void *plan;
  int dir;
  int dim;
} AthenaFFTPlan;
#endif // FFT

class MeshBlock;
class ParameterInput;
class Gravity;

class AthenaFFTIndex{
friend class AthenaFFT;
public:
  AthenaFFTIndex(int dim, MeshBlock *pmb);
  ~AthenaFFTIndex();

  AthenaFFTIndex(const AthenaFFTIndex *psrc);
  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

// mesh size
  int *Nx;
// MPI decomposition
  int *np, *ip;
// local size and indices
  int *nx, *is, *ie;

  int iloc[3],ploc[3];

  void SetLocalIndex();

  void SwapAxis(int ref_axis);
  void PermuteAxis(int npermute);
  void SwapProc(int ref_axis);
  void PermuteProc(int npermute);
  void RemapAxis(int dir);
  void RemapProc(int dir);
private:
  int dim_;
  int npermute_, swap_;
  void Permute_(int loc[], int npermute);
  void Swap_(int loc[], int ref_axis);
  void RemapArray_(int arr[], int loc[], int dir);
  void PrintIndex(void);
};



//! \class AthenaFFT
//  \brief 

class AthenaFFT {
friend class MeshBlock;
friend class Gravity;
public:
  AthenaFFT(MeshBlock *pmb);
  ~AthenaFFT();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

  int Nx[3],nx[3],disp[3];
  int knx[3],kdisp[3];
  Real dkx[3];

  AthenaFFTPlan *QuickCreatePlan(enum AthenaFFTDirection dir,AthenaFFTComplex *work);
  AthenaFFTPlan *CreatePlan(int nfast, AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(int nfast, int nslow, AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(int nfast, int nmid, int nslow, 
                            AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  void MpiInitialize();
  void Execute(AthenaFFTPlan *plan);
  void Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data);
  void Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in_data, AthenaFFTComplex *out_data);

  long int GetIndex(const int i, const int j, const int k);
  long int GetIndex(const int i, const int j, const int k, AthenaFFTIndex *pidx);

  long int GetGlobalIndex(const int i, const int j, const int k);

  long int cnt,gcnt;
  AthenaFFTPlan *fplan,*bplan;
  AthenaFFTComplex *in, *out;
  AthenaFFTIndex *orig_idx;
  AthenaFFTIndex *f_in,*f_out,*b_in,*b_out;
//  AthenaArray<AthenaFFTComplex> *work;
  
private:
#ifdef MPI_PARALLEL
  int permute0_, permute1_, permute2_;
  bool swap1_,swap2_;
  unsigned int decomp_, pdim_;
#endif
  unsigned int dim_;
  int nthreads_;
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

#endif // ATHENA_FFT_HPP
