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

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

// FFTW header

enum AthenaFFTDirection { AthenaFFTForward = -1, AthenaFFTBackward = 1 };

#include "mpifft.hpp"

class MeshBlock;
class ParameterInput;
class Gravity;

typedef struct AthenaFFTIndex{
// mesh size
  int N[3];
// decomposition
  int np[3],ip[3];
// local size and indices
  int nx[3],is[3],ie[3];
} AthenaFFTIndex;

//! \class AthenaFFT
//  \brief 

class AthenaFFT {
friend class MeshBlock;
friend class Gravity;
public:
  AthenaFFT(MeshBlock *pmb);
  ~AthenaFFT();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

  int dim;
  int idisp,jdisp,kdisp;
  int idisp_k,jdisp_k,kdisp_k;
  int nx1,nx2,nx3;
  int knx1,knx2,knx3;
  Real dkx,dky,dkz; 

  int gnfast, gnmid, gnslow;
  AthenaFFTIndex f_in,f_out,b_in,b_out;
  int permute0, permute1, permute2;
  bool swap1,swap2;
  unsigned int decomp; 

  AthenaFFTPlan *QuickCreatePlan(enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(AthenaFFTInt nfast, AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(AthenaFFTInt nfast, AthenaFFTInt nslow, AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(AthenaFFTInt nfast, AthenaFFTInt nmid, AthenaFFTInt nslow, 
                            AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  void Initialize();
  void MpiCleanup();
  void Execute(AthenaFFTPlan *plan);
  void Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data);
  void Execute(AthenaFFTPlan *plan, AthenaFFTComplex *in , AthenaFFTComplex *out);
  void CompatabilityCheck(int verbose);

  long int GetIndex(const int i, const int j, const int k);
  long int GetIndex(const int i, const int j, const int k, bool swap);
  long int GetFreq(const int i, const int j, const int k);
  long int GetFreq(const int i, const int j, const int k, bool swap);

  long int GetGlobalIndex(const int i, const int j, const int k);

  long int cnt,gcnt;
  AthenaFFTPlan *fplan,*bplan;
  AthenaFFTComplex *work;
//  AthenaArray<AthenaFFTComplex> *work;
  
private:
#ifdef MPI_PARALLEL
  int decomp_; // decomposition type
  MPI_Comm cart_2d_comm_; // communicator; only used for pencil decomp. lik ACCFFT
#endif
  int nthreads_;
  int np1_, np2_, np3_;
  AthenaFFTInt gnx1_,gnx2_,gnx3_;
};

namespace DecompositionNames{
  const unsigned int x_decomp = 1<<0;
  const unsigned int y_decomp = 1<<1;
  const unsigned int z_decomp = 1<<2;
  const unsigned int xy_decomp = x_decomp | y_decomp;
  const unsigned int yz_decomp = y_decomp | z_decomp;
  const unsigned int xz_decomp = x_decomp | z_decomp;
  const unsigned int xyz_decomp = x_decomp | y_decomp | z_decomp;
};

#endif // ATHENA_FFT_HPP
