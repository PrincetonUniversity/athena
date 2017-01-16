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

// FFTW header

enum AthenaFFTDirection { AthenaFFTForward = -1, AthenaFFTBackward = 1 };

#ifdef MPI_PARALLEL
#include "mpifft.hpp"
#else
#include "fftw3.h"

typedef fftw_complex AthenaFFTComplex;
typedef int AthenaFFTInt;
typedef struct AthenaFFTPlan{
  fftw_plan plan;
  enum AthenaFFTDirection dir;
  int dim;
} AthenaFFTPlan;
#endif


class MeshBlock;
class ParameterInput;

//! \class AthenaFFT
//  \brief 

class AthenaFFT {
friend MeshBlock;
public:
  AthenaFFT(MeshBlock *pmb);
  ~AthenaFFT();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

  AthenaFFTPlan *QuickCreatePlan(enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(AthenaFFTInt nx1, AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(AthenaFFTInt nx1, AthenaFFTInt nx2, AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  AthenaFFTPlan *CreatePlan(AthenaFFTInt nx1, AthenaFFTInt nx2, AthenaFFTInt nx3, 
                            AthenaFFTComplex *data, 
                            enum AthenaFFTDirection dir);
  void Initialize();
  void MpiCleanup();
  void Execute(AthenaFFTPlan *plan);
  void Execute(AthenaFFTPlan *plan, AthenaFFTComplex *data);

  int GetIndex(const int i);
  int GetIndex(const int j, const int i);
  int GetIndex(const int k, const int j, const int i);

  long int cnt,gcnt;
  AthenaFFTPlan *fplan,*bplan;
  AthenaFFTComplex *work;
//  AthenaArray<AthenaFFTComplex> *work;
  
private:
#ifdef MPI_PARALLEL
  MPI_Comm comm_;
#endif
  int nthreads_;
  int np1_, np2_, np3_;
  int gis_,gie_;
  int gjs_,gje_;
  int gks_,gke_;
  int nx1_,nx2_,nx3_;
  AthenaFFTInt gnx1_,gnx2_,gnx3_;
};

#endif // ATHENA_FFT_HPP
