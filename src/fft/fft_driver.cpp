//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft_driver.cpp
//  \brief implementation of functions in class FFTDriver

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "athena_fft.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

FFTDriver::FFTDriver(Mesh *pm, ParameterInput *pin)
{
  pmy_mesh_=pm;

  if(pm->use_meshgen_fn_[X1DIR]==true || pm->use_meshgen_fn_[X2DIR]==true
  || pm->use_meshgen_fn_[X3DIR]==true) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTDriver::FFTDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  dim_ = 1;
  if(pm->mesh_size.nx2 > 1) dim_=2;
  if(pm->mesh_size.nx3 > 1) dim_=3;
 
  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  // *** we also need to construct another neighbor list for Multigrid ***

  ranklist_  = new int[pm->nbtotal];
  for(int n=0; n<pm->nbtotal; n++)
    ranklist_[n]=pm->ranklist[n];

  nranks_  = Globals::nranks;
  nslist_  = new int[nranks_];
  nblist_  = new int[nranks_];
#ifdef MPI_PARALLEL
  MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_FFT);
#endif
  for(int n=0; n<nranks_; n++) {
    nslist_[n]  = pm->nslist[n];
    nblist_[n]  = pm->nblist[n];
  }

  // There should be only one FFTBlock per processor.
  // MeshBlocks in the same processor will be gathered into FFTBlock of the processor.
  fft_loclist_ = new LogicalLocation[nranks_];
  if(pm->nbtotal == nranks_){
    // If there is only one MeshBlock per processor, just copy its logical location
    for(int n=0; n<nranks_; n++) {
      fft_loclist_[n] = pm->loclist[n];
    }
    npx1=pm->nrbx1;
    npx2=pm->nrbx2;
    npx3=pm->nrbx3;
  } else {
    // Will be implemented later.
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTDriver::FFTDriver" << std::endl
        << "Currently, FFT solver only support one MeshBlock per processor." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  fft_mesh_size_.nx1=pm->mesh_size.nx1;
  fft_mesh_size_.nx2=pm->mesh_size.nx2;
  fft_mesh_size_.nx3=pm->mesh_size.nx3;
  fft_block_size_.nx1=pm->mesh_size.nx1/npx1;
  fft_block_size_.nx2=pm->mesh_size.nx2/npx2;
  fft_block_size_.nx3=pm->mesh_size.nx3/npx3;


#ifdef MPI_PARALLEL
  {using namespace DecompositionNames;
  decomp_ = 0; pdim_ = 0;
  if(npx1 > 1){
    decomp_ = decomp_ | x_decomp;
    pdim_++;
  }
  if(npx2 > 1){
    decomp_ = decomp_ | y_decomp;
    pdim_++;
  }
  if(npx3 > 1){
    decomp_ = decomp_ | z_decomp;
    pdim_++;
  }}
#endif

  int igid=nslist_[Globals::my_rank];
  pmy_fb=new FFTBlock(this, fft_loclist_[igid], igid, fft_mesh_size_, fft_block_size_);
}

// destructor

FFTDriver::~FFTDriver()
{
  delete [] nslist_;
  delete [] nblist_;
  delete pmy_fb;
}

void FFTDriver::QuickCreatePlan(void){
  
  pmy_fb->fplan_=pmy_fb->QuickCreatePlan(pmy_fb->in_,AthenaFFTForward);
  pmy_fb->bplan_=pmy_fb->QuickCreatePlan(pmy_fb->in_,AthenaFFTBackward);

  return;
}
