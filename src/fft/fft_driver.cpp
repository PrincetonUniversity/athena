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

FFTDriver::FFTDriver(Mesh *pm, ParameterInput *pin) {
  pmy_mesh_=pm;

  if (pm->use_uniform_meshgen_fn_[X1DIR]==false
      || pm->use_uniform_meshgen_fn_[X2DIR]==false
      || pm->use_uniform_meshgen_fn_[X3DIR]==false) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTDriver::FFTDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  dim_ = 1;
  if (pm->mesh_size.nx2 > 1) dim_=2;
  if (pm->mesh_size.nx3 > 1) dim_=3;

  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  // *** we also need to construct another neighbor list for Multigrid ***

  ranklist_  = new int[pm->nbtotal];
  for (int n=0; n<pm->nbtotal; n++)
    ranklist_[n]=pm->ranklist[n];

  nranks_  = Globals::nranks;
  nslist_  = new int[nranks_];
  nblist_  = new int[nranks_];
#ifdef MPI_PARALLEL
  MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_FFT);
#endif
  for (int n=0; n<nranks_; n++) {
    nslist_[n]  = pm->nslist[n];
    nblist_[n]  = pm->nblist[n];
  }

  int ns = nslist_[Globals::my_rank];
  int ne = ns+nblist_[Globals::my_rank];

  int64_t &lx1min = pm->loclist[ns].lx1;
  int64_t &lx2min = pm->loclist[ns].lx2;
  int64_t &lx3min = pm->loclist[ns].lx3;
  int64_t lx1max = lx1min;
  int64_t lx2max = lx2min;
  int64_t lx3max = lx3min;

  for (int n=ns; n<ne; n++) {
    int64_t &lx1 = pm->loclist[n].lx1;
    int64_t &lx2 = pm->loclist[n].lx2;
    int64_t &lx3 = pm->loclist[n].lx3;
    lx1min = lx1min<lx1?lx1min:lx1;
    lx2min = lx2min<lx2?lx2min:lx2;
    lx3min = lx3min<lx3?lx3min:lx3;
    lx1max = lx1max>lx1?lx1min:lx1;
    lx2max = lx2max>lx2?lx2min:lx2;
    lx3max = lx3max>lx3?lx3min:lx3;
  }

  int nbx1 = static_cast<int>(lx1max-lx1min+1);
  int nbx2 = static_cast<int>(lx2max-lx2min+1);
  int nbx3 = static_cast<int>(lx3max-lx3min+1);

  nmb = nbx1*nbx2*nbx3; // number of mesh blocks to be loaded to the FFT block
  if (pm->nbtotal/nmb != nranks_) {
    // Will be implemented later.
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTDriver::FFTDriver" << std::endl
        << nmb << " MeshBlocks will be loaded to the FFT block."  << std::endl
        << "Number of FFT blocks " << pm->nbtotal/nmb << " are not matched with "
        << "Number of processors " << nranks_ << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // There should be only one FFTBlock per processor.
  // MeshBlocks in the same processor will be gathered into FFTBlock of the processor.
  fft_loclist_ = new LogicalLocation[nranks_];

  for (int n=0; n<nranks_; n++) {
    int ns = nslist_[n];
    fft_loclist_[n] = pm->loclist[ns];
    fft_loclist_[n].lx1 = fft_loclist_[n].lx1/nbx1;
    fft_loclist_[n].lx2 = fft_loclist_[n].lx2/nbx2;
    fft_loclist_[n].lx3 = fft_loclist_[n].lx3/nbx3;
  }
  npx1 = static_cast<int>(pm->nrbx1/nbx1);
  npx2 = static_cast<int>(pm->nrbx2/nbx2);
  npx3 = static_cast<int>(pm->nrbx3/nbx3);

  fft_mesh_size_=pm->mesh_size;

  RegionSize &bsize = (pm->pblock->block_size);

  fft_block_size_.nx1=pm->mesh_size.nx1/npx1;
  fft_block_size_.nx2=pm->mesh_size.nx2/npx2;
  fft_block_size_.nx3=pm->mesh_size.nx3/npx3;

  Real x1size=bsize.x1max-bsize.x1min;
  Real x2size=bsize.x2max-bsize.x2min;
  Real x3size=bsize.x3max-bsize.x3min;
  fft_block_size_.x1min=bsize.x1min;
  fft_block_size_.x1max=bsize.x1min+x1size*nbx1;
  fft_block_size_.x2min=bsize.x2min;
  fft_block_size_.x2max=bsize.x2min+x2size*nbx2;
  fft_block_size_.x3min=bsize.x3min;
  fft_block_size_.x3max=bsize.x3min+x3size*nbx3;

  gcnt_ = fft_mesh_size_.nx1*fft_mesh_size_.nx2*fft_mesh_size_.nx3;

#ifdef MPI_PARALLEL
  decomp_ = 0; pdim_ = 0;
  if (npx1 > 1) {
    decomp_ = decomp_ | DecompositionNames::x_decomp;
    pdim_++;
  }
  if (npx2 > 1) {
    decomp_ = decomp_ | DecompositionNames::y_decomp;
    pdim_++;
  }
  if (npx3 > 1) {
    decomp_ = decomp_ | DecompositionNames::z_decomp;
    pdim_++;
  }
#endif

}

// destructor

FFTDriver::~FFTDriver() {
  delete [] ranklist_;
  delete [] nslist_;
  delete [] nblist_;
  delete [] fft_loclist_;
  delete pmy_fb;
}

void FFTDriver::InitializeFFTBlock(bool set_norm) {
  int igid=Globals::my_rank;
  pmy_fb=new FFTBlock(this, fft_loclist_[igid], igid, fft_mesh_size_, fft_block_size_);
  if (set_norm) pmy_fb->SetNormFactor(1./gcnt_);
}

void FFTDriver::QuickCreatePlan(void) {

  pmy_fb->fplan_=pmy_fb->QuickCreatePlan(pmy_fb->in_,AthenaFFTForward);
  pmy_fb->bplan_=pmy_fb->QuickCreatePlan(pmy_fb->in_,AthenaFFTBackward);

  return;
}
