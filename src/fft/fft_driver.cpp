//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft_driver.cpp
//  \brief implementation of functions in class FFTDriver

// C headers

// C++ headers
#include <cmath>
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "athena_fft.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

FFTDriver::FFTDriver(Mesh *pm, ParameterInput *pin) : nranks_(Globals::nranks),
                                                      pmy_mesh_(pm), dim_(pm->ndim) {
  if (!(pm->use_uniform_meshgen_fn_[X1DIR])
      || !(pm->use_uniform_meshgen_fn_[X2DIR])
      || !(pm->use_uniform_meshgen_fn_[X3DIR])) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTDriver::FFTDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  // *** we also need to construct another neighbor list for Multigrid ***

  ranklist_  = new int[pm->nbtotal];
  for (int n=0; n<pm->nbtotal; n++)
    ranklist_[n] = pm->ranklist[n];

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
  int ne = ns + nblist_[Globals::my_rank];

  std::int64_t &lx1min = pm->loclist[ns].lx1;
  std::int64_t &lx2min = pm->loclist[ns].lx2;
  std::int64_t &lx3min = pm->loclist[ns].lx3;
  std::int64_t lx1max = lx1min;
  std::int64_t lx2max = lx2min;
  std::int64_t lx3max = lx3min;

  int current_level = pm->root_level;

  for (int n=ns; n<ne; n++) {
    std::int64_t &lx1 = pm->loclist[n].lx1;
    std::int64_t &lx2 = pm->loclist[n].lx2;
    std::int64_t &lx3 = pm->loclist[n].lx3;
    lx1min = lx1min < lx1 ? lx1min : lx1;
    lx2min = lx2min < lx2 ? lx2min : lx2;
    lx3min = lx3min < lx3 ? lx3min : lx3;
    lx1max = lx1max > lx1 ? lx1min : lx1;
    lx2max = lx2max > lx2 ? lx2min : lx2;
    lx3max = lx3max > lx3 ? lx3min : lx3;
    if(pm->loclist[n].level > current_level) current_level = pm->loclist[n].level;
  }
  int ref_lev = current_level - pm->root_level;
  // KGF: possible underflow from std::int64_t
  int nbx1 = static_cast<int>(lx1max-lx1min+1);
  int nbx2 = static_cast<int>(lx2max-lx2min+1);
  int nbx3 = static_cast<int>(lx3max-lx3min+1);

  nmb = nbx1*nbx2*nbx3; // number of MeshBlock to be loaded to the FFT block
  if (pm->nbtotal/nmb != nranks_) {
    // this restriction (to a single cuboid FFTBlock that covers the union of all
    // MeshBlocks owned by an MPI rank) should be relaxed in the future
    std::stringstream msg;
    msg << "### FATAL ERROR in FFTDriver::FFTDriver" << std::endl
        << nmb << " MeshBlocks will be loaded to the FFT block."  << std::endl
        << "Number of FFT blocks " << pm->nbtotal/nmb << " are not matched with "
        << "Number of processors " << nranks_ << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // There should be only one FFTBlock per processor.
  // MeshBlocks in the same processor will be gathered into FFTBlock of the processor.
  fft_loclist_ = new LogicalLocation[nranks_];

  for (int n=0; n<nranks_; n++) {
    int ns_inner = nslist_[n];
    fft_loclist_[n] = pm->loclist[ns_inner];
    fft_loclist_[n].lx1 = fft_loclist_[n].lx1/nbx1;
    fft_loclist_[n].lx2 = fft_loclist_[n].lx2/nbx2;
    fft_loclist_[n].lx3 = fft_loclist_[n].lx3/nbx3;
  }
  npx1 = (pm->nrbx1*(1 << ref_lev))/nbx1;
  npx2 = (pm->nrbx2*(1 << ref_lev))/nbx2;
  npx3 = (pm->nrbx3*(1 << ref_lev))/nbx3;

  fft_mesh_size_ = pm->mesh_size;

  fft_mesh_size_.nx1 = pm->mesh_size.nx1*(1<<ref_lev);
  fft_mesh_size_.nx2 = pm->mesh_size.nx2*(1<<ref_lev);
  fft_mesh_size_.nx3 = pm->mesh_size.nx3*(1<<ref_lev);

  fft_block_size_.nx1 = fft_mesh_size_.nx1/npx1;
  fft_block_size_.nx2 = fft_mesh_size_.nx2/npx2;
  fft_block_size_.nx3 = fft_mesh_size_.nx3/npx3;

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
  int igid = Globals::my_rank;
  pmy_fb = new FFTBlock(this, fft_loclist_[igid], igid, fft_mesh_size_, fft_block_size_);
  if (set_norm) pmy_fb->SetNormFactor(1./gcnt_);
}

void FFTDriver::QuickCreatePlan() {
  pmy_fb->fplan_ = pmy_fb->QuickCreatePlan(
      pmy_fb->in_, FFTBlock::AthenaFFTDirection::forward);
  pmy_fb->bplan_ = pmy_fb->QuickCreatePlan(
      pmy_fb->in_, FFTBlock::AthenaFFTDirection::backward);

  return;
}
