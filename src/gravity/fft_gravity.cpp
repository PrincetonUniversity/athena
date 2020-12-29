//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft_gravity.cpp
//! \brief implementation of functions in class FFTGravity

// C headers

// C++ headers
#include <cmath>
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../fft/athena_fft.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../task_list/fft_grav_task_list.hpp"
#include "fft_gravity.hpp"
#include "gravity.hpp"

//----------------------------------------------------------------------------------------
//! \fn FFTGravityDriver::FFTGravityDriver(Mesh *pm, ParameterInput *pin)
//! \brief FFTGravityDriver constructor

FFTGravityDriver::FFTGravityDriver(Mesh *pm, ParameterInput *pin)
    : FFTDriver(pm, pin) {
  four_pi_G_ = pmy_mesh_->four_pi_G_;
  if (four_pi_G_ == 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // initialize using FFTGravity

  int igid=Globals::my_rank;
  pmy_fb = new FFTGravity(this, fft_loclist_[igid], igid, fft_mesh_size_,
                          fft_block_size_);
  pmy_fb->SetNormFactor(four_pi_G_/gcnt_);

  QuickCreatePlan();

  gtlist_ = new FFTGravitySolverTaskList(pin, pm);
}

FFTGravityDriver::~FFTGravityDriver() {
  delete gtlist_;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTGravityDriver::Solve(int stage, int mode)
//! \brief load the data and solve

void FFTGravityDriver::Solve(int stage, int mode) {
  FFTBlock *pfb = pmy_fb;
  AthenaArray<Real> in;
  // Load the source
  int nbs = nslist_[Globals::my_rank];
  int nbe = nbs+nblist_[Globals::my_rank]-1;
  for (int igid=nbs; igid<=nbe; igid++) {
    MeshBlock *pmb = pmy_mesh_->FindMeshBlock(igid);
    if (pmb != nullptr) {
      in.InitWithShallowSlice(pmb->phydro->u,4,IDN,1);
      pfb->LoadSource(in, 0, NGHOST, pmb->loc, pmb->block_size);
    }
    //    else { // on another process
    //    }
  }

  pfb->ExecuteForward();
  pfb->ApplyKernel(mode);
  pfb->ExecuteBackward();

  // Return the result
  for (int igid=nbs; igid<=nbe; igid++) {
    MeshBlock *pmb = pmy_mesh_->FindMeshBlock(igid);
    if (pmb != nullptr) {
      pfb->RetrieveResult(pmb->pgrav->phi, 0, NGHOST,
                          pmb->loc, pmb->block_size);
    }
    //    else { // on another process
    //    }
  }

  gtlist_->DoTaskListOneStage(pmy_mesh_, stage);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTGravity::ApplyKernel(int mode)
//! \brief Apply kernel
void FFTGravity::ApplyKernel(int mode) {
  Real pcoeff(0.0);
  Real dx1sq = SQR(TWO_PI/(kNx[0]*dkx[0]));
  Real dx2sq = SQR(TWO_PI/(kNx[1]*dkx[1]));
  Real dx3sq = SQR(TWO_PI/(kNx[2]*dkx[2]));
  for (int k=0; k<knx[2]; k++) {
    for (int j=0; j<knx[1]; j++) {
      for (int i=0; i<knx[0]; i++) {
        std::int64_t gidx = GetGlobalIndex(i,j,k);
        if (gidx == 0) {
          pcoeff = 0.0;
        } else {
          Real kx = (i+kdisp[0]);
          Real ky = (j+kdisp[1]);
          Real kz = (k+kdisp[2]);
          if (kx > 0.5*kNx[0]) kx -= kNx[0];
          if (ky > 0.5*kNx[1]) ky -= kNx[1];
          if (kz > 0.5*kNx[2]) kz -= kNx[2];
          if (mode == 0) { // Discrete FT
            kx *= TWO_PI/static_cast<Real>(kNx[0]);
            ky *= TWO_PI/static_cast<Real>(kNx[1]);
            kz *= TWO_PI/static_cast<Real>(kNx[2]);
            pcoeff = ((2.0*std::cos(kx)-2.0)/dx1sq);
            if (dim_ > 1) pcoeff += ((2.0*std::cos(ky)-2.0)/dx2sq);
            if (dim_ > 2) pcoeff += ((2.0*std::cos(kz)-2.0)/dx3sq);
          } else if (mode == 1) { // Continous FT
            kx *= dkx[0];
            ky *= dkx[1];
            kz *= dkx[2];
            pcoeff = -kx*kx;
            if (dim_ > 1) pcoeff -= ky*ky;
            if (dim_ > 2) pcoeff -= kz*kz;
          }
          pcoeff = 1.0/pcoeff;
        }

        std::int64_t idx_in = GetIndex(i,j,k,b_in_);
        std::int64_t idx_out = GetIndex(i,j,k,f_out_);
        in_[idx_in] = pcoeff*out_[idx_out];
      }
    }
  }
  return;
}
