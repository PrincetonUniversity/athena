//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fftgravity.cpp
//  \brief implementation of functions in class FFTGravity

// Athena++ headers
#include "fftgravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../fft/athena_fft.hpp"
#include "../globals.hpp"

//----------------------------------------------------------------------------------------
//! \fn FFTGravityDriver::FFTGravityDriver(Mesh *pm, ParameterInput *pin)
//  \brief FFTGravityDriver constructor

FFTGravityDriver::FFTGravityDriver(Mesh *pm, ParameterInput *pin)
{
  four_pi_G_=pmy_mesh_->four_pi_G_;
  if(four_pi_G_==0.0) {
   std::stringstream msg;
   msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void GravityDriver::Solve(int step)
//  \brief load the data and solve

void FFTGravityDriver::Solve(int step)
{
  FFTBlock *pfb=pfb_;
  AthenaArray<Real> in;

  // Load the source 
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
  for(int igid=nbs;igid<=nbe;igid++){
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(gid);
    if(pmb!=NULL) {
      in.InitWithShallowCopy(pmb->phydro->u);
      pfb->LoadSource(in, IDN, NGHOST, pmb->loc, pmb->block_size);
    }
//    else { // on another process
//    }
  }

  pfblock->Execute(pfblock->fplan); 
  pfblock->ApplyKernel(mode);
  pfblock->Execute(pfblock->bplan); 

  // Return the result
  for(int gid=nbs;gid<=nbe;gid++){
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(gid);
    if(pmb!=NULL) {
      pfblock->RetrieveResult(pmb->pgrav->phi, 0, NGHOST, 
                              pmb->loc, pmb->block_size);
    }
//    else { // on another process
//    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTBlock::ApplyKernel(const AthenaArray<Real> &src, int ns)
//  \brief Apply kernel
void FFTGravity::ApplyKernel(int mode)
{
  int is, ie, js, je, ks, ke;
  Real pcoeff;
  Real dx1sq=rdx_*rdx_;
  Real dx2sq=rdy_*rdy_;
  Real dx3sq=rdz_*rdz_;
  for(int k=0, k<knx_[2]; k++) {
    for(int j=0, j<knx_[1]; j++) {
      for(int i=0, i<knx_[0]; i++) {
        if(mode == 0){ // fully periodic in all directions
          long int gidx = GetGlobalIndex(i,j,k);
          if(gidx == 0){ pcoeffi = 0.0;}
          else {
            pcoeff = ((2.0*std::cos((i+kdisp_[0])*dkx_[0])-2.0)/dx1sq);
            if(pfft->dim_ > 1)
              pcoeff += ((2.0*std::cos((j+kdisp_[1])*dkx_[1])-2.0)/dx2sq);
            if(pfft->dim_ > 2)
              pcoeff += ((2.0*std::cos((k+kdisp_[2])*dkx_[2])-2.0)/dx3sq);
            pcoeff = 1.0/pcoeff;
          }
 
          long int idx_in=GetIndex(i,j,k,b_in_);
          long int idx_out=GetIndex(i,j,k,f_out_);
          in_[idx_in][0] = out_[idx_out][0]*pcoeff;
          in_[idx_in][1] = out_[idx_out][1]*pcoeff;
        } 
//        else if (mode == 1){ 
// some other cases
//        }
      }
    }
  }
  return;
}
