//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.cpp
//  \brief functions that apply BCs for FACE_CENTERED variables

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "../bvals.hpp"
#include "bvals_fc.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor

FaceCenteredBoundaryVariable::FaceCenteredBoundaryVariable(
    MeshBlock *pmb, FaceField &var, EdgeField &var_flux)
    : BoundaryVariable(pmb) {
  // Is there no better way to copy a reference to a FaceField in a class data member?
  var_fc.x1f.InitWithShallowCopy(var.x1f);
  var_fc.x2f.InitWithShallowCopy(var.x2f);
  var_fc.x3f.InitWithShallowCopy(var.x3f);
  // src.x1f.InitWithShallowCopy(var_fc.x1f);
  // src.x2f.InitWithShallowCopy(var_fc.x2f);
  // src.x3f.InitWithShallowCopy(var_fc.x3f);
  // dst.x1f.InitWithShallowCopy(var_fc.x1f);
  // dst.x2f.InitWithShallowCopy(var_fc.x2f);
  // dst.x3f.InitWithShallowCopy(var_fc.x3f);

  e1.InitWithShallowCopy(var_flux.x1e);
  e2.InitWithShallowCopy(var_flux.x2e);
  e3.InitWithShallowCopy(var_flux.x3e);
  // KGF: Taken from 2x functions in flux_correction_fc.cpp
  // AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  // AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  // AthenaArray<Real> &e3=pmb->pfield->e.x3e;

  InitBoundaryData(bd_var_, BoundaryQuantity::fc);
  InitBoundaryData(bd_var_flcor_, BoundaryQuantity::fc_flcor);
  // pbd_var_flcor_ = &(bd_fc_flcor_);
#ifdef MPI_PARALLEL
  fc_phys_id_ = pbval_->ReserveTagVariableIDs(2);
  fc_flx_phys_id_ = fc_phys_id_ + 1;
  if (pbval_->num_north_polar_blocks_ > 0
      || pbval_->num_south_polar_blocks_ > 0) {
    fc_flx_pole_phys_id_ = pbval_->ReserveTagVariableIDs(1);
  }
#endif

  // KGF: was not in "if (MAGNETIC_FIELDS_ENABLED)" conditional in master
  if (pbval_->num_north_polar_blocks_ > 0) {
    emf_north_send_ = new Real *[pbval_->num_north_polar_blocks_];
    emf_north_recv_ = new Real *[pbval_->num_north_polar_blocks_];
    emf_north_flag_ = new BoundaryStatus[pbval_->num_north_polar_blocks_];
#ifdef MPI_PARALLEL
    req_emf_north_send_ = new MPI_Request[pbval_->num_north_polar_blocks_];
    req_emf_north_recv_ = new MPI_Request[pbval_->num_north_polar_blocks_];
#endif
    for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
      emf_north_send_[n] = nullptr;
      emf_north_recv_[n] = nullptr;
      emf_north_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      req_emf_north_send_[n] = MPI_REQUEST_NULL;
      req_emf_north_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }
  if (pbval_->num_south_polar_blocks_ > 0) {
    emf_south_send_ = new Real *[pbval_->num_south_polar_blocks_];
    emf_south_recv_ = new Real *[pbval_->num_south_polar_blocks_];
    emf_south_flag_ = new BoundaryStatus[pbval_->num_south_polar_blocks_];
#ifdef MPI_PARALLEL
    req_emf_south_send_ = new MPI_Request[pbval_->num_south_polar_blocks_];
    req_emf_south_recv_ = new MPI_Request[pbval_->num_south_polar_blocks_];
#endif
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
      emf_south_send_[n] = nullptr;
      emf_south_recv_[n] = nullptr;
      emf_south_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      req_emf_south_send_[n] = MPI_REQUEST_NULL;
      req_emf_south_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }

  // Allocate buffers for polar neighbor communication
  if (pbval_->num_north_polar_blocks_ > 0) {
    for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
      emf_north_send_[n] = new Real[pmb->block_size.nx1];
      emf_north_recv_[n] = new Real[pmb->block_size.nx1];
    }
  }
  if (pbval_->num_south_polar_blocks_ > 0) {
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
      emf_south_send_[n] = new Real[pmb->block_size.nx1];
      emf_south_recv_[n] = new Real[pmb->block_size.nx1];
    }
  }
}

// destructor

FaceCenteredBoundaryVariable::~FaceCenteredBoundaryVariable() {
  //MeshBlock *pmb=pmy_block_;
  DestroyBoundaryData(bd_var_);
  DestroyBoundaryData(bd_var_flcor_);

  // KGF: the following is in "if (MAGNETIC_FIELDS_ENABLED)" conditional in master, but
  // only the emf_south/north_send/recv_[] etc. variables are conditionally allocated
  if (pbval_->num_north_polar_blocks_ > 0) {
    for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
      delete[] emf_north_send_[n];
      delete[] emf_north_recv_[n];
#ifdef MPI_PARALLEL
      if (req_emf_north_send_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_north_send_[n]);
      if (req_emf_north_recv_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_north_recv_[n]);
#endif
    }
    delete[] emf_north_send_;
    delete[] emf_north_recv_;
    delete[] emf_north_flag_;
#ifdef MPI_PARALLEL
    delete[] req_emf_north_send_;
    delete[] req_emf_north_recv_;
#endif
  }
  if (pbval_->num_south_polar_blocks_ > 0) {
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
      delete[] emf_south_send_[n];
      delete[] emf_south_recv_[n];
#ifdef MPI_PARALLEL
      if (req_emf_south_send_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_south_send_[n]);
      if (req_emf_south_recv_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_south_recv_[n]);
#endif
    }
    delete[] emf_south_send_;
    delete[] emf_south_recv_;
    delete[] emf_south_flag_;
#ifdef MPI_PARALLEL
    delete[] req_emf_south_send_;
    delete[] req_emf_south_recv_;
#endif
  }
}


int FaceCenteredBoundaryVariable::ComputeVariableBufferSize(const NeighborIndexes& ni,
                                                            int cng) {
  MeshBlock *pmb=pmy_block_;
  int f2d = (pmb->block_size.nx2 > 1 ? 1 : 0);
  int f3d = (pmb->block_size.nx3 > 1 ? 1 : 0);
  int cng1, cng2, cng3;
  cng1=cng;
  cng2=cng*f2d;
  cng3=cng*f3d;

  int size1=((ni.ox1==0) ? (pmb->block_size.nx1+1):NGHOST)
            *((ni.ox2==0) ? (pmb->block_size.nx2):NGHOST)
            *((ni.ox3==0) ? (pmb->block_size.nx3):NGHOST);
  int size2=((ni.ox1==0) ? (pmb->block_size.nx1):NGHOST)
            *((ni.ox2==0) ? (pmb->block_size.nx2+f2d):NGHOST)
            *((ni.ox3==0) ? (pmb->block_size.nx3):NGHOST);
  int size3=((ni.ox1==0) ? (pmb->block_size.nx1):NGHOST)
            *((ni.ox2==0) ? (pmb->block_size.nx2):NGHOST)
            *((ni.ox3==0) ? (pmb->block_size.nx3+f3d):NGHOST);
  int size=size1+size2+size3;
  if (pmy_mesh_->multilevel) {
    if (ni.type!=NeighborConnect::face) {
      if (ni.ox1!=0) size1=size1/NGHOST*(NGHOST+1);
      if (ni.ox2!=0) size2=size2/NGHOST*(NGHOST+1);
      if (ni.ox3!=0) size3=size3/NGHOST*(NGHOST+1);
    }
    size=size1+size2+size3;
    int f2c1=((ni.ox1==0) ? ((pmb->block_size.nx1+1)/2+1):NGHOST)
             *((ni.ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                   *((ni.ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
    int f2c2=((ni.ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
             *((ni.ox2==0) ? ((pmb->block_size.nx2+1)/2+f2d)
               : NGHOST)
             *((ni.ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
    int f2c3=((ni.ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
             *((ni.ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
             *((ni.ox3==0) ? ((pmb->block_size.nx3+1)/2+f3d)
               : NGHOST);
    if (ni.type!=NeighborConnect::face) {
      if (ni.ox1!=0) f2c1=f2c1/NGHOST*(NGHOST+1);
      if (ni.ox2!=0) f2c2=f2c2/NGHOST*(NGHOST+1);
      if (ni.ox3!=0) f2c3=f2c3/NGHOST*(NGHOST+1);
    }
    int fsize=f2c1+f2c2+f2c3;
    int c2f1=
        ((ni.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1+1):cng+1)
        *((ni.ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
        *((ni.ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
    int c2f2=
        ((ni.ox1==0) ?((pmb->block_size.nx1+1)/2+cng1):cng)
        *((ni.ox2==0)?((pmb->block_size.nx2+1)/2+cng2+f2d):cng+1)
        *((ni.ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng);
    int c2f3=
        ((ni.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
        *((ni.ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng)
        *((ni.ox3==0) ?
          ((pmb->block_size.nx3+1)/2+cng3+f3d) : cng+1);
    int csize=c2f1+c2f2+c2f3;
    size=std::max(size,std::max(csize,fsize));
  }
  return size;
}

int FaceCenteredBoundaryVariable::ComputeFluxCorrectionBufferSize(
    const NeighborIndexes& ni, int cng) {
  MeshBlock *pmb=pmy_block_;
  int size=0;

  if (ni.type==NeighborConnect::face) {
    if (pmb->block_size.nx3>1) { // 3D
      if (ni.ox1!=0)
        size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                   +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
      else if (ni.ox2!=0)
        size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
             +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
      else
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                   +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
    } else if (pmb->block_size.nx2>1) { // 2D
      if (ni.ox1!=0)
        size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
      else
        size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
    } else { // 1D
      size=2;
          }
  } else if (ni.type==NeighborConnect::edge) {
    if (pmb->block_size.nx3>1) { // 3D
      if (ni.ox3==0) size=pmb->block_size.nx3;
      if (ni.ox2==0) size=pmb->block_size.nx2;
      if (ni.ox1==0) size=pmb->block_size.nx1;
    } else if (pmb->block_size.nx2>1) {
      size=1;
    }
  }
  return size;
}

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                               const NeighborBlock& nb)
//  \brief Set face-centered boundary buffers for sending to a block on the same level

int FaceCenteredBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int p=0;

  // bx1
  if (nb.ox1==0)     si=pmb->is,          ei=pmb->ie+1;
  else if (nb.ox1>0) si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else              si=pmb->is+1,        ei=pmb->is+NGHOST;
  if (nb.ox2==0)     sj=pmb->js,          ej=pmb->je;
  else if (nb.ox2>0) sj=pmb->je-NGHOST+1, ej=pmb->je;
  else              sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if (nb.ox3==0)     sk=pmb->ks,          ek=pmb->ke;
  else if (nb.ox3>0) sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else              sk=pmb->ks,          ek=pmb->ks+NGHOST-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmy_mesh_->multilevel==true && nb.type != NeighborConnect::face) {
    if (nb.ox1>0) ei++;
    else if (nb.ox1<0) si--;
  }
  BufferUtility::Pack3DData(var_fc.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0)      si=pmb->is,          ei=pmb->ie;
  else if (nb.ox1>0)  si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else               si=pmb->is,          ei=pmb->is+NGHOST-1;
  if (pmb->block_size.nx2==1) sj=pmb->js,  ej=pmb->je;
  else if (nb.ox2==0) sj=pmb->js,          ej=pmb->je+1;
  else if (nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js+1,        ej=pmb->js+NGHOST;
  if (pmy_mesh_->multilevel==true && nb.type != NeighborConnect::face) {
    if (nb.ox2>0) ej++;
    else if (nb.ox2<0) sj--;
  }
  BufferUtility::Pack3DData(var_fc.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0)      sj=pmb->js,          ej=pmb->je;
  else if (nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if (pmb->block_size.nx3==1) sk=pmb->ks,  ek=pmb->ke;
  else if (nb.ox3==0) sk=pmb->ks,          ek=pmb->ke+1;
  else if (nb.ox3>0)  sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else               sk=pmb->ks+1,        ek=pmb->ks+NGHOST;
  if (pmy_mesh_->multilevel==true && nb.type != NeighborConnect::face) {
    if (nb.ox3>0) ek++;
    else if (nb.ox3<0) sk--;
  }
  BufferUtility::Pack3DData(var_fc.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set face-centered boundary buffers for sending to a block on the coarser level

int FaceCenteredBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cng=NGHOST;
  int p=0;

  // bx1
  if (nb.ox1==0)     si=pmb->cis,       ei=pmb->cie+1;
  else if (nb.ox1>0) si=pmb->cie-cng+1, ei=pmb->cie;
  else              si=pmb->cis+1,     ei=pmb->cis+cng;
  if (nb.ox2==0)     sj=pmb->cjs,       ej=pmb->cje;
  else if (nb.ox2>0) sj=pmb->cje-cng+1, ej=pmb->cje;
  else              sj=pmb->cjs,       ej=pmb->cjs+cng-1;
  if (nb.ox3==0)     sk=pmb->cks,       ek=pmb->cke;
  else if (nb.ox3>0) sk=pmb->cke-cng+1, ek=pmb->cke;
  else              sk=pmb->cks,       ek=pmb->cks+cng-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NeighborConnect::face) {
    if (nb.ox1>0) ei++;
    else if (nb.ox1<0) si--;
  }
  pmr->RestrictFieldX1(var_fc.x1f, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek);
  BufferUtility::Pack3DData(pmr->coarse_b_.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0)      si=pmb->cis,       ei=pmb->cie;
  else if (nb.ox1>0)  si=pmb->cie-cng+1, ei=pmb->cie;
  else               si=pmb->cis,       ei=pmb->cis+cng-1;
  if (pmb->block_size.nx2==1) sj=pmb->cjs, ej=pmb->cje;
  else if (nb.ox2==0) sj=pmb->cjs,       ej=pmb->cje+1;
  else if (nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs+1,     ej=pmb->cjs+cng;
  if (nb.type != NeighborConnect::face) {
    if (nb.ox2>0) ej++;
    else if (nb.ox2<0) sj--;
  }
  pmr->RestrictFieldX2(var_fc.x2f, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek);
  if (pmb->block_size.nx2==1) { // 1D
    for (int i=si; i<=ei; i++)
      pmr->coarse_b_.x2f(sk,sj+1,i)=pmr->coarse_b_.x2f(sk,sj,i);
  }
  BufferUtility::Pack3DData(pmr->coarse_b_.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0)      sj=pmb->cjs,       ej=pmb->cje;
  else if (nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs,       ej=pmb->cjs+cng-1;
  if (pmb->block_size.nx3==1) sk=pmb->cks,  ek=pmb->cke;
  else if (nb.ox3==0) sk=pmb->cks,       ek=pmb->cke+1;
  else if (nb.ox3>0)  sk=pmb->cke-cng+1, ek=pmb->cke;
  else               sk=pmb->cks+1,     ek=pmb->cks+cng;
  if (nb.type != NeighborConnect::face) {
    if (nb.ox3>0) ek++;
    else if (nb.ox3<0) sk--;
  }
  pmr->RestrictFieldX3(var_fc.x3f, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek);
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; j++) {
      for (int i=si; i<=ei; i++)
        pmr->coarse_b_.x3f(sk+1,j,i)=pmr->coarse_b_.x3f(sk,j,i);
    }
  }
  BufferUtility::Pack3DData(pmr->coarse_b_.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set face-centered boundary buffers for sending to a block on the finer level

int FaceCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                            const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;
  int p=0;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  // bx1
  if (nb.ox1==0) {
    if (nb.fi1==1)   si=pmb->is+pmb->block_size.nx1/2-pmb->cnghost, ei=pmb->ie+1;
    else            si=pmb->is, ei=pmb->ie+1-pmb->block_size.nx1/2+pmb->cnghost;
  } else if (nb.ox1>0) { si=pmb->ie+1-pmb->cnghost, ei=pmb->ie+1;}
  else              si=pmb->is,                ei=pmb->is+pmb->cnghost;
  if (nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je-cn, ej=pmb->je;}
  else              sj=pmb->js,    ej=pmb->js+cn;
  if (nb.ox3==0) {
    sk=pmb->ks,    ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke-cn, ek=pmb->ke;}
  else              sk=pmb->ks,    ek=pmb->ks+cn;
  BufferUtility::Pack3DData(var_fc.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0) {
    if (nb.fi1==1)   si=pmb->is+pmb->block_size.nx1/2-pmb->cnghost, ei=pmb->ie;
    else            si=pmb->is, ei=pmb->ie-pmb->block_size.nx1/2+pmb->cnghost;
  } else if (nb.ox1>0) { si=pmb->ie-cn, ei=pmb->ie;}
  else              si=pmb->is,    ei=pmb->is+cn;
  if (nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      ej++;
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+1-pmb->cnghost, ej=pmb->je+1;}
  else              sj=pmb->js,                ej=pmb->js+pmb->cnghost;
  BufferUtility::Pack3DData(var_fc.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0) {
    sj=pmb->js,    ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
        else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je-cn, ej=pmb->je;}
  else              sj=pmb->js,    ej=pmb->js+cn;
  if (nb.ox3==0) {
    sk=pmb->ks,    ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      ek++;
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke+1-pmb->cnghost, ek=pmb->ke+1;}
  else              sk=pmb->ks,                ek=pmb->ks+pmb->cnghost;
  BufferUtility::Pack3DData(var_fc.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendBoundaryBuffers()
//  \brief Send face-centered boundary buffers

void FaceCenteredBoundaryVariable::SendBoundaryBuffers() {
  MeshBlock *pmb=pmy_block_;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    int ssize;
    if (nb.level==pmb->loc.level)
      // KGF: src
      ssize=LoadBoundaryBufferSameLevel(bd_var_.send[nb.bufid],nb);
    else if (nb.level<pmb->loc.level)
      // KGF: src
      ssize=LoadBoundaryBufferToCoarser(bd_var_.send[nb.bufid],nb);
    else
      // KGF: src
      ssize=LoadBoundaryBufferToFiner(bd_var_.send[nb.bufid], nb);
    if (nb.rank == Globals::my_rank) { // on the same process
      CopyVariableBufferSameProcess(nb, ssize);
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&(bd_var_.req_send[nb.bufid]));
#endif
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set face-centered boundary received from a block on the same level

void FaceCenteredBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  int p=0;
  // bx1
  // for uniform grid: face-neighbors take care of the overlapping faces
  if (nb.ox1==0)     si=pmb->is,        ei=pmb->ie+1;
  else if (nb.ox1>0) si=pmb->ie+2,      ei=pmb->ie+NGHOST+1;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if (nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if (nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if (nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmy_mesh_->multilevel==true && nb.type != NeighborConnect::face) {
    if (nb.ox1>0) si--;
    else if (nb.ox1<0) ei++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          var_fc.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, var_fc.x1f, si, ei, sj, ej, sk, ek, p);
  }

  // bx2
  if (nb.ox1==0)      si=pmb->is,         ei=pmb->ie;
  else if (nb.ox1>0)  si=pmb->ie+1,       ei=pmb->ie+NGHOST;
  else               si=pmb->is-NGHOST,  ei=pmb->is-1;
  if (pmb->block_size.nx2==1) sj=pmb->js, ej=pmb->je;
  else if (nb.ox2==0) sj=pmb->js,         ej=pmb->je+1;
  else if (nb.ox2>0)  sj=pmb->je+2,       ej=pmb->je+NGHOST+1;
  else               sj=pmb->js-NGHOST,  ej=pmb->js-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmy_mesh_->multilevel==true && nb.type != NeighborConnect::face) {
    if (nb.ox2>0) sj--;
    else if (nb.ox2<0) ej++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          var_fc.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, var_fc.x2f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx2==1) { // 1D
#pragma omp simd
    for (int i=si; i<=ei; ++i)
      var_fc.x2f(sk,sj+1,i) = var_fc.x2f(sk,sj,i);
  }

  // bx3
  if (nb.ox2==0)      sj=pmb->js,         ej=pmb->je;
  else if (nb.ox2>0)  sj=pmb->je+1,       ej=pmb->je+NGHOST;
  else               sj=pmb->js-NGHOST,  ej=pmb->js-1;
  if (pmb->block_size.nx3==1) sk=pmb->ks, ek=pmb->ke;
  else if (nb.ox3==0) sk=pmb->ks,         ek=pmb->ke+1;
  else if (nb.ox3>0)  sk=pmb->ke+2,       ek=pmb->ke+NGHOST+1;
  else               sk=pmb->ks-NGHOST,  ek=pmb->ks-1;
  // for SMR/AMR, always include the overlapping faces in edge and corner boundaries
  if (pmy_mesh_->multilevel==true && nb.type != NeighborConnect::face) {
    if (nb.ox3>0) sk--;
    else if (nb.ox3<0) ek++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          var_fc.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, var_fc.x3f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        var_fc.x3f(sk+1,j,i) = var_fc.x3f(sk,j,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set face-centered prolongation buffer received from a block on the same level

void FaceCenteredBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                          const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;
  int p=0;

  // bx1
  if (nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie+1;
    if ((pmb->loc.lx1 & 1LL)==0LL) ei+=cng;
    else             si-=cng;
  } else if (nb.ox1>0) {  si=pmb->cie+1,   ei=pmb->cie+1+cng;}
  else               si=pmb->cis-cng, ei=pmb->cis;
  if (nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL)==0LL) ej+=cng;
      else             sj-=cng;
    }
  } else if (nb.ox2>0) {  sj=pmb->cje+1,   ej=pmb->cje+cng;}
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if (nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL)==0LL) ek+=cng;
      else             sk-=cng;
    }
  } else if (nb.ox3>0) {  sk=pmb->cke+1,   ek=pmb->cke+cng;}
  else               sk=pmb->cks-cng, ek=pmb->cks-1;

  if (nb.polar) {
    Real sign = flip_across_pole_[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek, p);
  }

  // bx2
  if (nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie;
    if ((pmb->loc.lx1 & 1LL)==0LL) ei+=cng;
    else             si-=cng;
  } else if (nb.ox1>0) {  si=pmb->cie+1,   ei=pmb->cie+cng;}
  else               si=pmb->cis-cng, ei=pmb->cis-1;
  if (nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      ej++;
      if ((pmb->loc.lx2 & 1LL)==0LL) ej+=cng;
      else             sj-=cng;
    }
  } else if (nb.ox2>0) {  sj=pmb->cje+1,   ej=pmb->cje+1+cng;}
  else               sj=pmb->cjs-cng, ej=pmb->cjs;

  if (nb.polar) {
    Real sign = flip_across_pole_[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek, p);
    if (pmb->block_size.nx2 == 1) { // 1D
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        pmr->coarse_b_.x2f(sk,sj+1,i) = pmr->coarse_b_.x2f(sk,sj,i);
    }
  }

  // bx3
  if (nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL)==0LL) ej+=cng;
      else             sj-=cng;
    }
  } else if (nb.ox2>0) {  sj=pmb->cje+1,   ej=pmb->cje+cng;}
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if (nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      ek++;
      if ((pmb->loc.lx3 & 1LL)==0LL) ek+=cng;
      else             sk-=cng;
    }
  } else if (nb.ox3>0) {  sk=pmb->cke+1,   ek=pmb->cke+1+cng;}
  else               sk=pmb->cks-cng, ek=pmb->cks;

  if (nb.polar) {
    Real sign = flip_across_pole_[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek, p);
    if (pmb->block_size.nx3 == 1) { // 2D
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i)
          pmr->coarse_b_.x3f(sk+1,j,i) = pmr->coarse_b_.x3f(sk,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetFielBoundaryFromFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set face-centered boundary received from a block on the same level

void FaceCenteredBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;
  int p=0;

  // bx1
  if (nb.ox1==0) {
    si=pmb->is, ei=pmb->ie+1;
    if (nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  } else if (nb.ox1>0) { si=pmb->ie+2,      ei=pmb->ie+NGHOST+1;}
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NeighborConnect::face) {
    if (nb.ox1>0) si--;
    else if (nb.ox1<0) ei++;
  }
  if (nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+1,      ej=pmb->je+NGHOST;}
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke+1,      ek=pmb->ke+NGHOST;}
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  if (nb.polar) {
    Real sign = flip_across_pole_[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          var_fc.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, var_fc.x1f, si, ei, sj, ej, sk, ek, p);
  }

  // bx2
  if (nb.ox1==0) {
    si=pmb->is, ei=pmb->ie;
    if (nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  } else if (nb.ox1>0) { si=pmb->ie+1,      ei=pmb->ie+NGHOST;}
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if (nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      ej++;
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+2,      ej=pmb->je+NGHOST+1;}
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NeighborConnect::face) {
    if (nb.ox2>0) sj--;
    else if (nb.ox2<0) ej++;
  }

  if (nb.polar) {
    Real sign = flip_across_pole_[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          var_fc.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, var_fc.x2f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx2==1) { // 1D
#pragma omp simd
    for (int i=si; i<=ei; ++i)
      var_fc.x2f(sk,sj+1,i) = var_fc.x2f(sk,sj,i);
  }

  // bx3
  if (nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ox1!=0) {
        if (nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      } else {
        if (nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ox2>0) { sj=pmb->je+1,      ej=pmb->je+NGHOST;}
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if (nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      ek++;
      if (nb.ox1!=0 && nb.ox2!=0) {
        if (nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      } else {
        if (nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ox3>0) { sk=pmb->ke+2,      ek=pmb->ke+NGHOST+1;}
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  // include the overlapping faces in edge and corner boundaries
  if (nb.type != NeighborConnect::face) {
    if (nb.ox3>0) sk--;
    else if (nb.ox3<0) ek++;
  }

  if (nb.polar) {
    Real sign = flip_across_pole_[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          var_fc.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, var_fc.x3f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        var_fc.x3f(sk+1,j,i) = var_fc.x3f(sk,j,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool FaceCenteredBoundaryVariable::ReceiveBoundaryBuffers()
//  \brief receive the face-centered boundary data

bool FaceCenteredBoundaryVariable::ReceiveBoundaryBuffers() {
  bool bflag=true;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_var_.flag[nb.bufid]==BoundaryStatus::arrived) continue;
    if (bd_var_.flag[nb.bufid]==BoundaryStatus::waiting) {
      if (nb.rank==Globals::my_rank) { // on the same process
        bflag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&(bd_var_.req_recv[nb.bufid]),&test,MPI_STATUS_IGNORE);
        if (static_cast<bool>(test)==false) {
          bflag=false;
          continue;
        }
        bd_var_.flag[nb.bufid] = BoundaryStatus::arrived;
      }
#endif
    }
  }
  return bflag;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetBoundaries()
//  \brief set the face-centered boundary data

void FaceCenteredBoundaryVariable::SetBoundaries() {
  MeshBlock *pmb=pmy_block_;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.level==pmb->loc.level)
      // KGF: dst
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.level<pmb->loc.level)
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      // KGF: dst
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
    PolarBoundarySingleAzimuthalBlock();
    PolarBoundaryAverageField();
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait()
//  \brief receive and set the face-centered boundary data for initialization

void FaceCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait() {
  MeshBlock *pmb=pmy_block_;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank)
      MPI_Wait(&(bd_var_.req_recv[nb.bufid]),MPI_STATUS_IGNORE);
#endif
    if (nb.level==pmb->loc.level)
      // KGF: dst
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.level<pmb->loc.level)
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      // KGF: dst
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
    PolarBoundarySingleAzimuthalBlock();
    PolarBoundaryAverageField();
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void FaceCenteredBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
  MeshBlock *pmb=pmy_block_;
  if (pmb->loc.level == pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    if (pbval_->block_bcs[BoundaryFace::inner_x2]==BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = var_fc.x1f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            var_fc.x1f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = var_fc.x2f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            var_fc.x2f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k)
            pbval_->azimuthal_shift_(k) = var_fc.x3f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            var_fc.x3f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
    }

    if (pbval_->block_bcs[BoundaryFace::outer_x2]==BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = var_fc.x1f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            var_fc.x1f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->je+2; j<=pmb->je+NGHOST+1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = var_fc.x2f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            var_fc.x2f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k)
            pbval_->azimuthal_shift_(k) = var_fc.x3f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            var_fc.x3f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::PolarBoundaryAverageField()
//  \brief set theta-component of field along axis

void FaceCenteredBoundaryVariable::PolarBoundaryAverageField() {
  MeshBlock *pmb = pmy_block_;
  int il = pmb->is - NGHOST;
  int iu = pmb->ie + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmb->block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) {
    int j = pmb->js;
    for (int k=kl; k<=ku; ++k) {
      for (int i=il; i<=iu; ++i) {
        var_fc.x2f(k,j,i) = 0.5 * (var_fc.x2f(k,j-1,i) + var_fc.x2f(k,j+1,i));
      }
    }
  }
  if (pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
    int j = pmb->je + 1;
    for (int k=kl; k<=ku; ++k) {
      for (int i=il; i<=iu; ++i) {
        var_fc.x2f(k,j,i) = 0.5 * (var_fc.x2f(k,j-1,i) + var_fc.x2f(k,j+1,i));
      }
    }
  }
  return;
}

void FaceCenteredBoundaryVariable::CountFineEdges() {
  MeshBlock* pmb=pmy_block_;
  int &mylevel=pmb->loc.level;

  // KGF: shared logic of counting of neighbor MeshBlock coarse/same/refined
  // count the number of the fine meshblocks contacting on each edge
  int eid=0;
  if (pmb->block_size.nx2 > 1) {
    for (int ox2=-1; ox2<=1; ox2+=2) {
      for (int ox1=-1; ox1<=1; ox1+=2) {
        int nis, nie, njs, nje;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        int nf=0, fl=mylevel;
        for (int nj=njs; nj<=nje; nj++) {
          for (int ni=nis; ni<=nie; ni++) {
            if (pbval_->nblevel[1][nj+1][ni+1] > fl)
              fl++, nf=0;
            if (pbval_->nblevel[1][nj+1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }
  if (pmb->block_size.nx3 > 1) {
    for (int ox3=-1; ox3<=1; ox3+=2) {
      for (int ox1=-1; ox1<=1; ox1+=2) {
        int nis, nie, nks, nke;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for (int nk=nks; nk<=nke; nk++) {
          for (int ni=nis; ni<=nie; ni++) {
            if (pbval_->nblevel[nk+1][1][ni+1] > fl)
              fl++, nf=0;
            if (pbval_->nblevel[nk+1][1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
    for (int ox3=-1; ox3<=1; ox3+=2) {
      for (int ox2=-1; ox2<=1; ox2+=2) {
        int njs, nje, nks, nke;
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for (int nk=nks; nk<=nke; nk++) {
          for (int nj=njs; nj<=nje; nj++) {
            if (pbval_->nblevel[nk+1][nj+1][1] > fl)
              fl++, nf=0;
            if (pbval_->nblevel[nk+1][nj+1][1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }
  // end KGF: shared logic of counting of neighbor MeshBlock coarse/same/refined
}

// ------- KGF: move to a separate file?

void FaceCenteredBoundaryVariable::SetupPersistentMPI() {
  // KGF: can I move the following fn call somewhere else so only MPI actions occur?
  CountFineEdges();

#ifdef MPI_PARALLEL
  MeshBlock* pmb=pmy_block_;
  int &mylevel=pmb->loc.level;

  // KGF: a lot of repeated logic as in InitBoundaryData()
  int f2d=0, f3d=0;
  int cng, cng1, cng2, cng3;
  if (pmb->block_size.nx2 > 1) f2d=1;
  if (pmb->block_size.nx3 > 1) f3d=1;
  cng  = cng1 =pmb->cnghost;
  cng2 = cng*f2d;
  cng3 = cng*f3d;
  int ssize, rsize;
  int tag;
  // Initialize non-polar neighbor communications to other ranks
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank!=Globals::my_rank) {
      int size, csize, fsize;
      int size1=((nb.ox1==0) ? (pmb->block_size.nx1+1):NGHOST)
                *((nb.ox2==0) ? (pmb->block_size.nx2):NGHOST)
                *((nb.ox3==0) ? (pmb->block_size.nx3):NGHOST);
      int size2=((nb.ox1==0) ? (pmb->block_size.nx1):NGHOST)
                *((nb.ox2==0) ? (pmb->block_size.nx2+f2d):NGHOST)
                *((nb.ox3==0) ? (pmb->block_size.nx3):NGHOST);
      int size3=((nb.ox1==0) ? (pmb->block_size.nx1):NGHOST)
                *((nb.ox2==0) ? (pmb->block_size.nx2):NGHOST)
                *((nb.ox3==0) ? (pmb->block_size.nx3+f3d):NGHOST);
      size=size1+size2+size3;
      if (pmy_mesh_->multilevel==true) {
        if (nb.type!=NeighborConnect::face) {
          if (nb.ox1!=0) size1=size1/NGHOST*(NGHOST+1);
          if (nb.ox2!=0) size2=size2/NGHOST*(NGHOST+1);
          if (nb.ox3!=0) size3=size3/NGHOST*(NGHOST+1);
        }
        size=size1+size2+size3;
        int f2c1=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2+1):NGHOST)
                 *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                 *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
        int f2c2=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                 *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2+f2d):NGHOST)
                 *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
        int f2c3=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                 *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                 *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2+f3d):NGHOST);
        if (nb.type!=NeighborConnect::face) {
          if (nb.ox1!=0) f2c1=f2c1/NGHOST*(NGHOST+1);
          if (nb.ox2!=0) f2c2=f2c2/NGHOST*(NGHOST+1);
          if (nb.ox3!=0) f2c3=f2c3/NGHOST*(NGHOST+1);
        }
        fsize=f2c1+f2c2+f2c3;
        int c2f1=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1+1):cng+1)
                 *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
                 *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
        int c2f2=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
                 *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2+f2d):cng+1)
                 *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
        int c2f3=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
                 *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
                 *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3+f3d):cng+1);
        csize=c2f1+c2f2+c2f3;
      } // end of multilevel==true
      if (nb.level==mylevel) // same refinement level
        ssize=size, rsize=size;
      else if (nb.level<mylevel) // coarser
        ssize=fsize, rsize=csize;
      else // finer
        ssize=csize, rsize=fsize;

      // face-centered field: bd_var_
      tag=pbval_->CreateBvalsMPITag(nb.lid, nb.targetid, fc_phys_id_);
      if (bd_var_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_send[nb.bufid]);
      MPI_Send_init(bd_var_.send[nb.bufid],ssize,MPI_ATHENA_REAL,
                    nb.rank,tag,MPI_COMM_WORLD,&(bd_var_.req_send[nb.bufid]));
      tag=pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, fc_phys_id_);
      if (bd_var_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
        MPI_Request_free(&bd_var_.req_recv[nb.bufid]);
      MPI_Recv_init(bd_var_.recv[nb.bufid],rsize,MPI_ATHENA_REAL,
                    nb.rank,tag,MPI_COMM_WORLD,&(bd_var_.req_recv[nb.bufid]));

      // emf correction
      int f2csize;
      if (nb.type==NeighborConnect::face) { // face
        if (pmb->block_size.nx3 > 1) { // 3D
          if (nb.fid==BoundaryFace::inner_x1 || nb.fid==BoundaryFace::outer_x1) {
            size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                 +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
            f2csize=(pmb->block_size.nx2/2+1)*(pmb->block_size.nx3/2)
                    +(pmb->block_size.nx2/2)*(pmb->block_size.nx3/2+1);
          } else if (nb.fid==BoundaryFace::inner_x2 || nb.fid==BoundaryFace::outer_x2) {
            size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                 +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
            f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx3/2)
                    +(pmb->block_size.nx1/2)*(pmb->block_size.nx3/2+1);
          } else if (nb.fid==BoundaryFace::inner_x3 || nb.fid==BoundaryFace::outer_x3) {
            size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                 +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
            f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx2/2)
                    +(pmb->block_size.nx1/2)*(pmb->block_size.nx2/2+1);
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          if (nb.fid==BoundaryFace::inner_x1 || nb.fid==BoundaryFace::outer_x1) {
            size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
            f2csize=(pmb->block_size.nx2/2+1)+pmb->block_size.nx2/2;
          } else if (nb.fid==BoundaryFace::inner_x2 || nb.fid==BoundaryFace::outer_x2) {
            size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
            f2csize=(pmb->block_size.nx1/2+1)+pmb->block_size.nx1/2;
          }
        } else { // 1D
          size=f2csize=2;
        }
      } else if (nb.type==NeighborConnect::edge) { // edge
        if (pmb->block_size.nx3 > 1) { // 3D
          if (nb.eid>=0 && nb.eid<4) {
            size=pmb->block_size.nx3;
            f2csize=pmb->block_size.nx3/2;
          } else if (nb.eid>=4 && nb.eid<8) {
            size=pmb->block_size.nx2;
            f2csize=pmb->block_size.nx2/2;
          } else if (nb.eid>=8 && nb.eid<12) {
            size=pmb->block_size.nx1;
            f2csize=pmb->block_size.nx1/2;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          size=f2csize=1;
        }
      } else { // corner
        continue;
      }
      // field flux (emf) correction: bd_var_flcor_
      if (nb.level==mylevel) { // the same level
        if ((nb.type==NeighborConnect::face) || ((nb.type==NeighborConnect::edge)
                                         && (edge_flag_[nb.eid]==true))) {
          tag=pbval_->CreateBvalsMPITag(nb.lid, nb.targetid, fc_flx_phys_id_);
          if (bd_var_flcor_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
            MPI_Request_free(&bd_var_flcor_.req_send[nb.bufid]);
          MPI_Send_init(bd_var_flcor_.send[nb.bufid],size,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&(bd_var_flcor_.req_send[nb.bufid]));
          tag=pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, fc_flx_phys_id_);
          if (bd_var_flcor_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
            MPI_Request_free(&bd_var_flcor_.req_recv[nb.bufid]);
          MPI_Recv_init(bd_var_flcor_.recv[nb.bufid],size,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&(bd_var_flcor_.req_recv[nb.bufid]));
        }
      }
      if (nb.level>mylevel) { // finer neighbor
        tag=pbval_->CreateBvalsMPITag(pmb->lid, nb.bufid, fc_flx_phys_id_);
        if (bd_var_flcor_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
          MPI_Request_free(&bd_var_flcor_.req_recv[nb.bufid]);
        MPI_Recv_init(bd_var_flcor_.recv[nb.bufid],f2csize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&(bd_var_flcor_.req_recv[nb.bufid]));
      }
      if (nb.level<mylevel) { // coarser neighbor
        tag=pbval_->CreateBvalsMPITag(nb.lid, nb.targetid, fc_flx_phys_id_);
        if (bd_var_flcor_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
          MPI_Request_free(&bd_var_flcor_.req_send[nb.bufid]);
        MPI_Send_init(bd_var_flcor_.send[nb.bufid],f2csize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&(bd_var_flcor_.req_send[nb.bufid]));
      }
    } // neighbor block is on separate MPI process
  } // end loop over neighbors
  // Initialize polar neighbor communications to other ranks

  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pbval_->polar_neighbor_north[n];
    if (nb.rank != Globals::my_rank) {
      tag = pbval_->CreateBvalsMPITag(nb.lid, pmb->loc.lx3, fc_flx_pole_phys_id_);
      if (req_emf_north_send_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_north_send_[n]);
      MPI_Send_init(emf_north_send_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
                    nb.rank, tag, MPI_COMM_WORLD, &req_emf_north_send_[n]);
      tag = pbval_->CreateBvalsMPITag(pmb->lid, n, fc_flx_pole_phys_id_);
      if (req_emf_north_recv_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_north_recv_[n]);
      MPI_Recv_init(emf_north_recv_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
                    nb.rank, tag, MPI_COMM_WORLD, &req_emf_north_recv_[n]);
    }
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    const PolarNeighborBlock &nb = pbval_->polar_neighbor_south[n];
    if (nb.rank != Globals::my_rank) {
      tag = pbval_->CreateBvalsMPITag(nb.lid, pmb->loc.lx3, fc_flx_pole_phys_id_);
      if (req_emf_south_send_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_south_send_[n]);
      MPI_Send_init(emf_south_send_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
                    nb.rank, tag, MPI_COMM_WORLD, &req_emf_south_send_[n]);
      tag = pbval_->CreateBvalsMPITag(pmb->lid, n, fc_flx_pole_phys_id_);
      if (req_emf_south_recv_[n]!=MPI_REQUEST_NULL)
        MPI_Request_free(&req_emf_south_recv_[n]);
      MPI_Recv_init(emf_south_recv_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
                    nb.rank, tag, MPI_COMM_WORLD, &req_emf_south_recv_[n]);
    }
  }
#endif
  return;
}

void FaceCenteredBoundaryVariable::StartReceivingForInit(bool cons_and_field) {
#ifdef MPI_PARALLEL
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank != Globals::my_rank) {
      if (cons_and_field) {  // normal case
        MPI_Start(&(bd_var_.req_recv[nb.bufid]));
      }
    }
  }
#endif
  return;
}

void FaceCenteredBoundaryVariable::StartReceivingAll() {
  recv_flx_same_lvl_ = true;
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank!=Globals::my_rank) {
      MPI_Start(&(bd_var_.req_recv[nb.bufid]));
      if (nb.type==NeighborConnect::face || nb.type==NeighborConnect::edge) {
        if ((nb.level>mylevel) ||
            ((nb.level==mylevel) && ((nb.type==NeighborConnect::face)
                                     || ((nb.type==NeighborConnect::edge)
                                         && (edge_flag_[nb.eid]==true)))))
          MPI_Start(&(bd_var_flcor_.req_recv[nb.bufid]));
      }
    }
  }

  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = pbval_->polar_neighbor_north[n];
      if (nb.rank != Globals::my_rank) {
        MPI_Start(&req_emf_north_recv_[n]);
      }
    }
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = pbval_->polar_neighbor_south[n];
      if (nb.rank != Globals::my_rank) {
        MPI_Start(&req_emf_south_recv_[n]);
      }
    }
#endif
  return;
}

void FaceCenteredBoundaryVariable::ClearBoundaryForInit(bool cons_and_field) {
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank) {
      if (cons_and_field) {  // normal case
        MPI_Wait(&(bd_var_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      }
    }
#endif
  }
  return;
}

void FaceCenteredBoundaryVariable::ClearBoundaryAll() {
    // Clear non-polar boundary communications
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BoundaryStatus::waiting;
    if ((nb.type==NeighborConnect::face) || (nb.type==NeighborConnect::edge))
      bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    MeshBlock *pmb=pmy_block_;
    if (nb.rank!=Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      if (nb.type==NeighborConnect::face || nb.type==NeighborConnect::edge) {
        if (nb.level < pmb->loc.level)
          MPI_Wait(&(bd_var_flcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
        else if ((nb.level==pmb->loc.level)
                 && ((nb.type==NeighborConnect::face)
                     || ((nb.type==NeighborConnect::edge)
                         && (edge_flag_[nb.eid]==true))))
          MPI_Wait(&(bd_var_flcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      }
    } // KGF: end block on other MPI process
#endif
  } // KGF: end loop over pbval_->nneighbor
  // Clear polar boundary communications
  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    emf_north_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    PolarNeighborBlock &nb = pbval_->polar_neighbor_north[n];
    if (nb.rank != Globals::my_rank)
        MPI_Wait(&req_emf_north_send_[n], MPI_STATUS_IGNORE);
#endif
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    emf_south_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    PolarNeighborBlock &nb = pbval_->polar_neighbor_south[n];
    if (nb.rank != Globals::my_rank)
      MPI_Wait(&req_emf_south_send_[n], MPI_STATUS_IGNORE);
#endif
  }
  return;
}
