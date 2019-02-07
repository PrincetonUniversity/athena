//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.cpp
//  \brief functions that apply BCs for FACE_CENTERED variables

// C headers

// C++ headers
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
    MeshBlock *pmb, enum BoundaryType type, FaceField &var, EdgeField &var_flux)
    : BoundaryVariable(pmb, type) {
  // Is there no better way to copy a reference to a FaceField in a class data member?
  var_fc.x1f.InitWithShallowCopy(var.x1f);
  var_fc.x2f.InitWithShallowCopy(var.x2f);
  var_fc.x3f.InitWithShallowCopy(var.x3f);
  src.x1f.InitWithShallowCopy(var_fc.x1f);
  src.x2f.InitWithShallowCopy(var_fc.x2f);
  src.x3f.InitWithShallowCopy(var_fc.x3f);
  dst.x1f.InitWithShallowCopy(var_fc.x1f);
  dst.x2f.InitWithShallowCopy(var_fc.x2f);
  dst.x3f.InitWithShallowCopy(var_fc.x3f);

  e1.InitWithShallowCopy(var_flux.x1e);
  e2.InitWithShallowCopy(var_flux.x2e);
  e3.InitWithShallowCopy(var_flux.x3e);
  // KGF: Taken from 2x functions in flux_correction_fc.cpp
  // AthenaArray<Real> &e1=pmb->pfield->e.x1e;
  // AthenaArray<Real> &e2=pmb->pfield->e.x2e;
  // AthenaArray<Real> &e3=pmb->pfield->e.x3e;

  InitBoundaryData(bd_fc_flcor_, BNDRY_EMFCOR);
  pbd_var_flcor_ = &(bd_fc_flcor_);

  // KGF: was not in "if (MAGNETIC_FIELDS_ENABLED)" conditional in master
  if (pbval_->num_north_polar_blocks_ > 0) {
    emf_north_send_ = new Real *[pbval_->num_north_polar_blocks_];
    emf_north_recv_ = new Real *[pbval_->num_north_polar_blocks_];
    emf_north_flag_ = new enum BoundaryStatus[pbval_->num_north_polar_blocks_];
#ifdef MPI_PARALLEL
    req_emf_north_send_ = new MPI_Request[pbval_->num_north_polar_blocks_];
    req_emf_north_recv_ = new MPI_Request[pbval_->num_north_polar_blocks_];
#endif
    for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
      emf_north_send_[n] = nullptr;
      emf_north_recv_[n] = nullptr;
      emf_north_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
      req_emf_north_send_[n] = MPI_REQUEST_NULL;
      req_emf_north_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }
  if (pbval_->num_south_polar_blocks_ > 0) {
    emf_south_send_ = new Real *[pbval_->num_south_polar_blocks_];
    emf_south_recv_ = new Real *[pbval_->num_south_polar_blocks_];
    emf_south_flag_ = new enum BoundaryStatus[pbval_->num_south_polar_blocks_];
#ifdef MPI_PARALLEL
    req_emf_south_send_ = new MPI_Request[pbval_->num_south_polar_blocks_];
    req_emf_south_recv_ = new MPI_Request[pbval_->num_south_polar_blocks_];
#endif
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
      emf_south_send_[n] = nullptr;
      emf_south_recv_[n] = nullptr;
      emf_south_flag_[n] = BNDRY_WAITING;
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
  DestroyBoundaryData(bd_fc_flcor_);

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
  if (pmy_mesh_->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) ei++;
    else if (nb.ox1<0) si--;
  }
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0)      si=pmb->is,          ei=pmb->ie;
  else if (nb.ox1>0)  si=pmb->ie-NGHOST+1, ei=pmb->ie;
  else               si=pmb->is,          ei=pmb->is+NGHOST-1;
  if (pmb->block_size.nx2==1) sj=pmb->js,  ej=pmb->je;
  else if (nb.ox2==0) sj=pmb->js,          ej=pmb->je+1;
  else if (nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js+1,        ej=pmb->js+NGHOST;
  if (pmy_mesh_->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) ej++;
    else if (nb.ox2<0) sj--;
  }
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, sj, ej, sk, ek, p);

  // bx3
  if (nb.ox2==0)      sj=pmb->js,          ej=pmb->je;
  else if (nb.ox2>0)  sj=pmb->je-NGHOST+1, ej=pmb->je;
  else               sj=pmb->js,          ej=pmb->js+NGHOST-1;
  if (pmb->block_size.nx3==1) sk=pmb->ks,  ek=pmb->ke;
  else if (nb.ox3==0) sk=pmb->ks,          ek=pmb->ke+1;
  else if (nb.ox3>0)  sk=pmb->ke-NGHOST+1, ek=pmb->ke;
  else               sk=pmb->ks+1,        ek=pmb->ks+NGHOST;
  if (pmy_mesh_->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) ek++;
    else if (nb.ox3<0) sk--;
  }
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek, p);

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
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) ei++;
    else if (nb.ox1<0) si--;
  }
  pmr->RestrictFieldX1(src.x1f, pmr->coarse_b_.x1f, si, ei, sj, ej, sk, ek);
  BufferUtility::Pack3DData(pmr->coarse_b_.x1f, buf, si, ei, sj, ej, sk, ek, p);

  // bx2
  if (nb.ox1==0)      si=pmb->cis,       ei=pmb->cie;
  else if (nb.ox1>0)  si=pmb->cie-cng+1, ei=pmb->cie;
  else               si=pmb->cis,       ei=pmb->cis+cng-1;
  if (pmb->block_size.nx2==1) sj=pmb->cjs, ej=pmb->cje;
  else if (nb.ox2==0) sj=pmb->cjs,       ej=pmb->cje+1;
  else if (nb.ox2>0)  sj=pmb->cje-cng+1, ej=pmb->cje;
  else               sj=pmb->cjs+1,     ej=pmb->cjs+cng;
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) ej++;
    else if (nb.ox2<0) sj--;
  }
  pmr->RestrictFieldX2(src.x2f, pmr->coarse_b_.x2f, si, ei, sj, ej, sk, ek);
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
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) ek++;
    else if (nb.ox3<0) sk--;
  }
  pmr->RestrictFieldX3(src.x3f, pmr->coarse_b_.x3f, si, ei, sj, ej, sk, ek);
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
  BufferUtility::Pack3DData(src.x1f, buf, si, ei, sj, ej, sk, ek, p);

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
  BufferUtility::Pack3DData(src.x2f, buf, si, ei, sj, ej, sk, ek, p);

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
  BufferUtility::Pack3DData(src.x3f, buf, si, ei, sj, ej, sk, ek, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendBoundaryBuffers(void)
//  \brief Send face-centered boundary buffers

void FaceCenteredBoundaryVariable::SendBoundaryBuffers(void) {
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
      CopyBufferSameProcess(nb, ssize);
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
  if (pmy_mesh_->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox1>0) si--;
    else if (nb.ox1<0) ei++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_[IB1] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x1f, si, ei, sj, ej, sk, ek, p);
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
  if (pmy_mesh_->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) sj--;
    else if (nb.ox2<0) ej++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx2==1) { // 1D
#pragma omp simd
    for (int i=si; i<=ei; ++i)
      dst.x2f(sk,sj+1,i) = dst.x2f(sk,sj,i);
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
  if (pmy_mesh_->multilevel==true && nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) sk--;
    else if (nb.ox3<0) ek++;
  }
  if (nb.polar) {
    Real sign = flip_across_pole_[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(sk+1,j,i) = dst.x3f(sk,j,i);
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
  if (nb.type != NEIGHBOR_FACE) {
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
          dst.x1f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x1f, si, ei, sj, ej, sk, ek, p);
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
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox2>0) sj--;
    else if (nb.ox2<0) ej++;
  }

  if (nb.polar) {
    Real sign = flip_across_pole_[IB2] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x2f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x2f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx2==1) { // 1D
#pragma omp simd
    for (int i=si; i<=ei; ++i)
      dst.x2f(sk,sj+1,i) = dst.x2f(sk,sj,i);
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
  if (nb.type != NEIGHBOR_FACE) {
    if (nb.ox3>0) sk--;
    else if (nb.ox3<0) ek++;
  }

  if (nb.polar) {
    Real sign = flip_across_pole_[IB3] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
        for (int i=si; i<=ei; ++i)
          dst.x3f(k,j,i) = sign*buf[p++];
      }
    }
  } else {
    BufferUtility::Unpack3DData(buf, dst.x3f, si, ei, sj, ej, sk, ek, p);
  }
  if (pmb->block_size.nx3==1) { // 1D or 2D
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst.x3f(sk+1,j,i) = dst.x3f(sk,j,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool FaceCenteredBoundaryVariable::ReceiveBoundaryBuffers(void)
//  \brief receive the face-centered boundary data

bool FaceCenteredBoundaryVariable::ReceiveBoundaryBuffers(void) {
  bool bflag=true;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_var_.flag[nb.bufid]==BNDRY_ARRIVED) continue;
    if (bd_var_.flag[nb.bufid]==BNDRY_WAITING) {
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
        bd_var_.flag[nb.bufid] = BNDRY_ARRIVED;
      }
#endif
    }
  }
  return bflag;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetBoundaries(void)
//  \brief set the face-centered boundary data

void FaceCenteredBoundaryVariable::SetBoundaries(void) {
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
    bd_var_.flag[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if (pbval_->block_bcs[INNER_X2] == POLAR_BNDRY
      || pbval_->block_bcs[OUTER_X2] == POLAR_BNDRY) {
    PolarBoundarySingleAzimuthalBlock();
    PolarBoundaryAverageField();
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait(void)
//  \brief receive and set the face-centered boundary data for initialization

void FaceCenteredBoundaryVariable::ReceiveAndSetBoundariesWithWait(void) {
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
    bd_var_.flag[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if (pbval_->block_bcs[INNER_X2] == POLAR_BNDRY
      || pbval_->block_bcs[OUTER_X2] == POLAR_BNDRY) {
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
    if (pbval_->block_bcs[INNER_X2]==POLAR_BNDRY) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = dst.x1f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x1f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = dst.x2f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x2f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k)
            pbval_->azimuthal_shift_(k) = dst.x3f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x3f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
    }

    if (pbval_->block_bcs[OUTER_X2]==POLAR_BNDRY) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST+1; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = dst.x1f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x1f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->je+2; j<=pmb->je+NGHOST+1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k)
            pbval_->azimuthal_shift_(k) = dst.x2f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x2f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
          }
        }
      }
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k)
            pbval_->azimuthal_shift_(k) = dst.x3f(k,j,i);
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST+1; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
            dst.x3f(k,j,i) = pbval_->azimuthal_shift_(k_shift);
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
  if (pbval_->block_bcs[INNER_X2] == POLAR_BNDRY) {
    int j = pmb->js;
    for (int k=kl; k<=ku; ++k) {
      for (int i=il; i<=iu; ++i) {
        dst.x2f(k,j,i) = 0.5 * (dst.x2f(k,j-1,i) + dst.x2f(k,j+1,i));
      }
    }
  }
  if (pbval_->block_bcs[OUTER_X2] == POLAR_BNDRY) {
    int j = pmb->je + 1;
    for (int k=kl; k<=ku; ++k) {
      for (int i=il; i<=iu; ++i) {
        dst.x2f(k,j,i) = 0.5 * (dst.x2f(k,j-1,i) + dst.x2f(k,j+1,i));
      }
    }
  }
  return;
}

// ------- KGF: move to a separate file?

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

void FaceCenteredBoundaryVariable::StartReceivingAll(const Real time) {
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank!=Globals::my_rank) {
      MPI_Start(&(bd_field_.req_recv[nb.bufid]));
      if (nb.type==NEIGHBOR_FACE || nb.type==NEIGHBOR_EDGE) {
        if ((nb.level>mylevel) ||
            ((nb.level==mylevel) && ((nb.type==NEIGHBOR_FACE)
                                     || ((nb.type==NEIGHBOR_EDGE)
                                         && (edge_flag_[nb.eid]==true)))))
          MPI_Start(&(bd_emfcor_.req_recv[nb.bufid]));
      }
    }
  }

  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = polar_neighbor_north[n];
      if (nb.rank != Globals::my_rank) {
        MPI_Start(&req_emf_north_recv_[n]);
      }
    }
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = polar_neighbor_south[n];
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
    bd_var_.flag[nb.bufid] = BNDRY_WAITING;
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

void FaceCenteredBoundaryVariable::ClearBoundaryAll(void) {
    // Clear non-polar boundary communications
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_var_.flag[nb.bufid] = BNDRY_WAITING;
    if ((nb.type==NEIGHBOR_FACE) || (nb.type==NEIGHBOR_EDGE))
      bd_fc_flcor_.flag[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    MeshBlock *pmb=pmy_block_;
    if (nb.rank!=Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_var_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      if (nb.type==NEIGHBOR_FACE || nb.type==NEIGHBOR_EDGE) {
        if (nb.level < pmb->loc.level)
          MPI_Wait(&(bd_emfcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
        else if ((nb.level==pmb->loc.level)
                 && ((nb.type==NEIGHBOR_FACE)
                     || ((nb.type==NEIGHBOR_EDGE) && (edge_flag_[nb.eid]==true))))
          MPI_Wait(&(bd_emfcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      }
    } // KGF: end block on other MPI process
#endif
  } // KGF: end loop over pbval_->nneighbor
  // Clear polar boundary communications
  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    emf_north_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    PolarNeighborBlock &nb = polar_neighbor_north[n];
    if (nb.rank != Globals::my_rank)
        MPI_Wait(&req_emf_north_send_[n], MPI_STATUS_IGNORE);
#endif
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    emf_south_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    PolarNeighborBlock &nb = polar_neighbor_south[n];
    if (nb.rank != Globals::my_rank)
      MPI_Wait(&req_emf_south_send_[n], MPI_STATUS_IGNORE);
#endif
  }
  return;
}
