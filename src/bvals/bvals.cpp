//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.cpp
//  \brief constructor/destructor and utility functions for BoundaryValues class

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/mggravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"


// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

BoundaryValues::BoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs,
                               ParameterInput *pin)
 : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs) {
  pmy_block_=pmb;
  for (int i=0; i<6; i++)
    BoundaryFunction_[i]=NULL;

// Set BC functions for each of the 6 boundaries in turn ---------------------------------
  // Inner x1
  nface_=2; nedge_=0;
  switch(block_bcs[INNER_X1]) {
    case REFLECTING_BNDRY:
      BoundaryFunction_[INNER_X1] = ReflectInnerX1;
      break;
    case OUTFLOW_BNDRY:
      BoundaryFunction_[INNER_X1] = OutflowInnerX1;
      break;
    case BLOCK_BNDRY: // block boundary
    case PERIODIC_BNDRY: // periodic boundary
      BoundaryFunction_[INNER_X1] = NULL;
      break;
    case SHEAR_PERIODIC_BNDRY: // shearing periodic boundary
      if (!SHEARING_BOX) block_bcs[INNER_X1]=PERIODIC_BNDRY;
      BoundaryFunction_[INNER_X1] = NULL;
      break;
    case USER_BNDRY: // user-enrolled BCs
      BoundaryFunction_[INNER_X1] = pmy_mesh_->BoundaryFunction_[INNER_X1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ix1_bc=" << block_bcs[INNER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
      break;
   }

  // Outer x1
  switch(block_bcs[OUTER_X1]) {
    case REFLECTING_BNDRY:
      BoundaryFunction_[OUTER_X1] = ReflectOuterX1;
      break;
    case OUTFLOW_BNDRY:
      BoundaryFunction_[OUTER_X1] = OutflowOuterX1;
      break;
    case BLOCK_BNDRY: // block boundary
    case PERIODIC_BNDRY: // periodic boundary
      BoundaryFunction_[OUTER_X1] = NULL;
      break;
    case SHEAR_PERIODIC_BNDRY: // shearing periodic boundary
      if (!SHEARING_BOX) block_bcs[OUTER_X1]=PERIODIC_BNDRY;
      BoundaryFunction_[OUTER_X1] = NULL;
      break;
    case USER_BNDRY: // user-enrolled BCs
      BoundaryFunction_[OUTER_X1] = pmy_mesh_->BoundaryFunction_[OUTER_X1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ox1_bc=" << block_bcs[OUTER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  if (pmb->block_size.nx2 > 1) {
    nface_=4; nedge_=4;
    // Inner x2
    switch(block_bcs[INNER_X2]) {
      case REFLECTING_BNDRY:
        BoundaryFunction_[INNER_X2] = ReflectInnerX2;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[INNER_X2] = OutflowInnerX2;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
      case POLAR_BNDRY: // polar boundary
        BoundaryFunction_[INNER_X2] = NULL;
        break;
      case POLAR_BNDRY_WEDGE: //polar boundary with a wedge
        BoundaryFunction_[INNER_X2] = PolarWedgeInnerX2;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[INNER_X2] = pmy_mesh_->BoundaryFunction_[INNER_X2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix2_bc=" << block_bcs[INNER_X2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

    // Outer x2
    switch(block_bcs[OUTER_X2]) {
      case REFLECTING_BNDRY:
        BoundaryFunction_[OUTER_X2] = ReflectOuterX2;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[OUTER_X2] = OutflowOuterX2;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
      case POLAR_BNDRY: // polar boundary
        BoundaryFunction_[OUTER_X2] = NULL;
        break;
      case POLAR_BNDRY_WEDGE: // polar boundary with a wedge
        BoundaryFunction_[OUTER_X2] = PolarWedgeOuterX2;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[OUTER_X2] = pmy_mesh_->BoundaryFunction_[OUTER_X2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox2_bc=" << block_bcs[OUTER_X2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  if (pmb->block_size.nx3 > 1) {
    nface_=6; nedge_=12;
    // Inner x3
    switch(block_bcs[INNER_X3]) {
      case REFLECTING_BNDRY:
        BoundaryFunction_[INNER_X3] = ReflectInnerX3;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[INNER_X3] = OutflowInnerX3;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
        BoundaryFunction_[INNER_X3] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[INNER_X3] = pmy_mesh_->BoundaryFunction_[INNER_X3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix3_bc=" << block_bcs[INNER_X3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

    // Outer x3
    switch(block_bcs[OUTER_X3]) {
      case REFLECTING_BNDRY:
        BoundaryFunction_[OUTER_X3] = ReflectOuterX3;
        break;
      case OUTFLOW_BNDRY:
        BoundaryFunction_[OUTER_X3] = OutflowOuterX3;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
        BoundaryFunction_[OUTER_X3] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        BoundaryFunction_[OUTER_X3] = pmy_mesh_->BoundaryFunction_[OUTER_X3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox3_bc=" << block_bcs[OUTER_X3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  // Count number of blocks wrapping around pole
  if (block_bcs[INNER_X2] == POLAR_BNDRY || block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    if (pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Number of MeshBlocks around the pole must be 1 or even." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    int level = pmb->loc.level - pmy_mesh_->root_level;
    // possible loss of precision to 32 bit int, if int64_t nrbx3 is large
    num_north_polar_blocks_ = static_cast<int>(pmy_mesh_->nrbx3 * (1 << level));
  } else {
    num_north_polar_blocks_ = 0;
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY || block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
    if (pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Number of MeshBlocks around the pole must be 1 or even." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    int level = pmb->loc.level - pmy_mesh_->root_level;
    // possible loss of precision to 32 bit int, if int64_t nrbx3 is large
    num_south_polar_blocks_ = static_cast<int>(pmy_mesh_->nrbx3 * (1 << level));
  } else {
    num_south_polar_blocks_ = 0;
  }

  InitBoundaryData(bd_hydro_, BNDRY_HYDRO);
  if (pmy_mesh_->multilevel==true) // SMR or AMR
    InitBoundaryData(bd_flcor_, BNDRY_FLCOR);
  if (MAGNETIC_FIELDS_ENABLED) {
    InitBoundaryData(bd_field_, BNDRY_FIELD);
    InitBoundaryData(bd_emfcor_, BNDRY_EMFCOR);
  }

  if (num_north_polar_blocks_ > 0) {
    emf_north_send_ = new Real *[num_north_polar_blocks_];
    emf_north_recv_ = new Real *[num_north_polar_blocks_];
    emf_north_flag_ = new enum BoundaryStatus[num_north_polar_blocks_];
#ifdef MPI_PARALLEL
    req_emf_north_send_ = new MPI_Request[num_north_polar_blocks_];
    req_emf_north_recv_ = new MPI_Request[num_north_polar_blocks_];
#endif
    for (int n = 0; n < num_north_polar_blocks_; ++n) {
      emf_north_send_[n] = NULL;
      emf_north_recv_[n] = NULL;
      emf_north_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
      req_emf_north_send_[n] = MPI_REQUEST_NULL;
      req_emf_north_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }
  if (num_south_polar_blocks_ > 0) {
    emf_south_send_ = new Real *[num_south_polar_blocks_];
    emf_south_recv_ = new Real *[num_south_polar_blocks_];
    emf_south_flag_ = new enum BoundaryStatus[num_south_polar_blocks_];
#ifdef MPI_PARALLEL
    req_emf_south_send_ = new MPI_Request[num_south_polar_blocks_];
    req_emf_south_recv_ = new MPI_Request[num_south_polar_blocks_];
#endif
    for (int n = 0; n < num_south_polar_blocks_; ++n) {
      emf_south_send_[n] = NULL;
      emf_south_recv_[n] = NULL;
      emf_south_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
      req_emf_south_send_[n] = MPI_REQUEST_NULL;
      req_emf_south_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }

  // Allocate buffers for polar neighbor communication
  if (MAGNETIC_FIELDS_ENABLED) {
    if (num_north_polar_blocks_ > 0) {
      for (int n = 0; n < num_north_polar_blocks_; ++n) {
        emf_north_send_[n] = new Real[pmb->block_size.nx1];
        emf_north_recv_[n] = new Real[pmb->block_size.nx1];
      }
    }
    if (num_south_polar_blocks_ > 0) {
      for (int n = 0; n < num_south_polar_blocks_; ++n) {
        emf_south_send_[n] = new Real[pmb->block_size.nx1];
        emf_south_recv_[n] = new Real[pmb->block_size.nx1];
      }
    }
  }

 /* single CPU in the azimuthal direction with the polar boundary*/
  if (pmb->loc.level == pmy_mesh_->root_level &&
     pmy_mesh_->nrbx3 == 1 &&
     (block_bcs[INNER_X2]==POLAR_BNDRY||block_bcs[OUTER_X2]==POLAR_BNDRY||
      block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE||block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE))
       exc_.NewAthenaArray(pmb->ke+NGHOST+2);

// set parameters for shearing box bc and allocate buffers
  if (SHEARING_BOX) {
    Mesh *pmy_mesh = pmb->pmy_mesh;
    Omega_0_ = pin->GetOrAddReal("problem","Omega0",0.001);
    qshear_  = pin->GetOrAddReal("problem","qshear",1.5);
    ShBoxCoord_ = pin->GetOrAddInteger("problem","shboxcoord",1);
    x1size_ = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
    x2size_ = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
    x3size_ = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;
    int level = pmb->loc.level - pmy_mesh->root_level;
    int64_t nrbx1 = pmy_mesh->nrbx1*(1L << level);
    int64_t nrbx2 = pmy_mesh->nrbx2*(1L << level);

    shbb_.outer = false;
    shbb_.inner = false;

    if (ShBoxCoord_ == 1) {
      int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
      int ncells3 = pmb->block_size.nx3;
      if (pmy_mesh->mesh_size.nx3>1) ncells3 += 2*NGHOST;
      ssize_ = NGHOST*ncells3;

      if (pmb->loc.lx1 == 0) { // if true for shearing inner blocks
        if (block_bcs[INNER_X1] != SHEAR_PERIODIC_BNDRY) {
          block_bcs[INNER_X1] = SHEAR_PERIODIC_BNDRY;
          BoundaryFunction_[INNER_X1] = NULL;
        }
        shboxvar_inner_hydro_.NewAthenaArray(NHYDRO,ncells3,ncells2,NGHOST);
        flx_inner_hydro_.NewAthenaArray(ncells2);
        if (MAGNETIC_FIELDS_ENABLED) {
          shboxvar_inner_field_.x1f.NewAthenaArray(ncells3,ncells2,NGHOST);
          shboxvar_inner_field_.x2f.NewAthenaArray(ncells3,ncells2+1,NGHOST);
          shboxvar_inner_field_.x3f.NewAthenaArray(ncells3+1,ncells2,NGHOST);
          flx_inner_field_.x1f.NewAthenaArray(ncells2);
          flx_inner_field_.x2f.NewAthenaArray(ncells2+1);
          flx_inner_field_.x3f.NewAthenaArray(ncells2);
          shboxvar_inner_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
          shboxvar_inner_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
          shboxmap_inner_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
          shboxmap_inner_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
          flx_inner_emf_.x2e.NewAthenaArray(ncells2);
          flx_inner_emf_.x3e.NewAthenaArray(ncells2+1);
        }
        shbb_.inner = true;
        shbb_.igidlist=new int[nrbx2];
        shbb_.ilidlist=new int[nrbx2];
        shbb_.irnklist=new int[nrbx2];
        shbb_.ilevlist=new int[nrbx2];
        // attach corner cells from L/R side
        int size = (pmb->block_size.nx2+NGHOST)*ssize_*NHYDRO;
        int bsize=0, esize=0;
        if (MAGNETIC_FIELDS_ENABLED) {
          // extra cell in azimuth/vertical
          bsize = (pmb->block_size.nx2+NGHOST+1)*(ssize_+NGHOST)*NFIELD;
          // face plus edge for EMF
          esize = 2*(pmb->block_size.nx2+NGHOST)*pmb->block_size.nx3
                +pmb->block_size.nx2+pmb->block_size.nx3+NGHOST;
        }
        for (int n=0; n<2; n++) {
          send_innerbuf_hydro_[n] = new Real[size];
          recv_innerbuf_hydro_[n] = new Real[size];
          shbox_inner_hydro_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
          rq_innersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_innerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_innerbuf_field_[n] = new Real[bsize];
            recv_innerbuf_field_[n] = new Real[bsize];
            shbox_inner_field_flag_[n]=BNDRY_WAITING;
            send_innerbuf_emf_[n] = new Real[esize];
            recv_innerbuf_emf_[n] = new Real[esize];
            shbox_inner_emf_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
            rq_innersend_field_[n] = MPI_REQUEST_NULL;
            rq_innerrecv_field_[n] = MPI_REQUEST_NULL;
            rq_innersend_emf_[n] = MPI_REQUEST_NULL;
            rq_innerrecv_emf_[n] = MPI_REQUEST_NULL;
#endif
          }
        }
        size = NGHOST*ssize_*NHYDRO;// corner cells only
        if (MAGNETIC_FIELDS_ENABLED) {
            bsize = NGHOST*(ssize_+NGHOST)*NFIELD;
            esize = 2*NGHOST*pmb->block_size.nx3+NGHOST;
        }
        for (int n=2; n<4; n++) {
          send_innerbuf_hydro_[n] = new Real[size];
          recv_innerbuf_hydro_[n] = new Real[size];
          shbox_inner_hydro_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
          rq_innersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_innerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_innerbuf_field_[n] = new Real[bsize];
            recv_innerbuf_field_[n] = new Real[bsize];
            shbox_inner_field_flag_[n]=BNDRY_WAITING;
            send_innerbuf_emf_[n] = new Real[esize];
            recv_innerbuf_emf_[n] = new Real[esize];
            shbox_inner_emf_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
            rq_innersend_field_[n] = MPI_REQUEST_NULL;
            rq_innerrecv_field_[n] = MPI_REQUEST_NULL;
            rq_innersend_emf_[n] = MPI_REQUEST_NULL;
            rq_innerrecv_emf_[n] = MPI_REQUEST_NULL;
#endif
          }
        }
      }

      if (pmb->loc.lx1 == (nrbx1-1)) { // if true for shearing outer blocks
        if (block_bcs[OUTER_X1] != SHEAR_PERIODIC_BNDRY) {
          block_bcs[OUTER_X1] = SHEAR_PERIODIC_BNDRY;
          BoundaryFunction_[OUTER_X1] = NULL;
        }
        shboxvar_outer_hydro_.NewAthenaArray(NHYDRO,ncells3,ncells2,NGHOST);
        flx_outer_hydro_.NewAthenaArray(ncells2);
        if (MAGNETIC_FIELDS_ENABLED) {
          shboxvar_outer_field_.x1f.NewAthenaArray(ncells3,ncells2,NGHOST);
          shboxvar_outer_field_.x2f.NewAthenaArray(ncells3,ncells2+1,NGHOST);
          shboxvar_outer_field_.x3f.NewAthenaArray(ncells3+1,ncells2,NGHOST);
          flx_outer_field_.x1f.NewAthenaArray(ncells2);
          flx_outer_field_.x2f.NewAthenaArray(ncells2+1);
          flx_outer_field_.x3f.NewAthenaArray(ncells2);
          shboxvar_outer_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
          shboxvar_outer_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
          shboxmap_outer_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
          shboxmap_outer_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
          flx_outer_emf_.x2e.NewAthenaArray(ncells2);
          flx_outer_emf_.x3e.NewAthenaArray(ncells2+1);
        }
        shbb_.outer = true;
        shbb_.ogidlist=new int[nrbx2];
        shbb_.olidlist=new int[nrbx2];
        shbb_.ornklist=new int[nrbx2];
        shbb_.olevlist=new int[nrbx2];
        // attach corner cells from L/R side
        int size = (pmb->block_size.nx2+NGHOST)*ssize_*NHYDRO;
        int bsize=0, esize=0;
        if (MAGNETIC_FIELDS_ENABLED) {
          // extra cell in azimuth/vertical
          bsize = (pmb->block_size.nx2+NGHOST+1)*(ssize_+NGHOST)*NFIELD;
          // face plus edge for EMF
          esize = 2*(pmb->block_size.nx2+NGHOST)*pmb->block_size.nx3
                +pmb->block_size.nx2+pmb->block_size.nx3+NGHOST;
        }
        for (int n=0; n<2; n++) {
          send_outerbuf_hydro_[n] = new Real[size];
          recv_outerbuf_hydro_[n] = new Real[size];
          shbox_outer_hydro_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
          rq_outersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_outerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_outerbuf_field_[n] = new Real[bsize];
            recv_outerbuf_field_[n] = new Real[bsize];
            shbox_outer_field_flag_[n]=BNDRY_WAITING;
            send_outerbuf_emf_[n] = new Real[esize];
            recv_outerbuf_emf_[n] = new Real[esize];
            shbox_outer_emf_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
            rq_outersend_field_[n] = MPI_REQUEST_NULL;
            rq_outerrecv_field_[n] = MPI_REQUEST_NULL;
            rq_outersend_emf_[n] = MPI_REQUEST_NULL;
            rq_outerrecv_emf_[n] = MPI_REQUEST_NULL;
#endif
          }
        }
        size = NGHOST*ssize_*NHYDRO;// corner cells only
        if (MAGNETIC_FIELDS_ENABLED) {
          bsize = NGHOST*(ssize_+NGHOST)*NFIELD;
          esize = 2*NGHOST*pmb->block_size.nx3+NGHOST;
        }
        for (int n=2; n<4; n++) {
          send_outerbuf_hydro_[n] = new Real[size];
          recv_outerbuf_hydro_[n] = new Real[size];
          shbox_outer_hydro_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
          rq_outersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_outerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_outerbuf_field_[n] = new Real[bsize];
            recv_outerbuf_field_[n] = new Real[bsize];
            shbox_outer_field_flag_[n]=BNDRY_WAITING;
            send_outerbuf_emf_[n] = new Real[esize];
            recv_outerbuf_emf_[n] = new Real[esize];
            shbox_outer_emf_flag_[n]=BNDRY_WAITING;
#ifdef MPI_PARALLEL
            rq_outersend_field_[n] = MPI_REQUEST_NULL;
            rq_outerrecv_field_[n] = MPI_REQUEST_NULL;
            rq_outersend_emf_[n] = MPI_REQUEST_NULL;
            rq_outerrecv_emf_[n] = MPI_REQUEST_NULL;
#endif
          }
        }
      }
    }
  } // shearing box
}

// destructor

BoundaryValues::~BoundaryValues() {
  MeshBlock *pmb=pmy_block_;

  DestroyBoundaryData(bd_hydro_);
  if (pmy_mesh_->multilevel==true) // SMR or AMR
    DestroyBoundaryData(bd_flcor_);
  if (MAGNETIC_FIELDS_ENABLED) {
    DestroyBoundaryData(bd_field_);
    DestroyBoundaryData(bd_emfcor_);
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    if (num_north_polar_blocks_ > 0) {
      for (int n = 0; n < num_north_polar_blocks_; ++n) {
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
    if (num_south_polar_blocks_ > 0) {
      for (int n = 0; n < num_south_polar_blocks_; ++n) {
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
  if (pmb->loc.level == pmy_mesh_->root_level &&
     pmy_mesh_->nrbx3 == 1 &&
     (block_bcs[INNER_X2]==POLAR_BNDRY||block_bcs[OUTER_X2]==POLAR_BNDRY||
      block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE||block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE))
       exc_.DeleteAthenaArray();

  if (SHEARING_BOX) {
    int level = pmb->loc.level - pmb->pmy_mesh->root_level;
    int64_t nrbx1 = pmb->pmy_mesh->nrbx1*(1L << level);
    if (pmb->loc.lx1 == 0) { // if true for shearing inner blocks
      shboxvar_inner_hydro_.DeleteAthenaArray();
      flx_inner_hydro_.DeleteAthenaArray();
      for (int n=0; n<4; n++) {
        delete[] send_innerbuf_hydro_[n];
        delete[] recv_innerbuf_hydro_[n];
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        shboxvar_inner_field_.x1f.DeleteAthenaArray();
        shboxvar_inner_field_.x2f.DeleteAthenaArray();
        shboxvar_inner_field_.x3f.DeleteAthenaArray();
        flx_inner_field_.x1f.DeleteAthenaArray();
        flx_inner_field_.x2f.DeleteAthenaArray();
        flx_inner_field_.x3f.DeleteAthenaArray();
        shboxvar_inner_emf_.x2e.DeleteAthenaArray();
        shboxvar_inner_emf_.x3e.DeleteAthenaArray();
        flx_inner_emf_.x2e.DeleteAthenaArray();
        flx_inner_emf_.x3e.DeleteAthenaArray();
        for (int n=0; n<4; n++) {
          delete[] send_innerbuf_field_[n];
          delete[] recv_innerbuf_field_[n];
          delete[] send_innerbuf_emf_[n];
          delete[] recv_innerbuf_emf_[n];
        }
      }
    }
    if (pmb->loc.lx1 == (nrbx1-1)) { // if true for shearing outer blocks
      shboxvar_outer_hydro_.DeleteAthenaArray();
      flx_outer_hydro_.DeleteAthenaArray();
      for (int n=0; n<4; n++) {
        delete[] send_outerbuf_hydro_[n];
        delete[] recv_outerbuf_hydro_[n];
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        shboxvar_outer_field_.x1f.DeleteAthenaArray();
        shboxvar_outer_field_.x2f.DeleteAthenaArray();
        shboxvar_outer_field_.x3f.DeleteAthenaArray();
        flx_outer_field_.x1f.DeleteAthenaArray();
        flx_outer_field_.x2f.DeleteAthenaArray();
        flx_outer_field_.x3f.DeleteAthenaArray();
        shboxvar_outer_emf_.x2e.DeleteAthenaArray();
        shboxvar_outer_emf_.x3e.DeleteAthenaArray();
        flx_outer_emf_.x2e.DeleteAthenaArray();
        flx_outer_emf_.x3e.DeleteAthenaArray();
        for (int n=0; n<4; n++) {
          delete[] send_outerbuf_field_[n];
          delete[] recv_outerbuf_field_[n];
          delete[] send_outerbuf_emf_[n];
          delete[] recv_outerbuf_emf_[n];
        }
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::InitBoundaryData(BoundaryData &bd, enum BoundaryType type)
//  \brief Initialize BoundaryData structure
void BoundaryValues::InitBoundaryData(BoundaryData &bd, enum BoundaryType type) {
  MeshBlock *pmb=pmy_block_;
  bool multilevel=pmy_mesh_->multilevel;
  int f2d=0, f3d=0;
  int cng, cng1, cng2, cng3;
  if (pmb->block_size.nx2 > 1) f2d=1;
  if (pmb->block_size.nx3 > 1) f3d=1;
  cng=cng1=pmb->cnghost;
  cng2=cng*f2d;
  cng3=cng*f3d;
  int size;
  bd.nbmax=maxneighbor_;
  if (type==BNDRY_FLCOR || type==BNDRY_EMFCOR) {
    for (bd.nbmax=0; BoundaryValues::ni[bd.nbmax].type==NEIGHBOR_FACE; bd.nbmax++);
  }
  if (type==BNDRY_EMFCOR) {
    for (          ; BoundaryValues::ni[bd.nbmax].type==NEIGHBOR_EDGE; bd.nbmax++);
  }
  for (int n=0;n<bd.nbmax;n++) {
    // Clear flags and requests
    bd.flag[n]=BNDRY_WAITING;
    bd.send[n]=NULL;
    bd.recv[n]=NULL;
#ifdef MPI_PARALLEL
    bd.req_send[n]=MPI_REQUEST_NULL;
    bd.req_recv[n]=MPI_REQUEST_NULL;
#endif

    // Allocate buffers
    // calculate the buffer size
    switch(type) {
      case BNDRY_HYDRO: {
        size=((BoundaryValues::ni[n].ox1==0)?pmb->block_size.nx1:NGHOST)
            *((BoundaryValues::ni[n].ox2==0)?pmb->block_size.nx2:NGHOST)
            *((BoundaryValues::ni[n].ox3==0)?pmb->block_size.nx3:NGHOST);
        if (multilevel) {
          int f2c=((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                 *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                 *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
          int c2f=((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
                 *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
                 *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
          size=std::max(size,c2f);
          size=std::max(size,f2c);
        }
        size*=NHYDRO;
      }
      break;
      case BNDRY_FIELD: {
        int size1=((BoundaryValues::ni[n].ox1==0) ? (pmb->block_size.nx1+1):NGHOST)
                 *((BoundaryValues::ni[n].ox2==0) ? (pmb->block_size.nx2):NGHOST)
                 *((BoundaryValues::ni[n].ox3==0) ? (pmb->block_size.nx3):NGHOST);
        int size2=((BoundaryValues::ni[n].ox1==0) ? (pmb->block_size.nx1):NGHOST)
                 *((BoundaryValues::ni[n].ox2==0) ? (pmb->block_size.nx2+f2d):NGHOST)
                 *((BoundaryValues::ni[n].ox3==0) ? (pmb->block_size.nx3):NGHOST);
        int size3=((BoundaryValues::ni[n].ox1==0) ? (pmb->block_size.nx1):NGHOST)
                 *((BoundaryValues::ni[n].ox2==0) ? (pmb->block_size.nx2):NGHOST)
                 *((BoundaryValues::ni[n].ox3==0) ? (pmb->block_size.nx3+f3d):NGHOST);
        size=size1+size2+size3;
        if (multilevel) {
          if (BoundaryValues::ni[n].type!=NEIGHBOR_FACE) {
            if (BoundaryValues::ni[n].ox1!=0) size1=size1/NGHOST*(NGHOST+1);
            if (BoundaryValues::ni[n].ox2!=0) size2=size2/NGHOST*(NGHOST+1);
            if (BoundaryValues::ni[n].ox3!=0) size3=size3/NGHOST*(NGHOST+1);
          }
          size=size1+size2+size3;
          int f2c1=((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+1):NGHOST)
                  *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                  *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
          int f2c2=((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                  *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+f2d)
                    : NGHOST)
                  *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
          int f2c3=((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
                  *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
                  *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+f3d)
                    : NGHOST);
          if (BoundaryValues::ni[n].type!=NEIGHBOR_FACE) {
            if (BoundaryValues::ni[n].ox1!=0) f2c1=f2c1/NGHOST*(NGHOST+1);
            if (BoundaryValues::ni[n].ox2!=0) f2c2=f2c2/NGHOST*(NGHOST+1);
            if (BoundaryValues::ni[n].ox3!=0) f2c3=f2c3/NGHOST*(NGHOST+1);
          }
          int fsize=f2c1+f2c2+f2c3;
          int c2f1=
            ((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1+1):cng+1)
           *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
           *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
          int c2f2=
            ((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
           *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2+f2d):cng+1)
           *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng);
          int c2f3=
            ((BoundaryValues::ni[n].ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng)
           *((BoundaryValues::ni[n].ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng)
           *((BoundaryValues::ni[n].ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3+f3d):cng+1);
          int csize=c2f1+c2f2+c2f3;
          size=std::max(size,std::max(csize,fsize));
        }
      }
      break;
      case BNDRY_FLCOR: {
        if (BoundaryValues::ni[n].ox1!=0)
          size=(pmb->block_size.nx2+1)/2*(pmb->block_size.nx3+1)/2*NHYDRO;
        if (BoundaryValues::ni[n].ox2!=0)
          size=(pmb->block_size.nx1+1)/2*(pmb->block_size.nx3+1)/2*NHYDRO;
        if (BoundaryValues::ni[n].ox3!=0)
          size=(pmb->block_size.nx1+1)/2*(pmb->block_size.nx2+1)/2*NHYDRO;
      }
      break;
      case BNDRY_EMFCOR: {
        if (BoundaryValues::ni[n].type==NEIGHBOR_FACE) {
          if (pmb->block_size.nx3>1) { // 3D
            if (BoundaryValues::ni[n].ox1!=0)
              size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                  +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
            else if (BoundaryValues::ni[n].ox2!=0)
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                  +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
            else
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                  +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
          } else if (pmb->block_size.nx2>1) { // 2D
            if (BoundaryValues::ni[n].ox1!=0)
              size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
            else
              size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
          } else { // 1D
            size=2;
          }
        } else if (BoundaryValues::ni[n].type==NEIGHBOR_EDGE) {
          if (pmb->block_size.nx3>1) { // 3D
            if (BoundaryValues::ni[n].ox3==0) size=pmb->block_size.nx3;
            if (BoundaryValues::ni[n].ox2==0) size=pmb->block_size.nx2;
            if (BoundaryValues::ni[n].ox1==0) size=pmb->block_size.nx1;
          } else if (pmb->block_size.nx2>1) {
            size=1;
          }
        }
      }
      break;
      default: {
        std::stringstream msg;
        msg << "### FATAL ERROR in InitBoundaryData" << std::endl
            << "Invalid boundary type is specified." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      break;
    }
    bd.send[n]=new Real[size];
    bd.recv[n]=new Real[size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::DestroyBoundaryData(BoundaryData &bd)
//  \brief Destroy BoundaryData structure
void BoundaryValues::DestroyBoundaryData(BoundaryData &bd) {
  for (int n=0;n<bd.nbmax;n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
#ifdef MPI_PARALLEL
    if (bd.req_send[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_send[n]);
    if (bd.req_recv[n]!=MPI_REQUEST_NULL)
      MPI_Request_free(&bd.req_recv[n]);
#endif
  }
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::Initialize(void)
//  \brief Initialize MPI requests

void BoundaryValues::Initialize(void) {
  MeshBlock* pmb=pmy_block_;
  int myox1, myox2, myox3;
  int tag;
  int f2d=0, f3d=0;
  int cng, cng1, cng2, cng3;
  if (pmb->block_size.nx2 > 1) f2d=1;
  if (pmb->block_size.nx3 > 1) f3d=1;
  cng=cng1=pmb->cnghost;
  cng2=cng*f2d;
  cng3=cng*f3d;
  int ssize, rsize;
  int64_t &lx1=pmb->loc.lx1;
  int64_t &lx2=pmb->loc.lx2;
  int64_t &lx3=pmb->loc.lx3;
  int &mylevel=pmb->loc.level;
  myox1=(static_cast<int>(lx1&1L));
  myox2=(static_cast<int>(lx2&1L));
  myox3=(static_cast<int>(lx3&1L));

  // count the number of the fine meshblocks contacting on each edge
  int eid=0;
  if (pmb->block_size.nx2 > 1) {
    for (int ox2=-1;ox2<=1;ox2+=2) {
      for (int ox1=-1;ox1<=1;ox1+=2) {
        int nis, nie, njs, nje;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        int nf=0, fl=mylevel;
        for (int nj=njs; nj<=nje; nj++) {
          for (int ni=nis; ni<=nie; ni++) {
            if (nblevel[1][nj+1][ni+1] > fl)
              fl++, nf=0;
            if (nblevel[1][nj+1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }
  if (pmb->block_size.nx3 > 1) {
    for (int ox3=-1;ox3<=1;ox3+=2) {
      for (int ox1=-1;ox1<=1;ox1+=2) {
        int nis, nie, nks, nke;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for (int nk=nks; nk<=nke; nk++) {
          for (int ni=nis; ni<=nie; ni++) {
            if (nblevel[nk+1][1][ni+1] > fl)
              fl++, nf=0;
            if (nblevel[nk+1][1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
    for (int ox3=-1;ox3<=1;ox3+=2) {
      for (int ox2=-1;ox2<=1;ox2+=2) {
        int njs, nje, nks, nke;
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for (int nk=nks; nk<=nke; nk++) {
          for (int nj=njs; nj<=nje; nj++) {
            if (nblevel[nk+1][nj+1][1] > fl)
              fl++, nf=0;
            if (nblevel[nk+1][nj+1][1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }

#ifdef MPI_PARALLEL
  // Initialize non-polar neighbor communications to other ranks
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.rank!=Globals::my_rank) {
      if (nb.level==mylevel) { // same
        ssize=rsize=((nb.ox1==0)?pmb->block_size.nx1:NGHOST)
                   *((nb.ox2==0)?pmb->block_size.nx2:NGHOST)
                   *((nb.ox3==0)?pmb->block_size.nx3:NGHOST);
      } else if (nb.level<mylevel) { // coarser
        ssize=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
             *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
             *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
        rsize=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng1)
             *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng2)
             *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng3);
      } else { // finer
        ssize=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2+cng1):cng1)
             *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2+cng2):cng2)
             *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2+cng3):cng3);
        rsize=((nb.ox1==0) ? ((pmb->block_size.nx1+1)/2):NGHOST)
             *((nb.ox2==0) ? ((pmb->block_size.nx2+1)/2):NGHOST)
             *((nb.ox3==0) ? ((pmb->block_size.nx3+1)/2):NGHOST);
      }
      ssize*=NHYDRO; rsize*=NHYDRO;
      // specify the offsets in the view point of the target block: flip ox? signs
      tag=CreateBvalsMPITag(nb.lid, TAG_HYDRO, nb.targetid);
      if (bd_hydro_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
        MPI_Request_free(&bd_hydro_.req_send[nb.bufid]);
      MPI_Send_init(bd_hydro_.send[nb.bufid],ssize,MPI_ATHENA_REAL,
                    nb.rank,tag,MPI_COMM_WORLD,&(bd_hydro_.req_send[nb.bufid]));
      tag=CreateBvalsMPITag(pmb->lid, TAG_HYDRO, nb.bufid);
      if (bd_hydro_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
        MPI_Request_free(&bd_hydro_.req_recv[nb.bufid]);
      MPI_Recv_init(bd_hydro_.recv[nb.bufid],rsize,MPI_ATHENA_REAL,
                    nb.rank,tag,MPI_COMM_WORLD,&(bd_hydro_.req_recv[nb.bufid]));

      // flux correction
      if (pmy_mesh_->multilevel==true && nb.type==NEIGHBOR_FACE) {
        int size;
        if (nb.fid==0 || nb.fid==1)
          size=((pmb->block_size.nx2+1)/2)*((pmb->block_size.nx3+1)/2);
        else if (nb.fid==2 || nb.fid==3)
          size=((pmb->block_size.nx1+1)/2)*((pmb->block_size.nx3+1)/2);
        else // (nb.fid==4 || nb.fid==5)
          size=((pmb->block_size.nx1+1)/2)*((pmb->block_size.nx2+1)/2);
        size*=NHYDRO;
        if (nb.level<mylevel) { // send to coarser
          tag=CreateBvalsMPITag(nb.lid, TAG_HYDFLX, nb.targetid);
          if (bd_flcor_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
            MPI_Request_free(&bd_flcor_.req_send[nb.bufid]);
          MPI_Send_init(bd_flcor_.send[nb.bufid],size,MPI_ATHENA_REAL,
              nb.rank,tag,MPI_COMM_WORLD,&(bd_flcor_.req_send[nb.bufid]));
        } else if (nb.level>mylevel) { // receive from finer
          tag=CreateBvalsMPITag(pmb->lid, TAG_HYDFLX, nb.bufid);
          if (bd_flcor_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
            MPI_Request_free(&bd_flcor_.req_recv[nb.bufid]);
          MPI_Recv_init(bd_flcor_.recv[nb.bufid],size,MPI_ATHENA_REAL,
              nb.rank,tag,MPI_COMM_WORLD,&(bd_flcor_.req_recv[nb.bufid]));
        }
      }

      if (MAGNETIC_FIELDS_ENABLED) {
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
          if (nb.type!=NEIGHBOR_FACE) {
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
          if (nb.type!=NEIGHBOR_FACE) {
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
        }
        if (nb.level==mylevel) // same
          ssize=size, rsize=size;
        else if (nb.level<mylevel) // coarser
          ssize=fsize, rsize=csize;
        else // finer
          ssize=csize, rsize=fsize;

        tag=CreateBvalsMPITag(nb.lid, TAG_FIELD, nb.targetid);
        if (bd_field_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
          MPI_Request_free(&bd_field_.req_send[nb.bufid]);
        MPI_Send_init(bd_field_.send[nb.bufid],ssize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&(bd_field_.req_send[nb.bufid]));
        tag=CreateBvalsMPITag(pmb->lid, TAG_FIELD, nb.bufid);
        if (bd_field_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
          MPI_Request_free(&bd_field_.req_recv[nb.bufid]);
        MPI_Recv_init(bd_field_.recv[nb.bufid],rsize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&(bd_field_.req_recv[nb.bufid]));
        // EMF correction
        int fi1, fi2, f2csize;
        if (nb.type==NEIGHBOR_FACE) { // face
          if (pmb->block_size.nx3 > 1) { // 3D
            if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
              size=(pmb->block_size.nx2+1)*(pmb->block_size.nx3)
                  +(pmb->block_size.nx2)*(pmb->block_size.nx3+1);
              f2csize=(pmb->block_size.nx2/2+1)*(pmb->block_size.nx3/2)
                  +(pmb->block_size.nx2/2)*(pmb->block_size.nx3/2+1);
            } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx3)
                  +(pmb->block_size.nx1)*(pmb->block_size.nx3+1);
              f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx3/2)
                  +(pmb->block_size.nx1/2)*(pmb->block_size.nx3/2+1);
            } else if (nb.fid==INNER_X3 || nb.fid==OUTER_X3) {
              size=(pmb->block_size.nx1+1)*(pmb->block_size.nx2)
                  +(pmb->block_size.nx1)*(pmb->block_size.nx2+1);
              f2csize=(pmb->block_size.nx1/2+1)*(pmb->block_size.nx2/2)
                  +(pmb->block_size.nx1/2)*(pmb->block_size.nx2/2+1);
            }
          } else if (pmb->block_size.nx2 > 1) { // 2D
            if (nb.fid==INNER_X1 || nb.fid==OUTER_X1) {
              size=(pmb->block_size.nx2+1)+pmb->block_size.nx2;
              f2csize=(pmb->block_size.nx2/2+1)+pmb->block_size.nx2/2;
            } else if (nb.fid==INNER_X2 || nb.fid==OUTER_X2) {
              size=(pmb->block_size.nx1+1)+pmb->block_size.nx1;
              f2csize=(pmb->block_size.nx1/2+1)+pmb->block_size.nx1/2;
            }
          } else { // 1D
            size=f2csize=2;
          }
        } else if (nb.type==NEIGHBOR_EDGE) { // edge
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

        if (nb.level==mylevel) { // the same level
          if ((nb.type==NEIGHBOR_FACE) || ((nb.type==NEIGHBOR_EDGE)
                                           && (edge_flag_[nb.eid]==true))) {
            tag=CreateBvalsMPITag(nb.lid, TAG_FLDFLX, nb.targetid);
            if (bd_emfcor_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
              MPI_Request_free(&bd_emfcor_.req_send[nb.bufid]);
            MPI_Send_init(bd_emfcor_.send[nb.bufid],size,MPI_ATHENA_REAL,
                          nb.rank,tag,MPI_COMM_WORLD,&(bd_emfcor_.req_send[nb.bufid]));
            tag=CreateBvalsMPITag(pmb->lid, TAG_FLDFLX, nb.bufid);
            if (bd_emfcor_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
              MPI_Request_free(&bd_emfcor_.req_recv[nb.bufid]);
            MPI_Recv_init(bd_emfcor_.recv[nb.bufid],size,MPI_ATHENA_REAL,
                          nb.rank,tag,MPI_COMM_WORLD,&(bd_emfcor_.req_recv[nb.bufid]));
          }
        }
        if (nb.level>mylevel) { // finer neighbor
          tag=CreateBvalsMPITag(pmb->lid, TAG_FLDFLX, nb.bufid);
          if (bd_emfcor_.req_recv[nb.bufid]!=MPI_REQUEST_NULL)
            MPI_Request_free(&bd_emfcor_.req_recv[nb.bufid]);
          MPI_Recv_init(bd_emfcor_.recv[nb.bufid],f2csize,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&(bd_emfcor_.req_recv[nb.bufid]));
        }
        if (nb.level<mylevel) { // coarser neighbor
          tag=CreateBvalsMPITag(nb.lid, TAG_FLDFLX, nb.targetid);
          if (bd_emfcor_.req_send[nb.bufid]!=MPI_REQUEST_NULL)
            MPI_Request_free(&bd_emfcor_.req_send[nb.bufid]);
          MPI_Send_init(bd_emfcor_.send[nb.bufid],f2csize,MPI_ATHENA_REAL,
                        nb.rank,tag,MPI_COMM_WORLD,&(bd_emfcor_.req_send[nb.bufid]));
        }
      }
    }
  }

  // Initialize polar neighbor communications to other ranks
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int n = 0; n < num_north_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = polar_neighbor_north[n];
      if (nb.rank != Globals::my_rank) {
        tag = CreateBvalsMPITag(nb.lid, TAG_FLDFLX_POLE, pmb->loc.lx3);
        if (req_emf_north_send_[n]!=MPI_REQUEST_NULL)
          MPI_Request_free(&req_emf_north_send_[n]);
        MPI_Send_init(emf_north_send_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
            nb.rank, tag, MPI_COMM_WORLD, &req_emf_north_send_[n]);
        tag = CreateBvalsMPITag(pmb->lid, TAG_FLDFLX_POLE, n);
        if (req_emf_north_recv_[n]!=MPI_REQUEST_NULL)
          MPI_Request_free(&req_emf_north_recv_[n]);
        MPI_Recv_init(emf_north_recv_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
            nb.rank, tag, MPI_COMM_WORLD, &req_emf_north_recv_[n]);
      }
    }
    for (int n = 0; n < num_south_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = polar_neighbor_south[n];
      if (nb.rank != Globals::my_rank) {
        tag = CreateBvalsMPITag(nb.lid, TAG_FLDFLX_POLE, pmb->loc.lx3);
        if (req_emf_south_send_[n]!=MPI_REQUEST_NULL)
          MPI_Request_free(&req_emf_south_send_[n]);
        MPI_Send_init(emf_south_send_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
           nb.rank, tag, MPI_COMM_WORLD, &req_emf_south_send_[n]);
        tag = CreateBvalsMPITag(pmb->lid, TAG_FLDFLX_POLE, n);
        if (req_emf_south_recv_[n]!=MPI_REQUEST_NULL)
          MPI_Request_free(&req_emf_south_recv_[n]);
        MPI_Recv_init(emf_south_recv_[n], pmb->block_size.nx1, MPI_ATHENA_REAL,
            nb.rank, tag, MPI_COMM_WORLD, &req_emf_south_recv_[n]);
      }
    }
  }
#endif

// initialize the shearing block lists
  if (SHEARING_BOX) {
    Mesh *pmesh = pmb->pmy_mesh;
    int level = pmb->loc.level - pmesh->root_level;
    int64_t nrbx1 = pmesh->nrbx1*(1L << level);
    int64_t nrbx2 = pmesh->nrbx2*(1L << level);
    int nbtotal = pmesh->nbtotal;
    int *ranklist = pmesh->ranklist;
    int *nslist = pmesh->nslist;
    LogicalLocation *loclist = pmesh->loclist;

    int count = 0;
    if (shbb_.inner) {
      for (int i=0;i<nbtotal;i++) {
        if (loclist[i].lx1 == 0 && loclist[i].lx3 == pmb->loc.lx3 &&
            loclist[i].level == pmb->loc.level) {
          shbb_.igidlist[count] = i;
          shbb_.ilidlist[count] = i - nslist[ranklist[i]];
          shbb_.irnklist[count] = ranklist[i];
          shbb_.ilevlist[count] = loclist[i].level;
          count++;
        }
      }
    }
    count = 0;
    if (shbb_.outer) {
      for (int i=0;i<nbtotal;i++) {
        if (loclist[i].lx1 == (nrbx1-1) && loclist[i].lx3 == pmb->loc.lx3 &&
          loclist[i].level == pmb->loc.level) {
          shbb_.ogidlist[count] = i;
          shbb_.olidlist[count] = i - nslist[ranklist[i]];
          shbb_.ornklist[count] = ranklist[i];
          shbb_.olevlist[count] = loclist[i].level;
          count++;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckBoundary(void)
//  \brief checks if the boundary conditions are correctly enrolled

void BoundaryValues::CheckBoundary(void) {
  MeshBlock *pmb=pmy_block_;
  for (int i=0;i<nface_;i++) {
    if (block_bcs[i]==USER_BNDRY) {
      if (BoundaryFunction_[i]==NULL) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the hydro boundary function "
            << "is not enrolled in direction " << i  << "." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingForInit(bool cons_and_field)
//  \brief initiate MPI_Irecv for initialization

void BoundaryValues::StartReceivingForInit(bool cons_and_field) {
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.rank!=Globals::my_rank) {
      if (cons_and_field) {  // normal case
        MPI_Start(&(bd_hydro_.req_recv[nb.bufid]));
        if (MAGNETIC_FIELDS_ENABLED)
          MPI_Start(&(bd_field_.req_recv[nb.bufid]));
      } else { // must be primitive initialization
        MPI_Start(&(bd_hydro_.req_recv[nb.bufid]));
      }
    }
  }
#endif
// find send_block_id and recv_block_id;
  if (SHEARING_BOX) {
    MeshBlock *pmb=pmy_block_;
    Mesh *pmesh = pmb->pmy_mesh;
    FindShearBlock(pmesh->time);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingAll(const Real time)
//  \brief initiate MPI_Irecv for all the sweeps

void BoundaryValues::StartReceivingAll(const Real time) {
  firsttime_=true;
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.rank!=Globals::my_rank) {
      MPI_Start(&(bd_hydro_.req_recv[nb.bufid]));
      if (nb.type==NEIGHBOR_FACE && nb.level>mylevel)
        MPI_Start(&(bd_flcor_.req_recv[nb.bufid]));
      if (MAGNETIC_FIELDS_ENABLED) {
        MPI_Start(&(bd_field_.req_recv[nb.bufid]));
        if (nb.type==NEIGHBOR_FACE || nb.type==NEIGHBOR_EDGE) {
          if ((nb.level>mylevel) || ((nb.level==mylevel) && ((nb.type==NEIGHBOR_FACE)
          || ((nb.type==NEIGHBOR_EDGE) && (edge_flag_[nb.eid]==true)))))
           MPI_Start(&(bd_emfcor_.req_recv[nb.bufid]));
        }
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int n = 0; n < num_north_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = polar_neighbor_north[n];
      if (nb.rank != Globals::my_rank) {
        MPI_Start(&req_emf_north_recv_[n]);
      }
    }
    for (int n = 0; n < num_south_polar_blocks_; ++n) {
      const PolarNeighborBlock &nb = polar_neighbor_south[n];
      if (nb.rank != Globals::my_rank) {
        MPI_Start(&req_emf_south_recv_[n]);
      }
    }
  }
#endif
// find send_block_id and recv_block_id; post non-blocking recv
  if (SHEARING_BOX) {
    MeshBlock *pmb=pmy_block_;
    Mesh *pmesh = pmb->pmy_mesh;
    FindShearBlock(time);
#ifdef MPI_PARALLEL
    int size,tag;
    if (shbb_.inner) { // inner boundary
      for (int n=0; n<4; n++) {
        if ((recv_inner_rank_[n]!=Globals::my_rank) &&
                          (recv_inner_rank_[n]!=-1)) {
          size = ssize_*NHYDRO*recv_innersize_hydro_[n];
          tag  = CreateBvalsMPITag(pmb->lid, TAG_SHBOX_HYDRO, n);
          MPI_Irecv(recv_innerbuf_hydro_[n],size,MPI_ATHENA_REAL,
                    recv_inner_rank_[n],tag,MPI_COMM_WORLD,
                    &rq_innerrecv_hydro_[n]);
          if (MAGNETIC_FIELDS_ENABLED) {
            size = recv_innersize_field_[n];
            tag  = CreateBvalsMPITag(pmb->lid, TAG_SHBOX_FIELD, n);
            MPI_Irecv(recv_innerbuf_field_[n],size,MPI_ATHENA_REAL,
                      recv_inner_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_innerrecv_field_[n]);
            size = recv_innersize_emf_[n];
            tag  = CreateBvalsMPITag(pmb->lid, TAG_SHBOX_EMF, n);
            MPI_Irecv(recv_innerbuf_emf_[n],size,MPI_ATHENA_REAL,
                      recv_inner_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_innerrecv_emf_[n]);
          }
        }
    }}

    if (shbb_.outer) { // outer boundary
      int offset=4;
      for (int n=0; n<4; n++) {
        if ((recv_outer_rank_[n]!=Globals::my_rank) &&
                          (recv_outer_rank_[n]!=-1)) {
          size = ssize_*NHYDRO*recv_outersize_hydro_[n];
          tag  = CreateBvalsMPITag(pmb->lid, TAG_SHBOX_HYDRO, n+offset);
          MPI_Irecv(recv_outerbuf_hydro_[n],size,MPI_ATHENA_REAL,
                    recv_outer_rank_[n],tag,MPI_COMM_WORLD,
                    &rq_outerrecv_hydro_[n]);
          if (MAGNETIC_FIELDS_ENABLED) {
            size = recv_outersize_field_[n];
            tag  = CreateBvalsMPITag(pmb->lid, TAG_SHBOX_FIELD, n+offset);
            MPI_Irecv(recv_outerbuf_field_[n],size,MPI_ATHENA_REAL,
                      recv_outer_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_outerrecv_field_[n]);
            size = recv_outersize_emf_[n];
            tag  = CreateBvalsMPITag(pmb->lid, TAG_SHBOX_EMF, n+offset);
            MPI_Irecv(recv_outerbuf_emf_[n],size,MPI_ATHENA_REAL,
                      recv_outer_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_outerrecv_emf_[n]);
          }
        }
    }}
#endif
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryForInit(void)
//  \brief clean up the boundary flags for initialization

void BoundaryValues::ClearBoundaryForInit(bool cons_and_field) {
  MeshBlock *pmb=pmy_block_;

  // Note step==0 corresponds to initial exchange of conserved variables, while step==1
  // corresponds to primitives sent only in the case of GR with refinement
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    bd_hydro_.flag[nb.bufid] = BNDRY_WAITING;
    if (MAGNETIC_FIELDS_ENABLED)
      bd_field_.flag[nb.bufid] = BNDRY_WAITING;
    if (GENERAL_RELATIVITY and pmy_mesh_->multilevel)
      bd_hydro_.flag[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank) {
      if (cons_and_field) {  // normal case
        // Wait for Isend
        MPI_Wait(&(bd_hydro_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
        if (MAGNETIC_FIELDS_ENABLED)
          MPI_Wait(&(bd_field_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      } else {  // must be primitive initialization
        if (GENERAL_RELATIVITY and pmy_mesh_->multilevel)
          MPI_Wait(&(bd_hydro_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      }
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundaryAll(void)
//  \brief clean up the boundary flags after each loop

void BoundaryValues::ClearBoundaryAll(void) {
  MeshBlock *pmb=pmy_block_;

  // Clear non-polar boundary communications
  for (int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    bd_hydro_.flag[nb.bufid] = BNDRY_WAITING;
    if (nb.type==NEIGHBOR_FACE)
      bd_flcor_.flag[nb.bufid] = BNDRY_WAITING;
    if (MAGNETIC_FIELDS_ENABLED) {
      bd_field_.flag[nb.bufid] = BNDRY_WAITING;
      if ((nb.type==NEIGHBOR_FACE) || (nb.type==NEIGHBOR_EDGE))
        bd_emfcor_.flag[nb.bufid] = BNDRY_WAITING;
    }
#ifdef MPI_PARALLEL
    if (nb.rank!=Globals::my_rank) {
      // Wait for Isend
      MPI_Wait(&(bd_hydro_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      if (nb.type==NEIGHBOR_FACE && nb.level<pmb->loc.level)
        MPI_Wait(&(bd_flcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
      if (MAGNETIC_FIELDS_ENABLED) {
        MPI_Wait(&(bd_field_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
        if (nb.type==NEIGHBOR_FACE || nb.type==NEIGHBOR_EDGE) {
          if (nb.level < pmb->loc.level)
            MPI_Wait(&(bd_emfcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
          else if ((nb.level==pmb->loc.level) && ((nb.type==NEIGHBOR_FACE)
              || ((nb.type==NEIGHBOR_EDGE) && (edge_flag_[nb.eid]==true))))
            MPI_Wait(&(bd_emfcor_.req_send[nb.bufid]),MPI_STATUS_IGNORE);
        }
      }
    }
#endif
  }

  // Clear polar boundary communications
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int n = 0; n < num_north_polar_blocks_; ++n) {
      PolarNeighborBlock &nb = polar_neighbor_north[n];
      emf_north_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
      if (nb.rank != Globals::my_rank)
        MPI_Wait(&req_emf_north_send_[n], MPI_STATUS_IGNORE);
#endif
    }
    for (int n = 0; n < num_south_polar_blocks_; ++n) {
      PolarNeighborBlock &nb = polar_neighbor_south[n];
      emf_south_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
      if (nb.rank != Globals::my_rank)
        MPI_Wait(&req_emf_south_send_[n], MPI_STATUS_IGNORE);
#endif
    }
  }
// clear shearingbox boundary communications
  if (SHEARING_BOX) {
    if (shbb_.inner == true) {
      for (int n=0; n<4; n++) {
        if (send_inner_rank_[n] == -1) continue;
        shbox_inner_hydro_flag_[n] = BNDRY_WAITING;
        if (MAGNETIC_FIELDS_ENABLED) {
          shbox_inner_field_flag_[n] = BNDRY_WAITING;
          shbox_inner_emf_flag_[n] = BNDRY_WAITING;
        }
#ifdef MPI_PARALLEL
        if (send_inner_rank_[n]!=Globals::my_rank) {
          MPI_Wait(&rq_innersend_hydro_[n],MPI_STATUS_IGNORE);
          if (MAGNETIC_FIELDS_ENABLED) {
            MPI_Wait(&rq_innersend_field_[n],MPI_STATUS_IGNORE);
            MPI_Wait(&rq_innersend_emf_[n],MPI_STATUS_IGNORE);
          }
        }
#endif
      }
    } // inner boundary

    if (shbb_.outer == true) {
      for (int n=0; n<4; n++) {
        if (send_outer_rank_[n] == -1) continue;
        shbox_outer_hydro_flag_[n] = BNDRY_WAITING;
        if (MAGNETIC_FIELDS_ENABLED) {
          shbox_outer_field_flag_[n] = BNDRY_WAITING;
        }
#ifdef MPI_PARALLEL
        if (send_outer_rank_[n]!=Globals::my_rank) {
          Mesh *pmesh = pmb->pmy_mesh;
          MPI_Wait(&rq_outersend_hydro_[n],MPI_STATUS_IGNORE);
          if (MAGNETIC_FIELDS_ENABLED) {
            MPI_Wait(&rq_outersend_field_[n],MPI_STATUS_IGNORE);
            MPI_Wait(&rq_outersend_emf_[n],MPI_STATUS_IGNORE);
          }
        }
#endif
      }
    }
  } // end shearing box
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ApplyPhysicalBoundaries(AthenaArray<Real> &pdst,
//           AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst,
//           const Real time, const Real dt)
//  \brief Apply all the physical boundary conditions for both hydro and field

void BoundaryValues::ApplyPhysicalBoundaries(AthenaArray<Real> &pdst,
                     AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst,
                     const Real time, const Real dt) {
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco=pmb->pcoord;
  int bis=pmb->is-NGHOST, bie=pmb->ie+NGHOST, bjs=pmb->js, bje=pmb->je,
      bks=pmb->ks, bke=pmb->ke;
  if (BoundaryFunction_[INNER_X2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
  if (BoundaryFunction_[OUTER_X2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
  if (BoundaryFunction_[INNER_X3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
  if (BoundaryFunction_[OUTER_X3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;
  // Apply boundary function on inner-x1
  if (BoundaryFunction_[INNER_X1] != NULL) {
    BoundaryFunction_[INNER_X1](pmb, pco, pdst, bfdst, time, dt,
                                pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST);
    if (MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
        pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
    }
    pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
      pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
  }

  // Apply boundary function on outer-x1
  if (BoundaryFunction_[OUTER_X1] != NULL) {
    BoundaryFunction_[OUTER_X1](pmb, pco, pdst, bfdst, time, dt,
                                pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST);
    if (MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
        pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
    }
    pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
      pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
  }

  if (pmb->block_size.nx2>1) { // 2D or 3D

    // Apply boundary function on inner-x2
    if (BoundaryFunction_[INNER_X2] != NULL) {
      BoundaryFunction_[INNER_X2](pmb, pco, pdst, bfdst, time, dt,
                                  bis, bie, pmb->js, pmb->je, bks, bke, NGHOST);
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
      }
      pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
    }

    // Apply boundary function on outer-x2
    if (BoundaryFunction_[OUTER_X2] != NULL) {
      BoundaryFunction_[OUTER_X2](pmb, pco, pdst, bfdst, time, dt,
                                  bis, bie, pmb->js, pmb->je, bks, bke, NGHOST);
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
      }
      pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
    }
  }

  if (pmb->block_size.nx3>1) { // 3D
    bjs=pmb->js-NGHOST;
    bje=pmb->je+NGHOST;

    // Apply boundary function on inner-x3
    if (BoundaryFunction_[INNER_X3] != NULL) {
      BoundaryFunction_[INNER_X3](pmb, pco, pdst, bfdst, time, dt,
                                  bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST);
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
      }
      pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
    }

    // Apply boundary function on outer-x3
    if (BoundaryFunction_[OUTER_X3] != NULL) {
      BoundaryFunction_[OUTER_X3](pmb, pco, pdst, bfdst, time, dt,
                                  bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST);
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pco,
          bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
      }
      pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pco,
        bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ProlongateBoundaries(AthenaArray<Real> &pdst,
//           AthenaArray<Real> &cdst, FaceField &bdst, AthenaArray<Real> &bcdst,
//           const Real time, const Real dt)
//  \brief Prolongate the level boundary using the coarse data

void BoundaryValues::ProlongateBoundaries(AthenaArray<Real> &pdst,
     AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst,
     const Real time, const Real dt) {
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int64_t &lx1=pmb->loc.lx1;
  int64_t &lx2=pmb->loc.lx2;
  int64_t &lx3=pmb->loc.lx3;
  int &mylevel=pmb->loc.level;

  for (int n=0; n<nneighbor; n++) {
    NeighborBlock& nb = neighbor[n];
    if (nb.level >= mylevel) continue;
    // fill the required ghost-ghost zone
    int nis, nie, njs, nje, nks, nke;
    nis=std::max(nb.ox1-1,-1), nie=std::min(nb.ox1+1,1);
    if (pmb->block_size.nx2==1) njs=0, nje=0;
    else njs=std::max(nb.ox2-1,-1), nje=std::min(nb.ox2+1,1);
    if (pmb->block_size.nx3==1) nks=0, nke=0;
    else nks=std::max(nb.ox3-1,-1), nke=std::min(nb.ox3+1,1);
    for (int nk=nks; nk<=nke; nk++) {
      for (int nj=njs; nj<=nje; nj++) {
        for (int ni=nis; ni<=nie; ni++) {
          int ntype=std::abs(ni)+std::abs(nj)+std::abs(nk);
          // skip myself or coarse levels; only the same level must be restricted
          if (ntype==0 || nblevel[nk+1][nj+1][ni+1]!=mylevel) continue;

          // this neighbor block is on the same level
          // and needs to be restricted for prolongation
          int ris, rie, rjs, rje, rks, rke;
          if (ni==0) {
            ris=pmb->cis, rie=pmb->cie;
            if (nb.ox1==1) ris=pmb->cie;
            else if (nb.ox1==-1) rie=pmb->cis;
          } else if (ni== 1) {
            ris=pmb->cie+1, rie=pmb->cie+1;
          } else { //(ni==-1)
            ris=pmb->cis-1, rie=pmb->cis-1;
          }
          if (nj==0) {
            rjs=pmb->cjs, rje=pmb->cje;
            if (nb.ox2==1) rjs=pmb->cje;
            else if (nb.ox2==-1) rje=pmb->cjs;
          } else if (nj== 1) {
            rjs=pmb->cje+1, rje=pmb->cje+1;
          } else { //(nj==-1)
            rjs=pmb->cjs-1, rje=pmb->cjs-1;
          }
          if (nk==0) {
            rks=pmb->cks, rke=pmb->cke;
            if (nb.ox3==1) rks=pmb->cke;
            else if (nb.ox3==-1) rke=pmb->cks;
          } else if (nk== 1) {
            rks=pmb->cke+1, rke=pmb->cke+1;
          } else { //(nk==-1)
            rks=pmb->cks-1, rke=pmb->cks-1;
          }

          pmb->pmr->RestrictCellCenteredValues(cdst, pmr->coarse_cons_, 0, NHYDRO-1,
                                               ris, rie, rjs, rje, rks, rke);
          if (GENERAL_RELATIVITY)
            pmb->pmr->RestrictCellCenteredValues(pdst, pmr->coarse_prim_, 0, NHYDRO-1,
                                                 ris, rie, rjs, rje, rks, rke);
          if (MAGNETIC_FIELDS_ENABLED) {
            int rs=ris, re=rie+1;
            if (rs==pmb->cis   && nblevel[nk+1][nj+1][ni  ]<mylevel) rs++;
            if (re==pmb->cie+1 && nblevel[nk+1][nj+1][ni+2]<mylevel) re--;
            pmr->RestrictFieldX1(bfdst.x1f, pmr->coarse_b_.x1f, rs, re, rjs, rje, rks,
                                 rke);
            if (pmb->block_size.nx2 > 1) {
              rs=rjs, re=rje+1;
              if (rs==pmb->cjs   && nblevel[nk+1][nj  ][ni+1]<mylevel) rs++;
              if (re==pmb->cje+1 && nblevel[nk+1][nj+2][ni+1]<mylevel) re--;
              pmr->RestrictFieldX2(bfdst.x2f, pmr->coarse_b_.x2f, ris, rie, rs, re, rks,
                                   rke);
            } else { // 1D
              pmr->RestrictFieldX2(bfdst.x2f, pmr->coarse_b_.x2f, ris, rie, rjs, rje, rks,
                                   rke);
              for (int i=ris; i<=rie; i++)
                pmr->coarse_b_.x2f(rks,rjs+1,i)=pmr->coarse_b_.x2f(rks,rjs,i);
            }
            if (pmb->block_size.nx3 > 1) {
              rs=rks, re=rke+1;
              if (rs==pmb->cks   && nblevel[nk  ][nj+1][ni+1]<mylevel) rs++;
              if (re==pmb->cke+1 && nblevel[nk+2][nj+1][ni+1]<mylevel) re--;
              pmr->RestrictFieldX3(bfdst.x3f, pmr->coarse_b_.x3f, ris, rie, rjs, rje, rs,
                                   re);
            } else { // 1D or 2D
              pmr->RestrictFieldX3(bfdst.x3f, pmr->coarse_b_.x3f, ris, rie, rjs, rje, rks,
                                   rke);
              for (int j=rjs; j<=rje; j++) {
                for (int i=ris; i<=rie; i++)
                  pmr->coarse_b_.x3f(rks+1,j,i)=pmr->coarse_b_.x3f(rks,j,i);
              }
            }
          }
        }
      }
    }

    // calculate the loop limits for the ghost zones
    int cn = pmb->cnghost-1;
    int si, ei, sj, ej, sk, ek, fsi, fei, fsj, fej, fsk, fek;
    if (nb.ox1==0) {
      si=pmb->cis, ei=pmb->cie;
      if ((lx1&1L)==0L) ei+=cn;
      else             si-=cn;
    } else if (nb.ox1>0) { si=pmb->cie+1,  ei=pmb->cie+cn;}
    else              si=pmb->cis-cn, ei=pmb->cis-1;
    if (nb.ox2==0) {
      sj=pmb->cjs, ej=pmb->cje;
      if (pmb->block_size.nx2 > 1) {
        if ((lx2&1L)==0L) ej+=cn;
        else             sj-=cn;
      }
    } else if (nb.ox2>0) { sj=pmb->cje+1,  ej=pmb->cje+cn;}
    else              sj=pmb->cjs-cn, ej=pmb->cjs-1;
    if (nb.ox3==0) {
      sk=pmb->cks, ek=pmb->cke;
      if (pmb->block_size.nx3 > 1) {
        if ((lx3&1L)==0L) ek+=cn;
        else             sk-=cn;
      }
    } else if (nb.ox3>0) { sk=pmb->cke+1,  ek=pmb->cke+cn;}
    else              sk=pmb->cks-cn, ek=pmb->cks-1;

    // convert the ghost zone and ghost-ghost zones into primitive variables
    // this includes cell-centered field calculation
    int f1m=0, f1p=0, f2m=0, f2p=0, f3m=0, f3p=0;
    if (nb.ox1==0) {
      if (nblevel[1][1][0]!=-1) f1m=1;
      if (nblevel[1][1][2]!=-1) f1p=1;
    } else {
      f1m=1;
      f1p=1;
    }
    if (pmb->block_size.nx2>1) {
      if (nb.ox2==0) {
        if (nblevel[1][0][1]!=-1) f2m=1;
        if (nblevel[1][2][1]!=-1) f2p=1;
      } else {
        f2m=1;
        f2p=1;
      }
    }
    if (pmb->block_size.nx3>1) {
      if (nb.ox3==0) {
        if (nblevel[0][1][1]!=-1) f3m=1;
        if (nblevel[2][1][1]!=-1) f3p=1;
      } else {
        f3m=1;
        f3p=1;
      }
    }

    pmb->peos->ConservedToPrimitive(pmr->coarse_cons_, pmr->coarse_prim_,
                 pmr->coarse_b_, pmr->coarse_prim_, pmr->coarse_bcc_, pmr->pcoarsec,
                 si-f1m, ei+f1p, sj-f2m, ej+f2p, sk-f3m, ek+f3p);

    // Apply physical boundaries
    if (nb.ox1==0) {
      if (BoundaryFunction_[INNER_X1]!=NULL) {
        BoundaryFunction_[INNER_X1](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                pmr->coarse_b_, time, dt, pmb->cis, pmb->cie, sj, ej, sk, ek, 1);
      }
      if (BoundaryFunction_[OUTER_X1]!=NULL) {
        BoundaryFunction_[OUTER_X1](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                pmr->coarse_b_, time, dt, pmb->cis, pmb->cie, sj, ej, sk, ek, 1);
      }
    }
    if (nb.ox2==0 && pmb->block_size.nx2 > 1) {
      if (BoundaryFunction_[INNER_X2]!=NULL) {
        BoundaryFunction_[INNER_X2](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                pmr->coarse_b_, time, dt, si, ei, pmb->cjs, pmb->cje, sk, ek, 1);
      }
      if (BoundaryFunction_[OUTER_X2]!=NULL) {
        BoundaryFunction_[OUTER_X2](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                pmr->coarse_b_, time, dt, si, ei, pmb->cjs, pmb->cje, sk, ek, 1);
      }
    }
    if (nb.ox3==0 && pmb->block_size.nx3 > 1) {
      if (BoundaryFunction_[INNER_X3]!=NULL) {
        BoundaryFunction_[INNER_X3](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                pmr->coarse_b_, time, dt, si, ei, sj, ej, pmb->cks, pmb->cke, 1);
      }
      if (BoundaryFunction_[OUTER_X3]!=NULL) {
        BoundaryFunction_[OUTER_X3](pmb, pmr->pcoarsec, pmr->coarse_prim_,
                pmr->coarse_b_, time, dt, si, ei, sj, ej, pmb->cks, pmb->cke, 1);
      }
    }

    // now that the ghost-ghost zones are filled
    // calculate the loop limits for the finer grid
    fsi=(si-pmb->cis)*2+pmb->is,   fei=(ei-pmb->cis)*2+pmb->is+1;
    if (pmb->block_size.nx2 > 1)
      fsj=(sj-pmb->cjs)*2+pmb->js, fej=(ej-pmb->cjs)*2+pmb->js+1;
    else fsj=pmb->js, fej=pmb->je;
    if (pmb->block_size.nx3 > 1)
      fsk=(sk-pmb->cks)*2+pmb->ks, fek=(ek-pmb->cks)*2+pmb->ks+1;
    else fsk=pmb->ks, fek=pmb->ke;

    // prolongate hydro variables using primitive
    pmr->ProlongateCellCenteredValues(pmr->coarse_prim_, pdst, 0, NHYDRO-1,
                                      si, ei, sj, ej, sk, ek);
    // prolongate magnetic fields
    if (MAGNETIC_FIELDS_ENABLED) {
      int il, iu, jl, ju, kl, ku;
      il=si, iu=ei+1;
      if ((nb.ox1>=0) && (nblevel[nb.ox3+1][nb.ox2+1][nb.ox1  ]>=mylevel)) il++;
      if ((nb.ox1<=0) && (nblevel[nb.ox3+1][nb.ox2+1][nb.ox1+2]>=mylevel)) iu--;
      if (pmb->block_size.nx2 > 1) {
        jl=sj, ju=ej+1;
        if ((nb.ox2>=0) && (nblevel[nb.ox3+1][nb.ox2  ][nb.ox1+1]>=mylevel)) jl++;
        if ((nb.ox2<=0) && (nblevel[nb.ox3+1][nb.ox2+2][nb.ox1+1]>=mylevel)) ju--;
      } else {
        jl=sj;
        ju=ej;
      }
      if (pmb->block_size.nx3 > 1) {
        kl=sk, ku=ek+1;
        if ((nb.ox3>=0) && (nblevel[nb.ox3  ][nb.ox2+1][nb.ox1+1]>=mylevel)) kl++;
        if ((nb.ox3<=0) && (nblevel[nb.ox3+2][nb.ox2+1][nb.ox1+1]>=mylevel)) ku--;
      } else {
        kl=sk;
        ku=ek;
      }

      // step 1. calculate x1 outer surface fields and slopes
      pmr->ProlongateSharedFieldX1(pmr->coarse_b_.x1f, bfdst.x1f, il, iu, sj, ej, sk, ek);
      // step 2. calculate x2 outer surface fields and slopes
      pmr->ProlongateSharedFieldX2(pmr->coarse_b_.x2f, bfdst.x2f, si, ei, jl, ju, sk, ek);
      // step 3. calculate x3 outer surface fields and slopes
      pmr->ProlongateSharedFieldX3(pmr->coarse_b_.x3f, bfdst.x3f, si, ei, sj, ej, kl, ku);
      // step 4. calculate the internal finer fields using the Toth & Roe method
      pmr->ProlongateInternalField(bfdst, si, ei, sj, ej, sk, ek);
      // Field prolongation completed, calculate cell centered fields
      pmb->pfield->CalculateCellCenteredField(bfdst, bcdst, pmb->pcoord,
                                              fsi, fei, fsj, fej, fsk, fek);
    }
    // calculate conservative variables
    pmb->peos->PrimitiveToConserved(pdst, bcdst, cdst, pmb->pcoord,
                                    fsi, fei, fsj, fej, fsk, fek);
  }
}
