//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.cpp
//  \brief constructor/destructor and utility functions for BoundaryValues class

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // std::memcpy
#include <iomanip>
#include <iostream>   // endl
#include <iterator>
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"
#include "bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// BoundaryValues constructor (the first object constructed inside the MeshBlock()
// constructor): sets functions for the appropriate boundary conditions at each of the 6
// dirs of a MeshBlock
BoundaryValues::BoundaryValues(MeshBlock *pmb, BoundaryFlag *input_bcs,
                               ParameterInput *pin)
    : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs) {
  pmy_block_ = pmb;
  // Check BC functions for each of the 6 boundaries in turn ---------------------
  for (int i=0; i<6; i++) {
    switch(block_bcs[i]) {
      case BoundaryFlag::reflect:
      case BoundaryFlag::outflow:
      case BoundaryFlag::user:
      case BoundaryFlag::polar_wedge:
        apply_bndry_fn_[i] = true;
        break;
      default: // already initialized to false in class
        break;
    }
  }
  // Inner x1
  nface_ = 2; nedge_ = 0;
  CheckBoundaryFlag(block_bcs[BoundaryFace::inner_x1], CoordinateDirection::X1DIR);
  CheckBoundaryFlag(block_bcs[BoundaryFace::outer_x1], CoordinateDirection::X1DIR);

  if (pmb->block_size.nx2 > 1) {
    nface_ = 4; nedge_ = 4;
    CheckBoundaryFlag(block_bcs[BoundaryFace::inner_x2], CoordinateDirection::X2DIR);
    CheckBoundaryFlag(block_bcs[BoundaryFace::outer_x2], CoordinateDirection::X2DIR);
  }

  if (pmb->block_size.nx3 > 1) {
    nface_=6; nedge_=12;
    CheckBoundaryFlag(block_bcs[BoundaryFace::inner_x3], CoordinateDirection::X3DIR);
    CheckBoundaryFlag(block_bcs[BoundaryFace::outer_x3], CoordinateDirection::X3DIR);
  }
  // Perform compatibilty checks of user selections of polar vs. polar_wedge boundaries
  if (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
      || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
      || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge) {
    CheckPolarBoundaries();
  }

  // polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range
  if ((pmb->loc.level == pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1)
      && (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
       || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
       || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
       || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
    azimuthal_shift_.NewAthenaArray(pmb->ke + NGHOST + 2);

  // prevent reallocation of contiguous memory space for each of 4x possible calls to
  // std::vector<BoundaryVariable *>.push_back() in Hydro, Field, Gravity, PassiveScalars
  bvars.reserve(3);
  // TOOD(KGF): rename to "bvars_time_int"? What about a std::vector for bvars_sts?
  bvars_main_int.reserve(2);
  // reserve phys=0 for former TAG_AMR=8; now hard-coded in Mesh::CreateAMRMPITag()
  bvars_next_phys_id_ = 1;

  // KGF: BVals constructor section only containing ALL shearing box-specific stuff
  // set parameters for shearing box bc and allocate buffers
  if (SHEARING_BOX) {
    Omega_0_ = pin->GetOrAddReal("problem","Omega0",0.001);
    qshear_  = pin->GetOrAddReal("problem","qshear",1.5);
    ShBoxCoord_ = pin->GetOrAddInteger("problem","shboxcoord",1);
    x1size_ = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
    x2size_ = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
    x3size_ = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
    int level = pmb->loc.level - pmy_mesh_->root_level;
    std::int64_t nrbx1 = pmy_mesh_->nrbx1*(1L << level);
    std::int64_t nrbx2 = pmy_mesh_->nrbx2*(1L << level);
    shbb_.outer = false;
    shbb_.inner = false;

    if (ShBoxCoord_ == 1) {
      int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
      int ncells3 = pmb->block_size.nx3;
      if (pmy_mesh_->mesh_size.nx3 > 1) ncells3 += 2*NGHOST;
      ssize_ = NGHOST*ncells3;

      if (pmb->loc.lx1 == 0) { // if true for shearing inner blocks
        if (block_bcs[BoundaryFace::inner_x1] != BoundaryFlag::shear_periodic) {
          block_bcs[BoundaryFace::inner_x1] = BoundaryFlag::shear_periodic;
          BoundaryFunction_[BoundaryFace::inner_x1] = nullptr;
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
        shbb_.igidlist = new int[nrbx2];
        shbb_.ilidlist = new int[nrbx2];
        shbb_.irnklist = new int[nrbx2];
        shbb_.ilevlist = new int[nrbx2];
        // attach corner cells from L/R side
        int size = (pmb->block_size.nx2 + NGHOST)*ssize_*NHYDRO;
        int bsize = 0, esize = 0;
        if (MAGNETIC_FIELDS_ENABLED) {
          // extra cell in azimuth/vertical
          bsize = (pmb->block_size.nx2 + NGHOST+1)*(ssize_ + NGHOST)*NFIELD;
          // face plus edge for EMF
          esize = 2*(pmb->block_size.nx2 + NGHOST)*pmb->block_size.nx3
                +pmb->block_size.nx2+pmb->block_size.nx3 + NGHOST;
        }
        for (int n=0; n<2; n++) {
          send_innerbuf_hydro_[n] = new Real[size];
          recv_innerbuf_hydro_[n] = new Real[size];
          shbox_inner_hydro_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
          rq_innersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_innerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_innerbuf_field_[n] = new Real[bsize];
            recv_innerbuf_field_[n] = new Real[bsize];
            shbox_inner_field_flag_[n] = BoundaryStatus::waiting;
            send_innerbuf_emf_[n] = new Real[esize];
            recv_innerbuf_emf_[n] = new Real[esize];
            shbox_inner_emf_flag_[n] = BoundaryStatus::waiting;
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
            bsize = NGHOST*(ssize_ + NGHOST)*NFIELD;
            esize = 2*NGHOST*pmb->block_size.nx3 + NGHOST;
        }
        for (int n=2; n<4; n++) {
          send_innerbuf_hydro_[n] = new Real[size];
          recv_innerbuf_hydro_[n] = new Real[size];
          shbox_inner_hydro_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
          rq_innersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_innerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_innerbuf_field_[n] = new Real[bsize];
            recv_innerbuf_field_[n] = new Real[bsize];
            shbox_inner_field_flag_[n] = BoundaryStatus::waiting;
            send_innerbuf_emf_[n] = new Real[esize];
            recv_innerbuf_emf_[n] = new Real[esize];
            shbox_inner_emf_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
            rq_innersend_field_[n] = MPI_REQUEST_NULL;
            rq_innerrecv_field_[n] = MPI_REQUEST_NULL;
            rq_innersend_emf_[n] = MPI_REQUEST_NULL;
            rq_innerrecv_emf_[n] = MPI_REQUEST_NULL;
#endif
          }
        }
      }

      if (pmb->loc.lx1 == (nrbx1 - 1)) { // if true for shearing outer blocks
        if (block_bcs[BoundaryFace::outer_x1] != BoundaryFlag::shear_periodic) {
          block_bcs[BoundaryFace::outer_x1] = BoundaryFlag::shear_periodic;
          BoundaryFunction_[BoundaryFace::outer_x1] = nullptr;
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
        shbb_.ogidlist = new int[nrbx2];
        shbb_.olidlist = new int[nrbx2];
        shbb_.ornklist = new int[nrbx2];
        shbb_.olevlist = new int[nrbx2];
        // attach corner cells from L/R side
        int size = (pmb->block_size.nx2 + NGHOST)*ssize_*NHYDRO;
        int bsize = 0, esize = 0;
        if (MAGNETIC_FIELDS_ENABLED) {
          // extra cell in azimuth/vertical
          bsize = (pmb->block_size.nx2 + NGHOST + 1)*(ssize_ + NGHOST)*NFIELD;
          // face plus edge for EMF
          esize = 2*(pmb->block_size.nx2 + NGHOST)*pmb->block_size.nx3
                + pmb->block_size.nx2+pmb->block_size.nx3 + NGHOST;
        }
        for (int n=0; n<2; n++) {
          send_outerbuf_hydro_[n] = new Real[size];
          recv_outerbuf_hydro_[n] = new Real[size];
          shbox_outer_hydro_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
          rq_outersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_outerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_outerbuf_field_[n] = new Real[bsize];
            recv_outerbuf_field_[n] = new Real[bsize];
            shbox_outer_field_flag_[n] = BoundaryStatus::waiting;
            send_outerbuf_emf_[n] = new Real[esize];
            recv_outerbuf_emf_[n] = new Real[esize];
            shbox_outer_emf_flag_[n] = BoundaryStatus::waiting;
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
          bsize = NGHOST*(ssize_ + NGHOST)*NFIELD;
          esize = 2*NGHOST*pmb->block_size.nx3 + NGHOST;
        }
        for (int n=2; n<4; n++) {
          send_outerbuf_hydro_[n] = new Real[size];
          recv_outerbuf_hydro_[n] = new Real[size];
          shbox_outer_hydro_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
          rq_outersend_hydro_[n] = MPI_REQUEST_NULL;
          rq_outerrecv_hydro_[n] = MPI_REQUEST_NULL;
#endif
          if (MAGNETIC_FIELDS_ENABLED) {
            send_outerbuf_field_[n] = new Real[bsize];
            recv_outerbuf_field_[n] = new Real[bsize];
            shbox_outer_field_flag_[n] = BoundaryStatus::waiting;
            send_outerbuf_emf_[n] = new Real[esize];
            recv_outerbuf_emf_[n] = new Real[esize];
            shbox_outer_emf_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
            rq_outersend_field_[n] = MPI_REQUEST_NULL;
            rq_outerrecv_field_[n] = MPI_REQUEST_NULL;
            rq_outersend_emf_[n] = MPI_REQUEST_NULL;
            rq_outerrecv_emf_[n] = MPI_REQUEST_NULL;
#endif
          }
        }
      }
    } // end KGF: if (ShBoxCoord_ == 1)
  } // end KGF: shearing box in BoundaryValues constructor
}

// destructor

BoundaryValues::~BoundaryValues() {
  MeshBlock *pmb = pmy_block_;

  // edge-case of single block across pole in MHD spherical polar coordinates
  if ((pmb->loc.level == pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1)
      && (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
          || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
          || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
          || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
    azimuthal_shift_.DeleteAthenaArray();

  // KGF: shearing box destructor
  if (SHEARING_BOX) {
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
    if (pmb->loc.lx1 == (nrbx1 - 1)) { // if true for shearing outer blocks
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
  } // KGF: end shearing box handling in destructor
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetupPersistentMPI()
//  \brief Setup persistent MPI requests to be reused throughout the entire simulation

void BoundaryValues::SetupPersistentMPI() {
  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    (*bvars_it)->SetupPersistentMPI();
  }

  // KGF: begin exclusive shearing-box section in BoundaryValues::SetupPersistentMPI()
  // initialize the shearing block lists
  if (SHEARING_BOX) {
    Mesh *pmesh = pmb->pmy_mesh;
    int level = pmb->loc.level - pmesh->root_level;
    std::int64_t nrbx1 = pmesh->nrbx1*(1L << level);
    // std::int64_t nrbx2 = pmesh->nrbx2*(1L << level); // unused variable
    int nbtotal = pmesh->nbtotal;
    int *ranklist = pmesh->ranklist;
    int *nslist = pmesh->nslist;
    LogicalLocation *loclist = pmesh->loclist;

    int count = 0;
    if (shbb_.inner) {
      for (int i=0; i<nbtotal; i++) {
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
      for (int i=0; i<nbtotal; i++) {
        if (loclist[i].lx1 == (nrbx1 - 1) && loclist[i].lx3 == pmb->loc.lx3 &&
            loclist[i].level == pmb->loc.level) {
          shbb_.ogidlist[count] = i;
          shbb_.olidlist[count] = i - nslist[ranklist[i]];
          shbb_.ornklist[count] = ranklist[i];
          shbb_.olevlist[count] = loclist[i].level;
          count++;
        }
      }
    }
  } // end KGF: exclusive shearing box portion of SetupPersistentMPI()
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckUserBoundaries()
//  \brief checks if the boundary functions are correctly enrolled (this compatibility
//  check is performed at the top of Mesh::Initialize(), after calling ProblemGenerator())

void BoundaryValues::CheckUserBoundaries() {
  for (int i=0; i<nface_; i++) {
    if (block_bcs[i] == BoundaryFlag::user) {
      if (pmy_mesh_->BoundaryFunction_[i] == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the actual boundary function "
            << "is not enrolled in direction " << i  << " (in [0,6])." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceiving(BoundaryCommSubset phase)
//  \brief initiate MPI_Irecv()

void BoundaryValues::StartReceiving(BoundaryCommSubset phase) {
  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    (*bvars_it)->StartReceiving(phase);
  }

  // KGF: begin shearing-box exclusive section of original StartReceivingForInit()
  // find send_block_id and recv_block_id;
  if (SHEARING_BOX) {
    MeshBlock *pmb = pmy_block_;
    Mesh *pmesh = pmb->pmy_mesh;
    FindShearBlock(pmesh->time);
  }
  // end KGF: shearing box

  // KGF: begin shearing-box exclusive section of original StartReceivingAll()
  // find send_block_id and recv_block_id; post non-blocking recv
  if (SHEARING_BOX) {
    FindShearBlock(time);
#ifdef MPI_PARALLEL
     //     Mesh *pmesh = pmb->pmy_mesh;
    int size,tag;
    if (shbb_.inner) { // inner boundary
      for (int n=0; n<4; n++) {
        if ((recv_inner_rank_[n]!=Globals::my_rank) &&
                          (recv_inner_rank_[n]!=-1)) {
          size = ssize_*NHYDRO*recv_innersize_hydro_[n];
          tag  = CreateBvalsMPITag(pmb->lid, n, AthenaTagMPI::shbox_hydro);
          MPI_Irecv(recv_innerbuf_hydro_[n],size,MPI_ATHENA_REAL,
                    recv_inner_rank_[n],tag,MPI_COMM_WORLD,
                    &rq_innerrecv_hydro_[n]);
          if (MAGNETIC_FIELDS_ENABLED) {
            size = recv_innersize_field_[n];
            tag  = CreateBvalsMPITag(pmb->lid, n, AthenaTagMPI::shbox_field);
            MPI_Irecv(recv_innerbuf_field_[n],size,MPI_ATHENA_REAL,
                      recv_inner_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_innerrecv_field_[n]);
            size = recv_innersize_emf_[n];
            tag  = CreateBvalsMPITag(pmb->lid, n, AthenaTagMPI::shbox_emf);
            MPI_Irecv(recv_innerbuf_emf_[n],size,MPI_ATHENA_REAL,
                      recv_inner_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_innerrecv_emf_[n]);
          }
        }
      }
    }

    if (shbb_.outer) { // outer boundary
      int offset=4;
      for (int n=0; n<4; n++) {
        if ((recv_outer_rank_[n]!=Globals::my_rank) &&
                          (recv_outer_rank_[n]!=-1)) {
          size = ssize_*NHYDRO*recv_outersize_hydro_[n];
          tag  = CreateBvalsMPITag(pmb->lid, n+offset, AthenaTagMPI::shbox_hydro);
          MPI_Irecv(recv_outerbuf_hydro_[n],size,MPI_ATHENA_REAL,
                    recv_outer_rank_[n],tag,MPI_COMM_WORLD,
                    &rq_outerrecv_hydro_[n]);
          if (MAGNETIC_FIELDS_ENABLED) {
            size = recv_outersize_field_[n];
            tag  = CreateBvalsMPITag(pmb->lid, n+offset, AthenaTagMPI::shbox_field);
            MPI_Irecv(recv_outerbuf_field_[n],size,MPI_ATHENA_REAL,
                      recv_outer_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_outerrecv_field_[n]);
            size = recv_outersize_emf_[n];
            tag  = CreateBvalsMPITag(pmb->lid, n+offset, AthenaTagMPI::shbox_emf);
            MPI_Irecv(recv_outerbuf_emf_[n],size,MPI_ATHENA_REAL,
                      recv_outer_rank_[n],tag,MPI_COMM_WORLD,
                      &rq_outerrecv_emf_[n]);
          }
        }
      }
    }
#endif
  } // end KGF: shearing-box exclusive section of StartReceivingAll
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundary(BoundaryCommSubset phase)
//  \brief clean up the boundary flags after each loop

void BoundaryValues::ClearBoundary(BoundaryCommSubset phase) {
  // Note BoundaryCommSubset::mesh_init corresponds to initial exchange of conserved fluid
  // variables and magentic fields, while BoundaryCommSubset::gr_amr corresponds to fluid
  // primitive variables sent only in the case of GR with refinement
  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    (*bvars_it)->ClearBoundary(phase);
  }

  // KGF: begin shearing box exclusive section of ClearBoundaryAll
  // clear shearing box boundary communications
  if (SHEARING_BOX) {
    if (shbb_.inner == true) {
      for (int n=0; n<4; n++) {
        if (send_inner_rank_[n] == -1) continue;
        shbox_inner_hydro_flag_[n] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          shbox_inner_field_flag_[n] = BoundaryStatus::waiting;
          shbox_inner_emf_flag_[n] = BoundaryStatus::waiting;
        }
#ifdef MPI_PARALLEL
        if (send_inner_rank_[n] != Globals::my_rank) {
          MPI_Wait(&rq_innersend_hydro_[n], MPI_STATUS_IGNORE);
          if (MAGNETIC_FIELDS_ENABLED) {
            MPI_Wait(&rq_innersend_field_[n], MPI_STATUS_IGNORE);
            MPI_Wait(&rq_innersend_emf_[n], MPI_STATUS_IGNORE);
          }
        }
#endif
      }
    } // inner boundary

    if (shbb_.outer == true) {
      for (int n=0; n<4; n++) {
        if (send_outer_rank_[n] == -1) continue;
        shbox_outer_hydro_flag_[n] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          shbox_outer_field_flag_[n] = BoundaryStatus::waiting;
        }
#ifdef MPI_PARALLEL
        if (send_outer_rank_[n] != Globals::my_rank) {
          MPI_Wait(&rq_outersend_hydro_[n], MPI_STATUS_IGNORE);
          if (MAGNETIC_FIELDS_ENABLED) {
            MPI_Wait(&rq_outersend_field_[n], MPI_STATUS_IGNORE);
            MPI_Wait(&rq_outersend_emf_[n], MPI_STATUS_IGNORE);
          }
        }
#endif
      }
    }
  } // end KGF: shearing box
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ApplyPhysicalBoundaries(const Real time, const Real dt)
//  \brief Apply all the physical boundary conditions for both hydro and field

void BoundaryValues::ApplyPhysicalBoundaries(const Real time, const Real dt) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int bis = pmb->is - NGHOST, bie = pmb->ie + NGHOST,
      bjs = pmb->js, bje = pmb->je,
      bks = pmb->ks, bke = pmb->ke;

  // Extend the transverse limits that correspond to periodic boundaries as they are
  // updated: x1, then x2, then x3
  if (!apply_bndry_fn_[BoundaryFace::inner_x2] && pmb->block_size.nx2 > 1)
    bjs = pmb->js - NGHOST;
  if (!apply_bndry_fn_[BoundaryFace::outer_x2] && pmb->block_size.nx2 > 1)
    bje = pmb->je + NGHOST;
  if (!apply_bndry_fn_[BoundaryFace::inner_x3] && pmb->block_size.nx3 > 1)
    bks = pmb->ks - NGHOST;
  if (!apply_bndry_fn_[BoundaryFace::outer_x3] && pmb->block_size.nx3 > 1)
    bke = pmb->ke + NGHOST;

  // KGF: temporarily hardcode Hydro and Field access for coupling in EOS U(W) + calc bcc
  // and when passed to user-defined boundary function stored in function pointer array

  // downcast BoundaryVariable ptrs to known derived class types: RTTI via dynamic_cast
  HydroBoundaryVariable *phbvar =
      dynamic_cast<HydroBoundaryVariable *>(bvars_main_int[0]);
  Hydro *ph = pmb->phydro;

  FaceCenteredBoundaryVariable *pfbvar = nullptr;
  Field *pf = nullptr;
  if (MAGNETIC_FIELDS_ENABLED) {
    pf = pmb->pfield;
    pfbvar = dynamic_cast<FaceCenteredBoundaryVariable *>(bvars_main_int[1]);
  }

  // Apply boundary function on inner-x1 and update W,bcc (if not periodic)
  if (apply_bndry_fn_[BoundaryFace::inner_x1]) {
    DispatchBoundaryFunctions(pmb, pco, time, dt,
                              pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST,
                              ph->w, pf->b, BoundaryFace::inner_x1);
    // KGF: COUPLING OF QUANTITIES (must be manually specified)
    if (MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                              pmb->is-NGHOST, pmb->is-1,
                                              bjs, bje, bks, bke);
    }
    pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                    pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
  }

  // Apply boundary function on outer-x1 and update W,bcc (if not periodic)
  if (apply_bndry_fn_[BoundaryFace::outer_x1]) {
    DispatchBoundaryFunctions(pmb, pco, time, dt,
                              pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST,
                              ph->w, pf->b, BoundaryFace::outer_x1);
    // KGF: COUPLING OF QUANTITIES (must be manually specified)
    if (MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                              pmb->ie+1, pmb->ie+NGHOST,
                                              bjs, bje, bks, bke);
    }
    pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                    pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
  }

  if (pmb->block_size.nx2 > 1) { // 2D or 3D
    // Apply boundary function on inner-x2 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::inner_x2]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, pmb->js, pmb->je, bks, bke, NGHOST,
                                ph->w, pf->b, BoundaryFace::inner_x2);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, pmb->js-NGHOST, pmb->js-1,
                                                bks, bke);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
    }

    // Apply boundary function on outer-x2 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::outer_x2]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, pmb->js, pmb->je, bks, bke, NGHOST,
                                ph->w, pf->b, BoundaryFace::outer_x2);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, pmb->je+1, pmb->je+NGHOST,
                                                bks, bke);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
    }
  }

  if (pmb->block_size.nx3 > 1) { // 3D
    bjs = pmb->js - NGHOST;
    bje = pmb->je + NGHOST;

    // Apply boundary function on inner-x3 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::inner_x3]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST,
                                ph->w, pf->b, BoundaryFace::inner_x3);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, bjs, bje,
                                                pmb->ks-NGHOST, pmb->ks-1);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
    }

    // Apply boundary function on outer-x3 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::outer_x3]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST,
                                ph->w, pf->b, BoundaryFace::outer_x3);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, bjs, bje,
                                                pmb->ke+1, pmb->ke+NGHOST);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
    }
  }
  return;
}

// Public function, allows calling sites other than BoundaryVariable objects
// E.g. if chemistry or radiation elects to communicate additional information with MPI
// outside the framework of the BoundaryVariable classes

int BoundaryValues::ReserveTagVariableIDs(int num_phys) {
  // TODO(felker): add safety checks? input, output are positive, obey <= 31= MAX_NUM_PHYS
  int start_id = bvars_next_phys_id_;
  bvars_next_phys_id_ += num_phys;
  return start_id;
}


// KGF: should "bvars_it" be fixed in this class member function? Or passed as argument?
void BoundaryValues::DispatchBoundaryFunctions(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh,
    AthenaArray<Real> &prim, FaceField &b, BoundaryFace face) {
  if (block_bcs[face] ==  BoundaryFlag::user) {  // user-enrolled BCs
    pmy_mesh_->BoundaryFunction_[face](pmb, pco, prim, b, time, dt,
                                       il, iu, jl, ju, kl, ku, NGHOST);
  }
  // KGF: this is only to silence the compiler -Wswitch warnings about not handling the
  // "undef" case when considering all possible BoundaryFace enumerator values. If "undef"
  // is actually passed to this function, it will likely die before that ATHENA_ERROR()
  // call, as it tries to access block_bcs[-1]
  std::stringstream msg;
  msg << "### FATAL ERROR in DispatchBoundaryFunctions" << std::endl
      << "face = BoundaryFace::undef passed to this function" << std::endl;

  // For any function in the BoundaryPhysics interface class, iterate over
  // BoundaryVariable pointers "enrolled"
  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    switch(block_bcs[face]) {
      case BoundaryFlag::user: // handled above, outside loop over BoundaryVariable objs
        break;
      case BoundaryFlag::reflect:
        switch(face) {
          case BoundaryFace::undef:
            ATHENA_ERROR(msg);
          case BoundaryFace::inner_x1:
            (*bvars_it)->ReflectInnerX1(time, dt, il, jl, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x1:
            (*bvars_it)->ReflectOuterX1(time, dt, iu, jl, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::inner_x2:
            (*bvars_it)->ReflectInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x2:
            (*bvars_it)->ReflectOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::inner_x3:
            (*bvars_it)->ReflectInnerX3(time, dt, il, iu, jl, ju, kl, NGHOST);
            break;
          case BoundaryFace::outer_x3:
            (*bvars_it)->ReflectOuterX3(time, dt, il, iu, jl, ju, ku, NGHOST);
            break;
        }
        break;
      case BoundaryFlag::outflow:
        switch(face) {
          case BoundaryFace::undef:
            ATHENA_ERROR(msg);
          case BoundaryFace::inner_x1:
            (*bvars_it)->OutflowInnerX1(time, dt, il, jl, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x1:
            (*bvars_it)->OutflowOuterX1(time, dt, iu, jl, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::inner_x2:
            (*bvars_it)->OutflowInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x2:
            (*bvars_it)->OutflowOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::inner_x3:
            (*bvars_it)->OutflowInnerX3(time, dt, il, iu, jl, ju, kl, NGHOST);
            break;
          case BoundaryFace::outer_x3:
            (*bvars_it)->OutflowOuterX3(time, dt, il, iu, jl, ju, ku, NGHOST);
            break;
        }
        break;
      case BoundaryFlag::polar_wedge:
        switch(face) {
          case BoundaryFace::undef:
            ATHENA_ERROR(msg);
          case BoundaryFace::inner_x2:
            (*bvars_it)->PolarWedgeInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x2:
            (*bvars_it)->PolarWedgeOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
            break;
          default:
            std::stringstream msg_polar;
            msg_polar << "### FATAL ERROR in DispatchBoundaryFunctions" << std::endl
                << "Attempting to call polar wedge boundary function on \n"
                << "MeshBlock boundary other than inner x2 or outer x2" << std::endl;
            ATHENA_ERROR(msg_polar);
        }
        break;
      default:
        std::stringstream msg_flag;
        msg_flag << "### FATAL ERROR in DispatchBoundaryFunctions" << std::endl
                 << "No BoundaryPhysics function associated with provided\n"
                 << "block_bcs[" << face << "] = BoundaryFlag::"
                 << GetBoundaryString(block_bcs[face]) << std::endl;
        ATHENA_ERROR(msg);
        break;
    } // end switch(block_bcs[face])
  } // end loop over BoundaryVariable *
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::FindShearBlock(const Real time)
//  \brief Calculate the following quantities:
//  send_gid recv_gid send_lid recv_lid send_rank recv_rank,
//  send_size_hydro  recv_size_hydro: for MPI_Irecv
//  eps_,joverlap_: for update the conservative

// TODO(felker): break up this ~400 line function:

void BoundaryValues::FindShearBlock(const Real time) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  Mesh *pmesh = pmb->pmy_mesh;

  int js = pmb->js; int je = pmb->je;

  int level = pmb->loc.level - pmesh->root_level;
  std::int64_t nrbx2 = pmesh->nrbx2*(1L << level);
  int nx2   = pmb->block_size.nx2;  // # of cells per meshblock
  int nx3   = pmb->block_size.nx3;  // # of cells per meshblock
  // KGF: for symmetry reasons, how can ncells3 but not ncells2 be used in this fn?
  // int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
  int ncells3 = pmb->block_size.nx3;
  if (pmesh->mesh_size.nx3 > 1) ncells3 += 2*NGHOST;

  Real qomL = qshear_*Omega_0_*x1size_;
  Real yshear = qomL*time;
  Real deltay = std::fmod(yshear, x2size_);
  int joffset = static_cast<int>(deltay/pco->dx2v(js)); // assumes uniform grid in azimuth
  int Ngrids  = static_cast<int>(joffset/nx2);
  joverlap_   = joffset - Ngrids*nx2;
  eps_ = (std::fmod(deltay, pco->dx2v(js)))/pco->dx2v(js);

  if (shbb_.inner == true) { // if inner block
    for (int n=0; n<4; n++) {
      send_inner_gid_[n]  = -1;
      send_inner_rank_[n] = -1;
      send_inner_lid_[n]  = -1;
      recv_inner_gid_[n]  = -1;
      recv_inner_rank_[n] = -1;
      recv_inner_lid_[n]  = -1;
      send_innersize_cc_[n] = 0;
      recv_innersize_cc_[n] = 0;
      shbox_inner_cc_flag_[n] = BoundaryStatus::completed;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_fc_[n] = 0;
        recv_innersize_fc_[n] = 0;
        shbox_inner_fc_flag_[n] = BoundaryStatus::completed;
        send_innersize_fc_flx_[n] = 0;
        recv_innersize_fc_flx_[n] = 0;
        shbox_inner_fc_flx_flag_[n] = BoundaryStatus::completed;
      }
    }
    int jblock = 0;
    for (int j=0; j<nrbx2; j++) {
      // index of current meshblock on the shearingboundary block list
      if (shbb_.igidlist[j] == pmb->gid)  jblock = j;
    }
    // send [js:je-joverlap] of the meshblock to other
    // attach [je-joverlap+1:MIN(je-joverlap+(NGHOST), je-js+1)]
    // to its right end.
    std::int64_t jtmp = jblock + Ngrids;
    if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    send_inner_gid_[1]  = shbb_.igidlist[jtmp];
    send_inner_rank_[1] = shbb_.irnklist[jtmp];
    send_inner_lid_[1]  = shbb_.ilidlist[jtmp];
    send_innersize_cc_[1] = std::min(je -js-joverlap_ + 1+NGHOST, je -js + 1);
    // recv [js+joverlap:je] from other
    // attach [je+1:MIN(je+NGHOST, je+joverlap)] to its right end.
    jtmp = jblock - Ngrids;
    if (jtmp < 0) jtmp += nrbx2;
    recv_inner_gid_[1]  = shbb_.igidlist[jtmp];
    recv_inner_rank_[1] = shbb_.irnklist[jtmp];
    recv_inner_lid_[1]  = shbb_.ilidlist[jtmp];
    recv_innersize_cc_[1] = send_innersize_cc_[1];
    shbox_inner_cc_flag_[1] = BoundaryStatus::waiting;
    if (MAGNETIC_FIELDS_ENABLED) {
      send_innersize_fc_[1] = send_innersize_cc_[1]
                                 *NGHOST*(NFIELD*ncells3 + 1)
                                 +NGHOST*ncells3;
      recv_innersize_fc_[1] = send_innersize_fc_[1];
      shbox_inner_fc_flag_[1] = BoundaryStatus::waiting;
      send_innersize_fc_flx_[1] = send_innersize_cc_[1]*(2*nx3 + 1) + nx3;
      recv_innersize_fc_flx_[1] = send_innersize_fc_flx_[1];
      shbox_inner_fc_flx_flag_[1] = BoundaryStatus::waiting;
    }


    // if there is overlap to next blocks
    if (joverlap_ != 0) {
      // send to the right
      jtmp = jblock + (Ngrids + 1);
      if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[0]  = shbb_.igidlist[jtmp];
      send_inner_rank_[0] = shbb_.irnklist[jtmp];
      send_inner_lid_[0]  = shbb_.ilidlist[jtmp];
      send_innersize_cc_[0] = std::min(joverlap_+NGHOST, je -js + 1);
      // receive from its left
      jtmp = jblock - (Ngrids + 1);
      if (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[0]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[0] = shbb_.irnklist[jtmp];
      recv_inner_lid_[0]  = shbb_.ilidlist[jtmp];
      recv_innersize_cc_[0] = send_innersize_cc_[0];
      shbox_inner_cc_flag_[0] = BoundaryStatus::waiting;// switch on if overlap
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_fc_[0] = send_innersize_cc_[0]
                                   *NGHOST*(NFIELD*ncells3 + 1)
                                   +NGHOST*ncells3;
        recv_innersize_fc_[0] = send_innersize_fc_[0];
        shbox_inner_fc_flag_[0] = BoundaryStatus::waiting;
        send_innersize_fc_flx_[0] = send_innersize_cc_[0]*(2*nx3 + 1)+nx3;
        recv_innersize_fc_flx_[0] = send_innersize_fc_flx_[0];
        shbox_inner_fc_flx_flag_[0] = BoundaryStatus::waiting;
      }
      // deal the left boundary cells with send[2]
      if (joverlap_ > (nx2 - NGHOST)) {
        // send to Right
        jtmp = jblock + (Ngrids + 2);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        send_inner_gid_[2]  = shbb_.igidlist[jtmp];
        send_inner_rank_[2] = shbb_.irnklist[jtmp];
        send_inner_lid_[2]  = shbb_.ilidlist[jtmp];
        send_innersize_cc_[2] = joverlap_-(nx2-NGHOST);
        // recv from Left
        jtmp = jblock - (Ngrids + 2);
        while (jtmp < 0) jtmp += nrbx2;
        recv_inner_gid_[2]  = shbb_.igidlist[jtmp];
        recv_inner_rank_[2] = shbb_.irnklist[jtmp];
        recv_inner_lid_[2]  = shbb_.ilidlist[jtmp];
        recv_innersize_cc_[2] = send_innersize_cc_[2];
        shbox_inner_cc_flag_[2] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_fc_[2] = send_innersize_cc_[2]*NGHOST*(NFIELD*ncells3 + 1);
          recv_innersize_fc_[2] = send_innersize_fc_[2];
          shbox_inner_fc_flag_[2] = BoundaryStatus::waiting;
          send_innersize_fc_flx_[2] = send_innersize_cc_[2]*(2*nx3 + 1);
          recv_innersize_fc_flx_[2] = send_innersize_fc_flx_[2];
          shbox_inner_fc_flx_flag_[2] = BoundaryStatus::waiting;
        }
      }
      // deal with the right boundary cells with send[3]
      if (joverlap_ < NGHOST) {
        jtmp = jblock + (Ngrids - 1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        send_inner_gid_[3]  = shbb_.igidlist[jtmp];
        send_inner_rank_[3] = shbb_.irnklist[jtmp];
        send_inner_lid_[3]  = shbb_.ilidlist[jtmp];
        send_innersize_cc_[3] = NGHOST-joverlap_;
        jtmp = jblock - (Ngrids-1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        recv_inner_gid_[3]  = shbb_.igidlist[jtmp];
        recv_inner_rank_[3] = shbb_.irnklist[jtmp];
        recv_inner_lid_[3]  = shbb_.ilidlist[jtmp];
        recv_innersize_cc_[3] = send_innersize_cc_[3];
        shbox_inner_cc_flag_[3] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_innersize_fc_[3] = send_innersize_cc_[3]*NGHOST
                                     *(NFIELD*ncells3 + 1);
          recv_innersize_fc_[3] = send_innersize_fc_[3];
          shbox_inner_fc_flag_[3] = BoundaryStatus::waiting;
          send_innersize_fc_flx_[3] = send_innersize_cc_[3]*(2*nx3 + 1);
          recv_innersize_fc_flx_[3] = send_innersize_fc_flx_[3];
          shbox_inner_fc_flx_flag_[3] = BoundaryStatus::waiting;
        }
      }
    } else { // joverlap_ == 0
      // send [je-(NGHOST-1):je] to Right
      jtmp = jblock + (Ngrids + 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      send_inner_gid_[2]  = shbb_.igidlist[jtmp];
      send_inner_rank_[2] = shbb_.irnklist[jtmp];
      send_inner_lid_[2]  = shbb_.ilidlist[jtmp];
      send_innersize_cc_[2] = NGHOST;
      // recv [js-NGHOST:js-1] from Left
      jtmp = jblock - (Ngrids + 1);
      while (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[2]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[2] = shbb_.irnklist[jtmp];
      recv_inner_lid_[2]  = shbb_.ilidlist[jtmp];
      recv_innersize_cc_[2] = send_innersize_cc_[2];
      shbox_inner_cc_flag_[2] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_fc_[2] = send_innersize_cc_[2]*NGHOST
                                   *(NFIELD*ncells3 + 1);
        recv_innersize_fc_[2] = send_innersize_fc_[2];
        shbox_inner_fc_flag_[2] = BoundaryStatus::waiting;
        send_innersize_fc_flx_[2] = send_innersize_cc_[2]*(2*nx3 + 1);
        recv_innersize_fc_flx_[2] = send_innersize_fc_flx_[2];
        shbox_inner_fc_flx_flag_[2] = BoundaryStatus::waiting;
      }

      // send [js:js+(NGHOST-1)] to Left
      jtmp = jblock + (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      send_inner_gid_[3]  = shbb_.igidlist[jtmp];
      send_inner_rank_[3] = shbb_.irnklist[jtmp];
      send_inner_lid_[3]  = shbb_.ilidlist[jtmp];
      send_innersize_cc_[3] = NGHOST;
      // recv [je + 1:je+(NGHOST-1)] from Right
      jtmp = jblock - (Ngrids-1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      recv_inner_gid_[3]  = shbb_.igidlist[jtmp];
      recv_inner_rank_[3] = shbb_.irnklist[jtmp];
      recv_inner_lid_[3]  = shbb_.ilidlist[jtmp];
      recv_innersize_cc_[3] = send_innersize_cc_[3];
      shbox_inner_cc_flag_[3] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_innersize_fc_[3] = send_innersize_cc_[3]*NGHOST
                                   *(NFIELD*ncells3 + 1);
        recv_innersize_fc_[3] = send_innersize_fc_[3];
        shbox_inner_fc_flag_[3] = BoundaryStatus::waiting;
        send_innersize_fc_flx_[3] = send_innersize_cc_[3]*(2*nx3 + 1);
        recv_innersize_fc_flx_[3] = send_innersize_fc_flx_[3];
        shbox_inner_fc_flx_flag_[3] = BoundaryStatus::waiting;
      }
    }
  } // inner bc

  if (shbb_.outer == true) { // if outer block
    for (int n=0; n<4; n++) {
      send_outer_gid_[n]  = -1;
      send_outer_rank_[n] = -1;
      send_outer_lid_[n]  = -1;
      recv_outer_gid_[n]  = -1;
      recv_outer_rank_[n] = -1;
      recv_outer_lid_[n]  = -1;
      send_outersize_cc_[n] = 0;
      recv_outersize_cc_[n] = 0;
      shbox_outer_cc_flag_[n] = BoundaryStatus::completed;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_fc_[n] = 0;
        recv_outersize_fc_[n] = 0;
        shbox_outer_fc_flag_[n] = BoundaryStatus::completed;
        send_outersize_fc_flx_[n] = 0;
        recv_outersize_fc_flx_[n] = 0;
        shbox_outer_fc_flx_flag_[n] = BoundaryStatus::completed;
      }
    }
    int jblock = 0;
    for (int j=0; j<nrbx2; j++) {
      // index of current meshblock on the shearingboundary block list
      if (shbb_.ogidlist[j] == pmb->gid) jblock = j;
    }
    // recv [js-NGHOST:je-joverlap] of the meshblock from other
    std::int64_t jtmp = jblock + Ngrids;
    if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
    recv_outer_gid_[1]  = shbb_.ogidlist[jtmp];
    recv_outer_rank_[1] = shbb_.ornklist[jtmp];
    recv_outer_lid_[1]  = shbb_.olidlist[jtmp];
    recv_outersize_cc_[1] = std::min(je -js - joverlap_ + 1 + NGHOST, je -js + 1);
    // send [js+joverlap-NGHOST:je] of the meshblock to other
    jtmp = jblock - Ngrids;
    if (jtmp < 0) jtmp += nrbx2;
    send_outer_gid_[1]  = shbb_.ogidlist[jtmp];
    send_outer_rank_[1] = shbb_.ornklist[jtmp];
    send_outer_lid_[1]  = shbb_.olidlist[jtmp];
    send_outersize_cc_[1] = recv_outersize_cc_[1];
    shbox_outer_cc_flag_[1] = BoundaryStatus::waiting;
    if (MAGNETIC_FIELDS_ENABLED) {
      send_outersize_fc_[1] = send_outersize_cc_[1]*NGHOST*(NFIELD*ncells3 + 1) +
                              NGHOST*ncells3;
      recv_outersize_fc_[1] = send_outersize_fc_[1];
      shbox_outer_fc_flag_[1] = BoundaryStatus::waiting;
      send_outersize_fc_flx_[1] = send_outersize_cc_[1]*(2*nx3 + 1)+nx3;
      recv_outersize_fc_flx_[1] = send_outersize_fc_flx_[1];
      shbox_outer_fc_flx_flag_[1] = BoundaryStatus::waiting;
    }

    // if there is overlap to next blocks
    if (joverlap_ != 0) {
      // recv from right
      jtmp = jblock + (Ngrids + 1);
      if (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[0]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[0] = shbb_.ornklist[jtmp];
      recv_outer_lid_[0]  = shbb_.olidlist[jtmp];
      recv_outersize_cc_[0] = std::min(joverlap_+NGHOST, je -js + 1);
      // send to left
      jtmp = jblock - (Ngrids + 1);
      if (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[0]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[0] = shbb_.ornklist[jtmp];
      send_outer_lid_[0]  = shbb_.olidlist[jtmp];
      send_outersize_cc_[0] = recv_outersize_cc_[0];
      shbox_outer_cc_flag_[0] = BoundaryStatus::waiting; // switch on if overlap
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_fc_[0] = send_outersize_cc_[0]
                                   *NGHOST*(NFIELD*ncells3 + 1)
                                   +NGHOST*ncells3;
        recv_outersize_fc_[0] = send_outersize_fc_[0];
        shbox_outer_fc_flag_[0] = BoundaryStatus::waiting;
        send_outersize_fc_flx_[0] = send_outersize_cc_[0]*(2*nx3 + 1)+nx3;
        recv_outersize_fc_flx_[0] = send_outersize_fc_flx_[0];
        shbox_outer_fc_flx_flag_[0] = BoundaryStatus::waiting;
      }
      // deal the left boundary cells with send[2]
      if (joverlap_ > (nx2 - NGHOST)) {
        // recv from left
        jtmp = jblock + (Ngrids + 2);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        recv_outer_gid_[2]  = shbb_.ogidlist[jtmp];
        recv_outer_rank_[2] = shbb_.ornklist[jtmp];
        recv_outer_lid_[2]  = shbb_.olidlist[jtmp];
        recv_outersize_cc_[2] = joverlap_-(nx2-NGHOST);
        // send to right
        jtmp = jblock - (Ngrids + 2);
        while (jtmp < 0) jtmp += nrbx2;
        send_outer_gid_[2]  = shbb_.ogidlist[jtmp];
        send_outer_rank_[2] = shbb_.ornklist[jtmp];
        send_outer_lid_[2]  = shbb_.olidlist[jtmp];
        send_outersize_cc_[2] = recv_outersize_cc_[2];
        shbox_outer_cc_flag_[2] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_fc_[2] = send_outersize_cc_[2]*NGHOST*(NFIELD*ncells3 + 1);
          recv_outersize_fc_[2] = send_outersize_fc_[2];
          shbox_outer_fc_flag_[2] = BoundaryStatus::waiting;
          send_outersize_fc_flx_[2] = send_outersize_cc_[2]*(2*nx3 + 1);
          recv_outersize_fc_flx_[2] = send_outersize_fc_flx_[2];
          shbox_outer_fc_flx_flag_[2] = BoundaryStatus::waiting;
        }
      }
      // deal the right boundary cells with send[3]
      if (joverlap_ < NGHOST) {
        jtmp = jblock + (Ngrids - 1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        recv_outer_gid_[3]  = shbb_.ogidlist[jtmp];
        recv_outer_rank_[3] = shbb_.ornklist[jtmp];
        recv_outer_lid_[3]  = shbb_.olidlist[jtmp];
        recv_outersize_cc_[3] = NGHOST-joverlap_;
        jtmp = jblock - (Ngrids-1);
        while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
        while (jtmp < 0) jtmp += nrbx2;
        send_outer_gid_[3]  = shbb_.ogidlist[jtmp];
        send_outer_rank_[3] = shbb_.ornklist[jtmp];
        send_outer_lid_[3]  = shbb_.olidlist[jtmp];
        send_outersize_cc_[3] = recv_outersize_cc_[3];
        shbox_outer_cc_flag_[3] = BoundaryStatus::waiting;
        if (MAGNETIC_FIELDS_ENABLED) {
          send_outersize_fc_[3] = send_outersize_cc_[3]*NGHOST
                                     *(NFIELD*ncells3 + 1);
          recv_outersize_fc_[3] = send_outersize_fc_[3];
          shbox_outer_fc_flag_[3] = BoundaryStatus::waiting;
          send_outersize_fc_flx_[3] = send_outersize_cc_[3]*(2*nx3 + 1);
          recv_outersize_fc_flx_[3] = send_outersize_fc_flx_[3];
          shbox_outer_fc_flx_flag_[3] = BoundaryStatus::waiting;
        }
      }
    } else { // joverlap_ == 0
      // recv [je + 1:je+NGHOST] from Left
      jtmp = jblock + (Ngrids + 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      recv_outer_gid_[2]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[2] = shbb_.ornklist[jtmp];
      recv_outer_lid_[2]  = shbb_.olidlist[jtmp];
      recv_outersize_cc_[2] = NGHOST;
      // send [js:js+NGHOST-1] to Right
      jtmp = jblock - (Ngrids + 1);
      while (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[2]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[2] = shbb_.ornklist[jtmp];
      send_outer_lid_[2]  = shbb_.olidlist[jtmp];
      send_outersize_cc_[2] = NGHOST;
      shbox_outer_cc_flag_[2] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_fc_[2] = send_outersize_cc_[2]
                                   *NGHOST*(NFIELD*ncells3 + 1);
        recv_outersize_fc_[2] = send_outersize_fc_[2];
        shbox_outer_fc_flag_[2] = BoundaryStatus::waiting;
        send_outersize_fc_flx_[2] = send_outersize_cc_[2]*(2*nx3 + 1);
        recv_outersize_fc_flx_[2] = send_outersize_fc_flx_[2];
        shbox_outer_fc_flx_flag_[2] = BoundaryStatus::waiting;
      }

      // recv [js-NGHOST:js-1] from Left
      jtmp = jblock + (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      recv_outer_gid_[3]  = shbb_.ogidlist[jtmp];
      recv_outer_rank_[3] = shbb_.ornklist[jtmp];
      recv_outer_lid_[3]  = shbb_.olidlist[jtmp];
      recv_outersize_cc_[3] = NGHOST;
      // send [je-(NGHOST-1):je] to Right
      jtmp = jblock - (Ngrids - 1);
      while (jtmp > (nrbx2 - 1)) jtmp -= nrbx2;
      while (jtmp < 0) jtmp += nrbx2;
      send_outer_gid_[3]  = shbb_.ogidlist[jtmp];
      send_outer_rank_[3] = shbb_.ornklist[jtmp];
      send_outer_lid_[3]  = shbb_.olidlist[jtmp];
      send_outersize_cc_[3] = NGHOST;
      shbox_outer_cc_flag_[3] = BoundaryStatus::waiting;
      if (MAGNETIC_FIELDS_ENABLED) {
        send_outersize_fc_[3] = send_outersize_cc_[3]
                                   *NGHOST*(NFIELD*ncells3 + 1);
        recv_outersize_fc_[3] = send_outersize_fc_[3];
        shbox_outer_fc_flag_[3] = BoundaryStatus::waiting;
        send_outersize_fc_flx_[3] = send_outersize_cc_[3]*(2*nx3 + 1);
        recv_outersize_fc_flx_[3] = send_outersize_fc_flx_[3];
        shbox_outer_fc_flx_flag_[3] = BoundaryStatus::waiting;
      }
    }
  }
  return;
}
