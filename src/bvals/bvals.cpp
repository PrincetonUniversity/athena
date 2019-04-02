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

namespace {
void BValFuncPlaceholder(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt,
    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // Free function temporarily used to replace all BValFunc function pointer targets other
  // than user-defined boundary functions
  return;
}
} // namespace

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

// called in MeshBlock() constructor

BoundaryValues::BoundaryValues(MeshBlock *pmb, BoundaryFlag *input_bcs,
                               ParameterInput *pin)
    : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs) {
  pmy_block_ = pmb;
  for (int i=0; i<6; i++)
    BoundaryFunction_[i] = nullptr;

  // Set BC functions for each of the 6 boundaries in turn ------------------------------
  // Inner x1
  nface_ = 2; nedge_ = 0;
  switch(block_bcs[BoundaryFace::inner_x1]) {
    case BoundaryFlag::reflect:
      BoundaryFunction_[BoundaryFace::inner_x1] = BValFuncPlaceholder;
      break;
    case BoundaryFlag::outflow:
      BoundaryFunction_[BoundaryFace::inner_x1] = BValFuncPlaceholder;
      break;
    case BoundaryFlag::block: // block boundary
    case BoundaryFlag::periodic: // periodic boundary
      BoundaryFunction_[BoundaryFace::inner_x1] = nullptr;
      break;
    // case BoundaryFlag::shear_periodic: // shearing periodic boundary
    //   if (!SHEARING_BOX) block_bcs[BoundaryFace::inner_x1]=BoundaryFlag::periodic;
    //   BoundaryFunction_[BoundaryFace::inner_x1] = nullptr;
    //   break;
    case BoundaryFlag::user: // user-enrolled BCs
      BoundaryFunction_[BoundaryFace::inner_x1] =
          pmy_mesh_->BoundaryFunction_[BoundaryFace::inner_x1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ix1_bc=" <<  static_cast<int>(block_bcs[BoundaryFace::inner_x1])
          << " not valid" << std::endl;
      ATHENA_ERROR(msg);
      break;
  }

  // Outer x1
  switch(block_bcs[BoundaryFace::outer_x1]) {
    case BoundaryFlag::reflect:
      BoundaryFunction_[BoundaryFace::outer_x1] = BValFuncPlaceholder;
      break;
    case BoundaryFlag::outflow:
      BoundaryFunction_[BoundaryFace::outer_x1] = BValFuncPlaceholder;
      break;
    case BoundaryFlag::block: // block boundary
    case BoundaryFlag::periodic: // periodic boundary
      BoundaryFunction_[BoundaryFace::outer_x1] = nullptr;
      break;
    // case BoundaryFlag::shear_periodic: // shearing periodic boundary
    //   if (!SHEARING_BOX) block_bcs[BoundaryFace::outer_x1]=BoundaryFlag::periodic;
    //   BoundaryFunction_[BoundaryFace::outer_x1] = nullptr;
    //   break;
    case BoundaryFlag::user: // user-enrolled BCs
      BoundaryFunction_[BoundaryFace::outer_x1] =
          pmy_mesh_->BoundaryFunction_[BoundaryFace::outer_x1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ox1_bc=" <<  static_cast<int>(block_bcs[BoundaryFace::outer_x1])
          << " not valid" << std::endl;
      ATHENA_ERROR(msg);
  }

  if (pmb->block_size.nx2 > 1) {
    nface_=4; nedge_=4;
    // Inner x2
    switch(block_bcs[BoundaryFace::inner_x2]) {
      case BoundaryFlag::reflect:
        BoundaryFunction_[BoundaryFace::inner_x2] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::outflow:
        BoundaryFunction_[BoundaryFace::inner_x2] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::block: // block boundary
      case BoundaryFlag::periodic: // periodic boundary
      case BoundaryFlag::polar: // polar boundary
        BoundaryFunction_[BoundaryFace::inner_x2] = nullptr;
        break;
      case BoundaryFlag::polar_wedge: //polar boundary with a wedge
        BoundaryFunction_[BoundaryFace::inner_x2] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::user: // user-enrolled BCs
        BoundaryFunction_[BoundaryFace::inner_x2] =
            pmy_mesh_->BoundaryFunction_[BoundaryFace::inner_x2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix2_bc=" <<  static_cast<int>(block_bcs[BoundaryFace::inner_x2])
            << " not valid" << std::endl;
        ATHENA_ERROR(msg);
    }

    // Outer x2
    switch(block_bcs[BoundaryFace::outer_x2]) {
      case BoundaryFlag::reflect:
        BoundaryFunction_[BoundaryFace::outer_x2] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::outflow:
        BoundaryFunction_[BoundaryFace::outer_x2] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::block: // block boundary
      case BoundaryFlag::periodic: // periodic boundary
      case BoundaryFlag::polar: // polar boundary
        BoundaryFunction_[BoundaryFace::outer_x2] = nullptr;
        break;
      case BoundaryFlag::polar_wedge: // polar boundary with a wedge
        BoundaryFunction_[BoundaryFace::outer_x2] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::user: // user-enrolled BCs
        BoundaryFunction_[BoundaryFace::outer_x2] =
            pmy_mesh_->BoundaryFunction_[BoundaryFace::outer_x2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox2_bc=" <<  static_cast<int>(block_bcs[BoundaryFace::outer_x2])
            << " not valid" << std::endl;
        ATHENA_ERROR(msg);
    }
  }

  if (pmb->block_size.nx3 > 1) {
    nface_=6; nedge_=12;
    // Inner x3
    switch(block_bcs[BoundaryFace::inner_x3]) {
      case BoundaryFlag::reflect:
        BoundaryFunction_[BoundaryFace::inner_x3] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::outflow:
        BoundaryFunction_[BoundaryFace::inner_x3] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::block: // block boundary
      case BoundaryFlag::periodic: // periodic boundary
        BoundaryFunction_[BoundaryFace::inner_x3] = nullptr;
        break;
      case BoundaryFlag::user: // user-enrolled BCs
        BoundaryFunction_[BoundaryFace::inner_x3] =
            pmy_mesh_->BoundaryFunction_[BoundaryFace::inner_x3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix3_bc=" <<  static_cast<int>(block_bcs[BoundaryFace::inner_x3])
            << " not valid" << std::endl;
        ATHENA_ERROR(msg);
    }

    // Outer x3
    switch(block_bcs[BoundaryFace::outer_x3]) {
      case BoundaryFlag::reflect:
        BoundaryFunction_[BoundaryFace::outer_x3] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::outflow:
        BoundaryFunction_[BoundaryFace::outer_x3] = BValFuncPlaceholder;
        break;
      case BoundaryFlag::block: // block boundary
      case BoundaryFlag::periodic: // periodic boundary
        BoundaryFunction_[BoundaryFace::outer_x3] = nullptr;
        break;
      case BoundaryFlag::user: // user-enrolled BCs
        BoundaryFunction_[BoundaryFace::outer_x3] =
            pmy_mesh_->BoundaryFunction_[BoundaryFace::outer_x3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox3_bc=" << static_cast<int>(block_bcs[BoundaryFace::outer_x3])
            << " not valid" << std::endl;
        ATHENA_ERROR(msg);
    }
  }
  // Perform compatibilty checks of user selections of polar vs. polar_wedge boundaries
  if (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
      || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
      || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge) {
    CheckPolarBoundaries();
  }

  // Count number of blocks wrapping around pole
  if (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge) {
    int level = pmb->loc.level - pmy_mesh_->root_level;
    // KGF: possible 32-bit int overflow, if level > 31 (or possibly less, if nrbx3>1 !)
    num_north_polar_blocks_ = static_cast<int>(pmy_mesh_->nrbx3 * (1 << level));
  } else {
    num_north_polar_blocks_ = 0;
  }
  if (block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
      || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge) {
    int level = pmb->loc.level - pmy_mesh_->root_level;
    // KGF: possible 32-bit int overflow, if level > 31 (or possibly less, if nrbx3>1 !)
    num_south_polar_blocks_ = static_cast<int>(pmy_mesh_->nrbx3 * (1 << level));
  } else {
    num_south_polar_blocks_ = 0;
  }
  // end KGF: shared logic of setting boundary functions and counting spherical blocks

  // polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range
  // KGF: (fixed by Z. Zhu on 2016-01-15 in ff7b4b1)
  // KGF: shouldn't this only be allocated for MHD?
  if (pmb->loc.level == pmy_mesh_->root_level &&
      pmy_mesh_->nrbx3 == 1 &&
      (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
       || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
       || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
       || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
    azimuthal_shift_.NewAthenaArray(pmb->ke + NGHOST + 2);
  // end KGF: special handling for spherical coordinates polar boundary when nrbx3=1

  // KGF: prevent reallocation of contiguous memory space for each of 3x current calls to
  // std::vector<BoundaryVariable *>.push_back() for Hydro, Field, Gravity
  bvars.reserve(3);
  // KGF: rename to "bvars_time_int"? what about sts?
  bvars_main_int.reserve(2);

  // reserve phys=0 for former TAG_AMR=8; now hard-coded in Mesh::CreateAMRMPITag()
  bvars_next_phys_id_ = 1;

  // KGF: BVals constructor section only containing ALL shearing box-specific stuff
  // set parameters for shearing box bc and allocate buffers
//   if (SHEARING_BOX) {
//     Omega_0_ = pin->GetOrAddReal("problem","Omega0",0.001);
//     qshear_  = pin->GetOrAddReal("problem","qshear",1.5);
//     ShBoxCoord_ = pin->GetOrAddInteger("problem","shboxcoord",1);
//     x1size_ = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
//     x2size_ = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
//     x3size_ = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
//     int level = pmb->loc.level - pmy_mesh_->root_level;
//     std::int64_t nrbx1 = pmy_mesh_->nrbx1*(1L << level);
//     std::int64_t nrbx2 = pmy_mesh_->nrbx2*(1L << level);
//     shbb_.outer = false;
//     shbb_.inner = false;

//     if (ShBoxCoord_ == 1) {
//       int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
//       int ncells3 = pmb->block_size.nx3;
//       if (pmy_mesh_->mesh_size.nx3>1) ncells3 += 2*NGHOST;
//       ssize_ = NGHOST*ncells3;

//       if (pmb->loc.lx1 == 0) { // if true for shearing inner blocks
//         if (block_bcs[BoundaryFace::inner_x1] != BoundaryFlag::shear_periodic) {
//           block_bcs[BoundaryFace::inner_x1] = BoundaryFlag::shear_periodic;
//           BoundaryFunction_[BoundaryFace::inner_x1] = nullptr;
//         }
//         shboxvar_inner_hydro_.NewAthenaArray(NHYDRO,ncells3,ncells2,NGHOST);
//         flx_inner_hydro_.NewAthenaArray(ncells2);
//         if (MAGNETIC_FIELDS_ENABLED) {
//           shboxvar_inner_field_.x1f.NewAthenaArray(ncells3,ncells2,NGHOST);
//           shboxvar_inner_field_.x2f.NewAthenaArray(ncells3,ncells2+1,NGHOST);
//           shboxvar_inner_field_.x3f.NewAthenaArray(ncells3+1,ncells2,NGHOST);
//           flx_inner_field_.x1f.NewAthenaArray(ncells2);
//           flx_inner_field_.x2f.NewAthenaArray(ncells2+1);
//           flx_inner_field_.x3f.NewAthenaArray(ncells2);
//           shboxvar_inner_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
//           shboxvar_inner_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
//           shboxmap_inner_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
//           shboxmap_inner_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
//           flx_inner_emf_.x2e.NewAthenaArray(ncells2);
//           flx_inner_emf_.x3e.NewAthenaArray(ncells2+1);
//         }
//         shbb_.inner = true;
//         shbb_.igidlist=new int[nrbx2];
//         shbb_.ilidlist=new int[nrbx2];
//         shbb_.irnklist=new int[nrbx2];
//         shbb_.ilevlist=new int[nrbx2];
//         // attach corner cells from L/R side
//         int size = (pmb->block_size.nx2 + NGHOST)*ssize_*NHYDRO;
//         int bsize=0, esize=0;
//         if (MAGNETIC_FIELDS_ENABLED) {
//           // extra cell in azimuth/vertical
//           bsize = (pmb->block_size.nx2 + NGHOST+1)*(ssize_ + NGHOST)*NFIELD;
//           // face plus edge for EMF
//           esize = 2*(pmb->block_size.nx2 + NGHOST)*pmb->block_size.nx3
//                 +pmb->block_size.nx2+pmb->block_size.nx3 + NGHOST;
//         }
//         for (int n=0; n<2; n++) {
//           send_innerbuf_hydro_[n] = new Real[size];
//           recv_innerbuf_hydro_[n] = new Real[size];
//           shbox_inner_hydro_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//           rq_innersend_hydro_[n] = MPI_REQUEST_NULL;
//           rq_innerrecv_hydro_[n] = MPI_REQUEST_NULL;
// #endif
//           if (MAGNETIC_FIELDS_ENABLED) {
//             send_innerbuf_field_[n] = new Real[bsize];
//             recv_innerbuf_field_[n] = new Real[bsize];
//             shbox_inner_field_flag_[n]=BoundaryStatus::waiting;
//             send_innerbuf_emf_[n] = new Real[esize];
//             recv_innerbuf_emf_[n] = new Real[esize];
//             shbox_inner_emf_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//             rq_innersend_field_[n] = MPI_REQUEST_NULL;
//             rq_innerrecv_field_[n] = MPI_REQUEST_NULL;
//             rq_innersend_emf_[n] = MPI_REQUEST_NULL;
//             rq_innerrecv_emf_[n] = MPI_REQUEST_NULL;
// #endif
//           }
//         }
//         size = NGHOST*ssize_*NHYDRO;// corner cells only
//         if (MAGNETIC_FIELDS_ENABLED) {
//             bsize = NGHOST*(ssize_ + NGHOST)*NFIELD;
//             esize = 2*NGHOST*pmb->block_size.nx3 + NGHOST;
//         }
//         for (int n=2; n<4; n++) {
//           send_innerbuf_hydro_[n] = new Real[size];
//           recv_innerbuf_hydro_[n] = new Real[size];
//           shbox_inner_hydro_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//           rq_innersend_hydro_[n] = MPI_REQUEST_NULL;
//           rq_innerrecv_hydro_[n] = MPI_REQUEST_NULL;
// #endif
//           if (MAGNETIC_FIELDS_ENABLED) {
//             send_innerbuf_field_[n] = new Real[bsize];
//             recv_innerbuf_field_[n] = new Real[bsize];
//             shbox_inner_field_flag_[n]=BoundaryStatus::waiting;
//             send_innerbuf_emf_[n] = new Real[esize];
//             recv_innerbuf_emf_[n] = new Real[esize];
//             shbox_inner_emf_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//             rq_innersend_field_[n] = MPI_REQUEST_NULL;
//             rq_innerrecv_field_[n] = MPI_REQUEST_NULL;
//             rq_innersend_emf_[n] = MPI_REQUEST_NULL;
//             rq_innerrecv_emf_[n] = MPI_REQUEST_NULL;
// #endif
//           }
//         }
//       }

//       if (pmb->loc.lx1 == (nrbx1-1)) { // if true for shearing outer blocks
//         if (block_bcs[BoundaryFace::outer_x1] != BoundaryFlag::shear_periodic) {
//           block_bcs[BoundaryFace::outer_x1] = BoundaryFlag::shear_periodic;
//           BoundaryFunction_[BoundaryFace::outer_x1] = nullptr;
//         }
//         shboxvar_outer_hydro_.NewAthenaArray(NHYDRO,ncells3,ncells2,NGHOST);
//         flx_outer_hydro_.NewAthenaArray(ncells2);
//         if (MAGNETIC_FIELDS_ENABLED) {
//           shboxvar_outer_field_.x1f.NewAthenaArray(ncells3,ncells2,NGHOST);
//           shboxvar_outer_field_.x2f.NewAthenaArray(ncells3,ncells2+1,NGHOST);
//           shboxvar_outer_field_.x3f.NewAthenaArray(ncells3+1,ncells2,NGHOST);
//           flx_outer_field_.x1f.NewAthenaArray(ncells2);
//           flx_outer_field_.x2f.NewAthenaArray(ncells2+1);
//           flx_outer_field_.x3f.NewAthenaArray(ncells2);
//           shboxvar_outer_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
//           shboxvar_outer_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
//           shboxmap_outer_emf_.x2e.NewAthenaArray(ncells3+1,ncells2);
//           shboxmap_outer_emf_.x3e.NewAthenaArray(ncells3,ncells2+1);
//           flx_outer_emf_.x2e.NewAthenaArray(ncells2);
//           flx_outer_emf_.x3e.NewAthenaArray(ncells2+1);
//         }
//         shbb_.outer = true;
//         shbb_.ogidlist=new int[nrbx2];
//         shbb_.olidlist=new int[nrbx2];
//         shbb_.ornklist=new int[nrbx2];
//         shbb_.olevlist=new int[nrbx2];
//         // attach corner cells from L/R side
//         int size = (pmb->block_size.nx2 + NGHOST)*ssize_*NHYDRO;
//         int bsize=0, esize=0;
//         if (MAGNETIC_FIELDS_ENABLED) {
//           // extra cell in azimuth/vertical
//           bsize = (pmb->block_size.nx2 + NGHOST+1)*(ssize_ + NGHOST)*NFIELD;
//           // face plus edge for EMF
//           esize = 2*(pmb->block_size.nx2 + NGHOST)*pmb->block_size.nx3
//                 +pmb->block_size.nx2+pmb->block_size.nx3 + NGHOST;
//         }
//         for (int n=0; n<2; n++) {
//           send_outerbuf_hydro_[n] = new Real[size];
//           recv_outerbuf_hydro_[n] = new Real[size];
//           shbox_outer_hydro_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//           rq_outersend_hydro_[n] = MPI_REQUEST_NULL;
//           rq_outerrecv_hydro_[n] = MPI_REQUEST_NULL;
// #endif
//           if (MAGNETIC_FIELDS_ENABLED) {
//             send_outerbuf_field_[n] = new Real[bsize];
//             recv_outerbuf_field_[n] = new Real[bsize];
//             shbox_outer_field_flag_[n]=BoundaryStatus::waiting;
//             send_outerbuf_emf_[n] = new Real[esize];
//             recv_outerbuf_emf_[n] = new Real[esize];
//             shbox_outer_emf_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//             rq_outersend_field_[n] = MPI_REQUEST_NULL;
//             rq_outerrecv_field_[n] = MPI_REQUEST_NULL;
//             rq_outersend_emf_[n] = MPI_REQUEST_NULL;
//             rq_outerrecv_emf_[n] = MPI_REQUEST_NULL;
// #endif
//           }
//         }
//         size = NGHOST*ssize_*NHYDRO;// corner cells only
//         if (MAGNETIC_FIELDS_ENABLED) {
//           bsize = NGHOST*(ssize_ + NGHOST)*NFIELD;
//           esize = 2*NGHOST*pmb->block_size.nx3 + NGHOST;
//         }
//         for (int n=2; n<4; n++) {
//           send_outerbuf_hydro_[n] = new Real[size];
//           recv_outerbuf_hydro_[n] = new Real[size];
//           shbox_outer_hydro_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//           rq_outersend_hydro_[n] = MPI_REQUEST_NULL;
//           rq_outerrecv_hydro_[n] = MPI_REQUEST_NULL;
// #endif
//           if (MAGNETIC_FIELDS_ENABLED) {
//             send_outerbuf_field_[n] = new Real[bsize];
//             recv_outerbuf_field_[n] = new Real[bsize];
//             shbox_outer_field_flag_[n]=BoundaryStatus::waiting;
//             send_outerbuf_emf_[n] = new Real[esize];
//             recv_outerbuf_emf_[n] = new Real[esize];
//             shbox_outer_emf_flag_[n]=BoundaryStatus::waiting;
// #ifdef MPI_PARALLEL
//             rq_outersend_field_[n] = MPI_REQUEST_NULL;
//             rq_outerrecv_field_[n] = MPI_REQUEST_NULL;
//             rq_outersend_emf_[n] = MPI_REQUEST_NULL;
//             rq_outerrecv_emf_[n] = MPI_REQUEST_NULL;
// #endif
//           }
//         }
//       }
//     }
//   } // end KGF: shearing box in BoundaryValues constructor
}

// destructor

BoundaryValues::~BoundaryValues() {
  MeshBlock *pmb = pmy_block_;

  // KGF: edge-case of single block across pole in MHD spherical polar coordinates
  // Note, this conditional is outside "if MAGNETIC_FIELDS_ENABLED" in master (also
  // true for its counterpart in constructor). Probably should be inside.
  if (pmb->loc.level == pmy_mesh_->root_level &&
      pmy_mesh_->nrbx3 == 1 &&
      (block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
       || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
       || block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
       || block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
    azimuthal_shift_.DeleteAthenaArray();
  // end KGF: edge-case...

  // end KGF: destructor counterpart of special handling of emf in spherical polar

  // KGF: shearing box destructor
  // if (SHEARING_BOX) {
  //   int level = pmb->loc.level - pmb->pmy_mesh->root_level;
  //   std::int64_t nrbx1 = pmb->pmy_mesh->nrbx1*(1L << level);
  //   if (pmb->loc.lx1 == 0) { // if true for shearing inner blocks
  //     shboxvar_inner_hydro_.DeleteAthenaArray();
  //     flx_inner_hydro_.DeleteAthenaArray();
  //     for (int n=0; n<4; n++) {
  //       delete[] send_innerbuf_hydro_[n];
  //       delete[] recv_innerbuf_hydro_[n];
  //     }
  //     if (MAGNETIC_FIELDS_ENABLED) {
  //       shboxvar_inner_field_.x1f.DeleteAthenaArray();
  //       shboxvar_inner_field_.x2f.DeleteAthenaArray();
  //       shboxvar_inner_field_.x3f.DeleteAthenaArray();
  //       flx_inner_field_.x1f.DeleteAthenaArray();
  //       flx_inner_field_.x2f.DeleteAthenaArray();
  //       flx_inner_field_.x3f.DeleteAthenaArray();
  //       shboxvar_inner_emf_.x2e.DeleteAthenaArray();
  //       shboxvar_inner_emf_.x3e.DeleteAthenaArray();
  //       flx_inner_emf_.x2e.DeleteAthenaArray();
  //       flx_inner_emf_.x3e.DeleteAthenaArray();
  //       for (int n=0; n<4; n++) {
  //         delete[] send_innerbuf_field_[n];
  //         delete[] recv_innerbuf_field_[n];
  //         delete[] send_innerbuf_emf_[n];
  //         delete[] recv_innerbuf_emf_[n];
  //       }
  //     }
  //   }
  //   if (pmb->loc.lx1 == (nrbx1-1)) { // if true for shearing outer blocks
  //     shboxvar_outer_hydro_.DeleteAthenaArray();
  //     flx_outer_hydro_.DeleteAthenaArray();
  //     for (int n=0; n<4; n++) {
  //       delete[] send_outerbuf_hydro_[n];
  //       delete[] recv_outerbuf_hydro_[n];
  //     }
  //     if (MAGNETIC_FIELDS_ENABLED) {
  //       shboxvar_outer_field_.x1f.DeleteAthenaArray();
  //       shboxvar_outer_field_.x2f.DeleteAthenaArray();
  //       shboxvar_outer_field_.x3f.DeleteAthenaArray();
  //       flx_outer_field_.x1f.DeleteAthenaArray();
  //       flx_outer_field_.x2f.DeleteAthenaArray();
  //       flx_outer_field_.x3f.DeleteAthenaArray();
  //       shboxvar_outer_emf_.x2e.DeleteAthenaArray();
  //       shboxvar_outer_emf_.x3e.DeleteAthenaArray();
  //       flx_outer_emf_.x2e.DeleteAthenaArray();
  //       flx_outer_emf_.x3e.DeleteAthenaArray();
  //       for (int n=0; n<4; n++) {
  //         delete[] send_outerbuf_field_[n];
  //         delete[] recv_outerbuf_field_[n];
  //         delete[] send_outerbuf_emf_[n];
  //         delete[] recv_outerbuf_emf_[n];
  //       }
  //     }
  //   }
  // } // KGF: end shearing box handling in destructor
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetupPersistentMPI()
//  \brief Setup persistent MPI requests
// TODO(felker): rename to a less generic name to avoid confusion with InitBoundaryData
// KGF: called in Mesh::Initialize(), after CheckBoundary()
//      and before StartReceivingForInit(true)
void BoundaryValues::SetupPersistentMPI() {
  // KGF: (although the counting is for 2x BoundaryVariables private members that are only
  // used for emf purposes in flux_correction_fc.cpp)
  // num_north_polar_blocks_, num_south_polar_blocks_, nedge_, nface_ are calculated in
  // BoundaryValues() constructor.

  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    (*bvars_it)->SetupPersistentMPI();
  }

  // KGF: begin exclusive shearing-box section in BoundaryValues::SetupPersistentMPI()

  // initialize the shearing block lists
  // if (SHEARING_BOX) {
  //   Mesh *pmesh = pmb->pmy_mesh;
  //   int level = pmb->loc.level - pmesh->root_level;
  //   std::int64_t nrbx1 = pmesh->nrbx1*(1L << level);
  //   // std::int64_t nrbx2 = pmesh->nrbx2*(1L << level); // unused variable
  //   int nbtotal = pmesh->nbtotal;
  //   int *ranklist = pmesh->ranklist;
  //   int *nslist = pmesh->nslist;
  //   LogicalLocation *loclist = pmesh->loclist;

  //   int count = 0;
  //   if (shbb_.inner) {
  //     for (int i=0;i<nbtotal;i++) {
  //       if (loclist[i].lx1 == 0 && loclist[i].lx3 == pmb->loc.lx3 &&
  //           loclist[i].level == pmb->loc.level) {
  //         shbb_.igidlist[count] = i;
  //         shbb_.ilidlist[count] = i - nslist[ranklist[i]];
  //         shbb_.irnklist[count] = ranklist[i];
  //         shbb_.ilevlist[count] = loclist[i].level;
  //         count++;
  //       }
  //     }
  //   }
  //   count = 0;
  //   if (shbb_.outer) {
  //     for (int i=0;i<nbtotal;i++) {
  //       if (loclist[i].lx1 == (nrbx1-1) && loclist[i].lx3 == pmb->loc.lx3 &&
  //         loclist[i].level == pmb->loc.level) {
  //         shbb_.ogidlist[count] = i;
  //         shbb_.olidlist[count] = i - nslist[ranklist[i]];
  //         shbb_.ornklist[count] = ranklist[i];
  //         shbb_.olevlist[count] = loclist[i].level;
  //         count++;
  //       }
  //     }
  //   }
  // } // end KGF: exclusive shearing box portion of SetupPersistentMPI()
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckBoundary()
//  \brief checks if the boundary conditions are correctly enrolled (and other boundary
//  values compatibility checks performed at the top of Mesh::Initialize())

void BoundaryValues::CheckBoundary() {
  for (int i=0; i<nface_; i++) {
    if (block_bcs[i]==BoundaryFlag::user) {
      if (BoundaryFunction_[i]==nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the hydro boundary function "
            << "is not enrolled in direction " << i  << "." << std::endl;
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
  // KGF: approach #1: call each fn of element of bvar vector inside #ifdef and loop
// #ifdef MPI_PARALLEL
//   MeshBlock *pmb = pmy_block_;
//   int mylevel = pmb->loc.level;
//   for (int n=0; n<nneighbor; n++) {
//     NeighborBlock& nb = neighbor[n];
//     if (nb.rank!=Globals::my_rank) {
//       for (auto bvars_it = bvars.begin(); bvars_it != bvars.end(); ++bvars_it) {
//         bvars_it->StartReceiving(BoundaryCommSubset phase);
//       }
//     }
//   }
// #endif

  // KGF: approach #2: make loop over bvar vector the outermost loop; separate,
  // independent loops over nneighbor
  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    (*bvars_it)->StartReceiving(phase);
  }

  // KGF: begin shearing-box exclusive section of original StartReceivingForInit()
  // find send_block_id and recv_block_id;
  // if (SHEARING_BOX) {
  //   MeshBlock *pmb = pmy_block_;
  //   Mesh *pmesh = pmb->pmy_mesh;
  //   FindShearBlock(pmesh->time);
  // }
  // end KGF: shearing box

  // KGF: begin shearing-box exclusive section of original StartReceivingAll()
  // find send_block_id and recv_block_id; post non-blocking recv
//   if (SHEARING_BOX) {
//     FindShearBlock(time);
// #ifdef MPI_PARALLEL
//      //     Mesh *pmesh = pmb->pmy_mesh;
//     int size,tag;
//     if (shbb_.inner) { // inner boundary
//       for (int n=0; n<4; n++) {
//         if ((recv_inner_rank_[n]!=Globals::my_rank) &&
//                           (recv_inner_rank_[n]!=-1)) {
//           size = ssize_*NHYDRO*recv_innersize_hydro_[n];
//           tag  = CreateBvalsMPITag(pmb->lid, n, AthenaTagMPI::shbox_hydro);
//           MPI_Irecv(recv_innerbuf_hydro_[n],size,MPI_ATHENA_REAL,
//                     recv_inner_rank_[n],tag,MPI_COMM_WORLD,
//                     &rq_innerrecv_hydro_[n]);
//           if (MAGNETIC_FIELDS_ENABLED) {
//             size = recv_innersize_field_[n];
//             tag  = CreateBvalsMPITag(pmb->lid, n, AthenaTagMPI::shbox_field);
//             MPI_Irecv(recv_innerbuf_field_[n],size,MPI_ATHENA_REAL,
//                       recv_inner_rank_[n],tag,MPI_COMM_WORLD,
//                       &rq_innerrecv_field_[n]);
//             size = recv_innersize_emf_[n];
//             tag  = CreateBvalsMPITag(pmb->lid, n, AthenaTagMPI::shbox_emf);
//             MPI_Irecv(recv_innerbuf_emf_[n],size,MPI_ATHENA_REAL,
//                       recv_inner_rank_[n],tag,MPI_COMM_WORLD,
//                       &rq_innerrecv_emf_[n]);
//           }
//         }
//       }
//     }

//     if (shbb_.outer) { // outer boundary
//       int offset=4;
//       for (int n=0; n<4; n++) {
//         if ((recv_outer_rank_[n]!=Globals::my_rank) &&
//                           (recv_outer_rank_[n]!=-1)) {
//           size = ssize_*NHYDRO*recv_outersize_hydro_[n];
//           tag  = CreateBvalsMPITag(pmb->lid, n+offset, AthenaTagMPI::shbox_hydro);
//           MPI_Irecv(recv_outerbuf_hydro_[n],size,MPI_ATHENA_REAL,
//                     recv_outer_rank_[n],tag,MPI_COMM_WORLD,
//                     &rq_outerrecv_hydro_[n]);
//           if (MAGNETIC_FIELDS_ENABLED) {
//             size = recv_outersize_field_[n];
//             tag  = CreateBvalsMPITag(pmb->lid, n+offset, AthenaTagMPI::shbox_field);
//             MPI_Irecv(recv_outerbuf_field_[n],size,MPI_ATHENA_REAL,
//                       recv_outer_rank_[n],tag,MPI_COMM_WORLD,
//                       &rq_outerrecv_field_[n]);
//             size = recv_outersize_emf_[n];
//             tag  = CreateBvalsMPITag(pmb->lid, n+offset, AthenaTagMPI::shbox_emf);
//             MPI_Irecv(recv_outerbuf_emf_[n],size,MPI_ATHENA_REAL,
//                       recv_outer_rank_[n],tag,MPI_COMM_WORLD,
//                       &rq_outerrecv_emf_[n]);
//           }
//         }
//       }
//     }
// #endif
//   } // end KGF: shearing-box exclusive section of StartReceivingAll
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

  // KGF: begin shearing-box exclusive section of ClearBoundaryAll
  // clear shearingbox boundary communications
//   if (SHEARING_BOX) {
//     if (shbb_.inner == true) {
//       for (int n=0; n<4; n++) {
//         if (send_inner_rank_[n] == -1) continue;
//         shbox_inner_hydro_flag_[n] = BoundaryStatus::waiting;
//         if (MAGNETIC_FIELDS_ENABLED) {
//           shbox_inner_field_flag_[n] = BoundaryStatus::waiting;
//           shbox_inner_emf_flag_[n] = BoundaryStatus::waiting;
//         }
// #ifdef MPI_PARALLEL
//         if (send_inner_rank_[n]!=Globals::my_rank) {
//           MPI_Wait(&rq_innersend_hydro_[n],MPI_STATUS_IGNORE);
//           if (MAGNETIC_FIELDS_ENABLED) {
//             MPI_Wait(&rq_innersend_field_[n],MPI_STATUS_IGNORE);
//             MPI_Wait(&rq_innersend_emf_[n],MPI_STATUS_IGNORE);
//           }
//         }
// #endif
//       }
//     } // inner boundary

//     if (shbb_.outer == true) {
//       for (int n=0; n<4; n++) {
//         if (send_outer_rank_[n] == -1) continue;
//         shbox_outer_hydro_flag_[n] = BoundaryStatus::waiting;
//         if (MAGNETIC_FIELDS_ENABLED) {
//           shbox_outer_field_flag_[n] = BoundaryStatus::waiting;
//         }
// #ifdef MPI_PARALLEL
//         if (send_outer_rank_[n]!=Globals::my_rank) {
//           MPI_Wait(&rq_outersend_hydro_[n],MPI_STATUS_IGNORE);
//           if (MAGNETIC_FIELDS_ENABLED) {
//             MPI_Wait(&rq_outersend_field_[n],MPI_STATUS_IGNORE);
//             MPI_Wait(&rq_outersend_emf_[n],MPI_STATUS_IGNORE);
//           }
//         }
// #endif
//       }
//     }
//   } // end KGF: shearing box
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ApplyPhysicalBoundaries(AthenaArray<Real> &pdst,
//           AthenaArray<Real> &cdst, FaceField &bfdst, AthenaArray<Real> &bcdst,
//           const Real time, const Real dt)
//  \brief Apply all the physical boundary conditions for both hydro and field

void BoundaryValues::ApplyPhysicalBoundaries(const Real time, const Real dt) {
  // AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
  // FaceField &bfdst, AthenaArray<Real> &bcdst,
  // const Real time, const Real dt) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int bis = pmb->is - NGHOST, bie = pmb->ie + NGHOST,
      bjs = pmb->js, bje = pmb->je,
      bks = pmb->ks, bke = pmb->ke;

  // Extend the transverse limits that correspond to periodic boundaries as they are
  // updated

  // x1, then x2, then x3
  if (BoundaryFunction_[BoundaryFace::inner_x2] == nullptr && pmb->block_size.nx2 > 1)
    bjs = pmb->js - NGHOST;
  if (BoundaryFunction_[BoundaryFace::outer_x2] == nullptr && pmb->block_size.nx2 > 1)
    bje = pmb->je + NGHOST;
  if (BoundaryFunction_[BoundaryFace::inner_x3] == nullptr && pmb->block_size.nx3 > 1)
    bks = pmb->ks - NGHOST;
  if (BoundaryFunction_[BoundaryFace::outer_x3] == nullptr && pmb->block_size.nx3 > 1)
    bke = pmb->ke + NGHOST;

  // KGF: temporarily hardcode Hydro and Field access for coupling in EOS, and when passed
  // to user-defined boundary function stored in function pointer array

  // downcast BoundaryVariable pointers to known derived class pointer types:
  // RTTI via dynamic_cast
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
  if (BoundaryFunction_[BoundaryFace::inner_x1] != nullptr) {
    switch(block_bcs[BoundaryFace::inner_x1]) {
      case BoundaryFlag::reflect:
        for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
             ++bvars_it) {
          (*bvars_it)->ReflectInnerX1(pmb, pco, time, dt, pmb->is, pmb->ie,
                                      bjs, bje, bks, bke, NGHOST);
        }
        break;
      case BoundaryFlag::outflow:
        for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
             ++bvars_it) {
          (*bvars_it)->OutflowInnerX1(pmb, pco, time, dt, pmb->is, pmb->ie,
                                      bjs, bje, bks, bke, NGHOST);
        }
        break;
      case BoundaryFlag::user: // user-enrolled BCs
        BoundaryFunction_[BoundaryFace::inner_x1] (
            pmb, pco, ph->w, pf->b, time, dt,
            pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST);
        break;
      default:
        break;
    }
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
  if (BoundaryFunction_[BoundaryFace::outer_x1] != nullptr) {
    switch(block_bcs[BoundaryFace::outer_x1]) {
      case BoundaryFlag::reflect:
        for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
             ++bvars_it) {
          (*bvars_it)->ReflectOuterX1(pmb, pco, time, dt, pmb->is, pmb->ie,
                                      bjs, bje, bks, bke, NGHOST);
        }
        break;
      case BoundaryFlag::outflow:
        for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
             ++bvars_it) {
          (*bvars_it)->OutflowOuterX1(pmb, pco, time, dt, pmb->is, pmb->ie,
                                      bjs, bje, bks, bke, NGHOST);
        }
        break;
      case BoundaryFlag::user: // user-enrolled BCs
        BoundaryFunction_[BoundaryFace::outer_x1](
            pmb, pco, ph->w, pf->b, time, dt,
            pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST);
        break;
      default:
        break;
    }
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
    if (BoundaryFunction_[BoundaryFace::inner_x2] != nullptr) {
      switch(block_bcs[BoundaryFace::inner_x2]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectInnerX2(pmb, pco, time, dt, bis, bie,
                                        pmb->js, pmb->je, bks, bke, NGHOST);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowInnerX2(pmb, pco, time, dt, bis, bie,
                                        pmb->js, pmb->je, bks, bke, NGHOST);
          }
          break;
        case BoundaryFlag::polar_wedge:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->PolarWedgeInnerX2(pmb, pco, time, dt, bis, bie,
                                           pmb->js, pmb->je, bks, bke, NGHOST);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          BoundaryFunction_[BoundaryFace::inner_x2](pmb, pco, ph->w, pf->b, time, dt,
                                                    bis, bie, pmb->js, pmb->je, bks, bke,
                                                    NGHOST);
          break;
        default:
          break;
      }
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
    if (BoundaryFunction_[BoundaryFace::outer_x2] != nullptr) {
      switch(block_bcs[BoundaryFace::outer_x2]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectOuterX2(pmb, pco, time, dt, bis, bie,
                                        pmb->js, pmb->je, bks, bke, NGHOST);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowOuterX2(pmb, pco, time, dt, bis, bie,
                                        pmb->js, pmb->je, bks, bke, NGHOST);
          }
          break;
        case BoundaryFlag::polar_wedge:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->PolarWedgeOuterX2(pmb, pco, time, dt, bis, bie,
                                           pmb->js, pmb->je, bks, bke, NGHOST);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          BoundaryFunction_[BoundaryFace::outer_x2](pmb, pco, ph->w, pf->b, time, dt,
                                                    bis, bie, pmb->js, pmb->je, bks, bke,
                                                    NGHOST);
          break;
        default:
          break;
      }
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
    if (BoundaryFunction_[BoundaryFace::inner_x3] != nullptr) {
      switch(block_bcs[BoundaryFace::inner_x3]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectInnerX3(pmb, pco, time, dt, bis, bie,
                                        bjs, bje, pmb->ks, pmb->ke, NGHOST);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowInnerX3(pmb, pco, time, dt, bis, bie,
                                        bjs, bje, pmb->ks, pmb->ke, NGHOST);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          BoundaryFunction_[BoundaryFace::inner_x3](
              pmb, pco, ph->w, pf->b, time, dt,
              bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST);
          break;
        default:
          break;
      }
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
    if (BoundaryFunction_[BoundaryFace::outer_x3] != nullptr) {
      switch(block_bcs[BoundaryFace::outer_x3]) {
        case BoundaryFlag::reflect:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->ReflectOuterX3(pmb, pco, time, dt, bis, bie,
                                        bjs, bje, pmb->ks, pmb->ke, NGHOST);
          }
          break;
        case BoundaryFlag::outflow:
          for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
               ++bvars_it) {
            (*bvars_it)->OutflowOuterX3(pmb, pco, time, dt, bis, bie,
                                        bjs, bje, pmb->ks, pmb->ke, NGHOST);
          }
          break;
        case BoundaryFlag::user: // user-enrolled BCs
          BoundaryFunction_[BoundaryFace::outer_x3](
              pmb, pco, ph->w, pf->b, time, dt,
              bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST);
          break;
        default:
          break;
      }
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
  // KGF: add safety checks? input, output are positive, obey <= 31= MAX_NUM_PHYS
  int start_id = bvars_next_phys_id_;
  bvars_next_phys_id_ += num_phys;
  return start_id;
}
