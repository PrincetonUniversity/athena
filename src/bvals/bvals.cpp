//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.cpp
//! \brief constructor/destructor and utility functions for BoundaryValues class

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
#include <utility>    // swap()
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cr/cr.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../multigrid/multigrid.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/buffer_utils.hpp"
#include "bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \brief BoundaryValues constructor
//!        (the first object constructed inside the MeshBlock() constructor)
//!
//! Sets functions for the appropriate boundary conditions at each of the 6
//! dirs of a MeshBlock
//!
//! At the end, there is a section containing ALL shearing box-specific stuff:
//! set parameters for shearing box bc and allocate buffers
BoundaryValues::BoundaryValues(MeshBlock *pmb, BoundaryFlag *input_bcs,
                               ParameterInput *pin)
    : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs), pmy_block_(pmb),
      is_shear{}, loc_shear{0, pmy_mesh_->nrbx1*(1L << (pmb->loc.level -
                                                        pmy_mesh_->root_level)) - 1},
      sb_data_{}, sb_flux_data_{} {
  // Check BC functions for each of the 6 boundaries in turn ---------------------
  for (int i=0; i<6; i++) {
    switch (block_bcs[i]) {
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
    nface_ = 6; nedge_ = 12;
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
  bvars_main_int.reserve(2);
  if (STS_ENABLED) {
    bvars_sts.reserve(1);
  }

  // Matches initial value of Mesh::next_phys_id_
  // reserve phys=0 for former TAG_AMR=8; now hard-coded in Mesh::CreateAMRMPITag()
  bvars_next_phys_id_ = 1;

  // BVals constructor section only containing ALL shearing box-specific stuff
  // set parameters for shearing box bc and allocate buffers
  shearing_box = 0;
  if (pmy_mesh_->shear_periodic) {
    // It is required to know the reconstruction scheme before pmb->precon is defined.
    // If a higher-order scheme than PPM is implemented, xgh_ must be larger than 2.
    std::string input_recon = pin->GetOrAddString("time", "xorder", "2");
    if (input_recon == "1") {
      xorder_ = 1;
    } else if ((input_recon == "2") || (input_recon == "2c")) {
      xorder_ = 2;
    } else if ((input_recon == "3") || (input_recon == "3c")) {
      xorder_ = 3;
    } else if ((input_recon == "4") || (input_recon == "4c")) {
      xorder_ = 4;
    }
    if (xorder_ <= 2) xgh_ = 1;
    else
      xgh_ = 2;

    shearing_box = pin->GetOrAddInteger("orbital_advection","shboxcoord",1);
    if (shearing_box != 1 && shearing_box != 2) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
          << "<orbital_advection> shboxcoord must be 1 or 2." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (pmb->block_size.nx3>1) { // 3D
      if (shearing_box == 2) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
            << "When using shear_periodic bondaries in 3D, "
            << "<orbital_advection> shboxcoord must be 1." << std::endl;
        ATHENA_ERROR(msg);
      }
    } else if (pmb->block_size.nx2==1) { // 1D
      if (shearing_box == 1) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
            << "When using shear_periodic bondaries in 1D, "
            << "<orbital_advection> shboxcoord must be 2." << std::endl;
        ATHENA_ERROR(msg);
      }
    }

    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
          << "shear_periodic boundaries work only in cartesian coordinates."<<std::endl;
      ATHENA_ERROR(msg);
    }

    if (!pmy_mesh_->use_uniform_meshgen_fn_[1]) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
          << "shear_periodic boundaries work only with uniform spacing "
          << "in the x2 direction." << std::endl
          << "Check <mesh> x2rat parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }

    std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
    if (integrator == "rk4" || integrator == "ssprk5_4") {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
          << "shear_periodic boundaries do not work with rk4 or ssprk5_4." << std::endl
          << "Check <time> integrator parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (pmy_mesh_->OrbitalVelocity_ != nullptr) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
          << "shear_periodic boundaries require the default orbital velocity profile."
          << std::endl;
      ATHENA_ERROR(msg);
    }

    if (RELATIVISTIC_DYNAMICS) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class."<<std::endl
          << "Neigher SR nor GR works with Shear Periodic."<<std::endl;
      ATHENA_ERROR(msg);
    }

    int level = pmb->loc.level - pmy_mesh_->root_level;
    // nblx2 is the only var used in the ctor for alllocating SimpleNeighborBlock arrays
    nblx1 = pmy_mesh_->nrbx1*(1L << level);
    nblx2 = pmy_mesh_->nrbx2*(1L << level);
    nblx3 = pmy_mesh_->nrbx3*(1L << level);

    if (shearing_box == 1) {
      if (NGHOST+xgh_ > pmb->block_size.nx2) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues Class."<<std::endl
            << "x2 block_size must be larger than NGHOST + "<< xgh_
            << " = "<< NGHOST+xgh_ << "with shear_periodic boundary."<<std::endl;
        ATHENA_ERROR(msg);
      }
      int pnum = pmb->block_size.nx2+2*NGHOST+1;
      if (MAGNETIC_FIELDS_ENABLED) pnum++;
      pflux_.NewAthenaArray(pnum);

      int nc3 = pmb->ncells3;
      ssize_ = NGHOST*nc3;
      //! \todo (felker):
      //! * much of this should be a part of InitBoundaryData()
      for (int upper=0; upper<2; upper++) {
        if (pmb->loc.lx1 == loc_shear[upper]) { // if true for shearing inner blocks
          is_shear[upper] = true;
          shbb_[upper] = new SimpleNeighborBlock[nblx2];
        } // end "if is a shearing boundary"
      } // end loop over inner, outer shearing boundaries
    } // end "if (shearing_box == 1)"
  } // end shearing box component of BoundaryValues ctor
}

//----------------------------------------------------------------------------------------
//! destructor

BoundaryValues::~BoundaryValues() {
  if (shearing_box == 1) {
    for (int upper=0; upper<2; upper++)
      if (is_shear[upper]) delete[] shbb_[upper];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetupPersistentMPI()
//! \brief Setup persistent MPI requests to be reused throughout the entire simulation
//!
//! Includes an exclusive shearing-box section that
//! initializes the shearing block lists

void BoundaryValues::SetupPersistentMPI() {
  // TODO(KGF): given the fn name, it is confusing that the fn contents are not entirely
  // wrapped in #ifdef MPI_PARALLEL
  for (auto bvars_it = bvars_main_int.begin(); bvars_it != bvars_main_int.end();
       ++bvars_it) {
    (*bvars_it)->SetupPersistentMPI();
  }
  // KGF: begin exclusive shearing-box section in BoundaryValues::SetupPersistentMPI()

  // TODO(KGF): the shearing box metadata is required regardless if MPI is used, so better
  // to encapsulate it in a separate function

  // initialize the shearing block lists
  if (shearing_box != 0) {
    MeshBlock *pmb = pmy_block_;
    int *ranklist = pmy_mesh_->ranklist;
    int *nslist = pmy_mesh_->nslist;
    LogicalLocation *loclist = pmy_mesh_->loclist;
    // this error needs to be downgraded to a non-fatal warning if the AMR+SMR workaround
    // to the MHD+shear+refinement restriction is employed (#569)--- or remove entirely,
    // hoping that RefinementCondition() does not violate below restrictions
    if (pmy_mesh_->adaptive) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
          << "shear_periodic boundaries do NOT work with AMR at present." << std::endl;
      ATHENA_ERROR(msg);
    } else {
      for (int upper=0; upper<2; upper++) {
        if (is_shear[upper]) {
          LogicalLocation loc;
          loc.level = pmb->loc.level;
          loc.lx1   = loc_shear[upper];
          loc.lx2   = pmb->loc.lx2;
          for (int64_t lx3=0; lx3<nblx3; lx3++) {
            loc.lx3   = lx3;
            MeshBlockTree *mbt = pmy_mesh_->tree.FindMeshBlock(loc);
            // see discussion #569: hydro+shear with SMR (not AMR) works with x1, x3
            // refinement but MHD with x3 refinement will NOT work without adjusting
            // corner E-field flux correction to also take into account adjacent x3 block
            // emf (currently just averages emf values between inner/outermost x1 MBs)
            if ((mbt == nullptr || mbt->GetGid() == -1) && MAGNETIC_FIELDS_ENABLED) {
              std::stringstream msg;
              msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
                  << "shear_periodic boundaries do NOT work with MHD"
                  << "if there is refinment in the x3 direction." << std::endl;
              ATHENA_ERROR(msg);
            }
          }
          loc.lx3   = pmb->loc.lx3;
          for (int64_t lx2=0; lx2<nblx2; lx2++) {
            loc.lx2 = lx2;
            MeshBlockTree *mbt = pmy_mesh_->tree.FindMeshBlock(loc);
            if (mbt == nullptr || mbt->GetGid() == -1) {
              std::stringstream msg;
              msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
                  << "shear_periodic boundaries do NOT work "
                  << "if there is refinment in the x2 direction." << std::endl;
              ATHENA_ERROR(msg);
            }
            int gid = mbt->GetGid();
            shbb_[upper][lx2].gid = gid;
            shbb_[upper][lx2].lid = gid - nslist[ranklist[gid]];
            shbb_[upper][lx2].rank = ranklist[gid];
            shbb_[upper][lx2].level = loclist[gid].level;
          }
          loc.lx1 = loc_shear[1-upper];
          loc.lx2   = pmb->loc.lx2;
          MeshBlockTree *mbt = pmy_mesh_->tree.FindMeshBlock(loc);
          if (mbt == nullptr || mbt->GetGid() == -1) {
            std::stringstream msg;
            msg << "### FATAL ERROR in BoundaryValues Class" << std::endl
                << "shear_periodic boundaries do NOT work if MeshBlocks contacting "
                << "the shear boundaries are on different levels." << std::endl;
            ATHENA_ERROR(msg);
          }
        }
      }
    }
    qomL_ = pmb->porb->OrbitalVelocity(pmb->porb,pmy_mesh_->mesh_size.x1min,0,0)
              - pmb->porb->OrbitalVelocity(pmb->porb,pmy_mesh_->mesh_size.x1max,0,0);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckUserBoundaries()
//! \brief checks if the boundary functions are correctly enrolled
//!
//! This compatibility check is performed at the top of Mesh::Initialize(),
//! after calling ProblemGenerator()

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

      if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
        if (pmy_mesh_->RadBoundaryFunc_[i] == nullptr) {
          std::stringstream msg;
          msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
              << "A user-defined boundary is specified but the actual RadBoundaryFunc_ "
              << "is not enrolled in direction " << i  << " (in [0,6])."
              << std::endl;
          ATHENA_ERROR(msg);
        }
      }
      if (CR_ENABLED) {
        if (pmy_mesh_->CRBoundaryFunc_[i] == nullptr) {
          std::stringstream msg;
          msg << "### FATAL ERROR in BoundaryValues::CheckBoundary" << std::endl
              << "A user-defined boundary is specified but the actual CRBoundaryFunc_ "
              << "is not enrolled in direction " << i  << " (in [0,6])."
              << std::endl;
          ATHENA_ERROR(msg);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingSubset(BoundaryCommSubset phase,
//!                                std::vector<BoundaryVariable *> bvars_subset)
//! \brief initiate MPI_Irecv()

void BoundaryValues::StartReceivingSubset(BoundaryCommSubset phase,
                                          std::vector<BoundaryVariable *> bvars_subset) {
  for (auto bvars_it = bvars_subset.begin(); bvars_it != bvars_subset.end();
       ++bvars_it) {
    (*bvars_it)->StartReceiving(phase);
  }

  // KGF: begin shearing-box exclusive section of original StartReceivingForInit()
  // find send_block_id and recv_block_id;
  if (shearing_box != 0) {
    switch (phase) {
      case BoundaryCommSubset::mesh_init:
        break;
      case BoundaryCommSubset::radiation:
        break;
      case BoundaryCommSubset::radhydro:
      case BoundaryCommSubset::all:
      case BoundaryCommSubset::orbital:
        for (auto bvars_it = bvars_subset.begin(); bvars_it != bvars_subset.end();
             ++bvars_it) {
          (*bvars_it)->StartReceivingShear(phase);
        }
        break;
      case BoundaryCommSubset::gr_amr:
        // shearing box is currently incompatible with both GR and AMR
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues::StartReceiving" << std::endl
            << "BoundaryCommSubset::gr_amr was passed as the 'phase' argument while\n"
            << "SHEARING_BOX=1 is enabled. Shearing box calculations are currently\n"
            << "incompatible with both AMR and GR" << std::endl;
        ATHENA_ERROR(msg);
        break;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundary(BoundaryCommSubset phase,
//!                                         std::vector<BoundaryVariable *> bvars_subset)
//! \brief clean up the boundary flags after each loop
//!
//! \note
//! BoundaryCommSubset::mesh_init corresponds to initial exchange of conserved fluid
//! variables and magentic fields, while BoundaryCommSubset::gr_amr corresponds to fluid
//! primitive variables sent only in the case of GR with refinement

void BoundaryValues::ClearBoundarySubset(BoundaryCommSubset phase,
                                         std::vector<BoundaryVariable *> bvars_subset) {
  for (auto bvars_it = bvars_subset.begin(); bvars_it != bvars_subset.end();
       ++bvars_it) {
    (*bvars_it)->ClearBoundary(phase);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ApplyPhysicalBoundaries(const Real time, const Real dt,
//!                                        std::vector<BoundaryVariable *> bvars_subset)
//! \brief Apply all the physical boundary conditions for both hydro and field
//!
//! \note
//! - temporarily hardcode Hydro and Field access for coupling in EOS U(W) + calc bcc
//!   and when passed to user-defined boundary function stored in function pointer array

void BoundaryValues::ApplyPhysicalBoundaries(const Real time, const Real dt,
                                   std::vector<BoundaryVariable *> bvars_subset) {
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

  // KGF: COUPLING OF QUANTITIES (must be manually specified)
  // downcast BoundaryVariable ptrs to known derived class types: RTTI via dynamic_cast
  // HydroBoundaryVariable *phbvar =
  //     dynamic_cast<HydroBoundaryVariable *>(bvars_main_int[0]);
  Hydro *ph = pmb->phydro;

  //! \todo (felker):
  //! - passing nullptrs (pf) if no MHD (coarse_* no longer in MeshRefinement)
  //!   (may be fine to unconditionally directly set to pmb->pfield. See bvals_refine.cpp)

  // FaceCenteredBoundaryVariable *pfbvar = nullptr;
  Field *pf = nullptr;
  if (MAGNETIC_FIELDS_ENABLED) {
    pf = pmb->pfield;
    // pfbvar = dynamic_cast<FaceCenteredBoundaryVariable *>(bvars_main_int[1]);
  }
  PassiveScalars *ps = nullptr;
  if (NSCALARS > 0) {
    ps = pmb->pscalars;
  }

  NRRadiation *prad = nullptr;
  RadBoundaryVariable *pradbvar = nullptr;
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    prad = pmb->pnrrad;
    pradbvar = &(prad->rad_bvar);
  }

  CosmicRay *pcr = nullptr;
  //CellCenteredBoundaryVariable *pcrbvar = nullptr;
  if (CR_ENABLED) {
    pcr = pmb->pcr;
    //pcrbvar = &(pcr->cr_bvar);
  }

  // Apply boundary function on inner-x1 and update W,bcc (if not periodic)
  if (apply_bndry_fn_[BoundaryFace::inner_x1]) {
    DispatchBoundaryFunctions(pmb, pco, time, dt,
                              pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST,
                              ph->w, pf->b, prad->ir, pcr->u_cr,
                              BoundaryFace::inner_x1, bvars_subset);
    // KGF: COUPLING OF QUANTITIES (must be manually specified)
    if (MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                              pmb->is-NGHOST, pmb->is-1,
                                              bjs, bje, bks, bke);
    }
    pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                    pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
    if (NSCALARS > 0) {
      pmb->peos->PassiveScalarPrimitiveToConserved(
          ps->r, ph->u, ps->s, pco, pmb->is-NGHOST, pmb->is-1, bjs, bje, bks, bke);
    }
  }

  // Apply boundary function on outer-x1 and update W,bcc (if not periodic)
  if (apply_bndry_fn_[BoundaryFace::outer_x1]) {
    DispatchBoundaryFunctions(pmb, pco, time, dt,
                              pmb->is, pmb->ie, bjs, bje, bks, bke, NGHOST,
                              ph->w, pf->b, prad->ir, pcr->u_cr,
                              BoundaryFace::outer_x1, bvars_subset);
    // KGF: COUPLING OF QUANTITIES (must be manually specified)
    if (MAGNETIC_FIELDS_ENABLED) {
      pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                              pmb->ie+1, pmb->ie+NGHOST,
                                              bjs, bje, bks, bke);
    }
    pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                    pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
    if (NSCALARS > 0) {
      pmb->peos->PassiveScalarPrimitiveToConserved(
          ps->r, ph->u, ps->s, pco, pmb->ie+1, pmb->ie+NGHOST, bjs, bje, bks, bke);
    }
  }

  if (pmb->block_size.nx2 > 1) { // 2D or 3D
    // Apply boundary function on inner-x2 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::inner_x2]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, pmb->js, pmb->je, bks, bke, NGHOST,
                                ph->w, pf->b, prad->ir, pcr->u_cr,
                                BoundaryFace::inner_x2, bvars_subset);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, pmb->js-NGHOST, pmb->js-1,
                                                bks, bke);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
      if (NSCALARS > 0) {
        pmb->peos->PassiveScalarPrimitiveToConserved(
            ps->r, ph->u, ps->s, pco, bis, bie, pmb->js-NGHOST, pmb->js-1, bks, bke);
      }
    }

    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) &&
           (block_bcs[BoundaryFace::inner_x2] != BoundaryFlag::block)) {
      if (prad->rotate_theta == 1) {
        pradbvar->RotateHPi_InnerX2(time, dt, bis, bie, pmb->js, bks, bke, NGHOST);
      }
    }// end radiation

    // Apply boundary function on outer-x2 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::outer_x2]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, pmb->js, pmb->je, bks, bke, NGHOST,
                                ph->w, pf->b, prad->ir, pcr->u_cr,
                                BoundaryFace::outer_x2, bvars_subset);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, pmb->je+1, pmb->je+NGHOST,
                                                bks, bke);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
      if (NSCALARS > 0) {
        pmb->peos->PassiveScalarPrimitiveToConserved(
            ps->r, ph->u, ps->s, pco, bis, bie, pmb->je+1, pmb->je+NGHOST, bks, bke);
      }
    }
  }

  if (pmb->block_size.nx3 > 1) { // 3D
    bjs = pmb->js - NGHOST;
    bje = pmb->je + NGHOST;

    // Apply boundary function on inner-x3 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::inner_x3]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST,
                                ph->w, pf->b, prad->ir, pcr->u_cr,
                                BoundaryFace::inner_x3, bvars_subset);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, bjs, bje,
                                                pmb->ks-NGHOST, pmb->ks-1);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
      if (NSCALARS > 0) {
        pmb->peos->PassiveScalarPrimitiveToConserved(
            ps->r, ph->u, ps->s, pco, bis, bie, bjs, bje, pmb->ks-NGHOST, pmb->ks-1);
      }
    }

    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) &&
           (block_bcs[BoundaryFace::inner_x3] != BoundaryFlag::block)) {
      if (prad->rotate_phi == 1) {
        pradbvar->RotateHPi_InnerX3(time, dt, bis, bie, bjs, bje, pmb->ks, NGHOST);
      } else if (prad->rotate_phi == 2) {
        pradbvar->RotatePi_InnerX3(time, dt, bis, bie, bjs, bje, pmb->ks,NGHOST);
      }
    }

    // Apply boundary function on outer-x3 and update W,bcc (if not periodic)
    if (apply_bndry_fn_[BoundaryFace::outer_x3]) {
      DispatchBoundaryFunctions(pmb, pco, time, dt,
                                bis, bie, bjs, bje, pmb->ks, pmb->ke, NGHOST,
                                ph->w, pf->b, prad->ir, pcr->u_cr,
                                BoundaryFace::outer_x3, bvars_subset);
      // KGF: COUPLING OF QUANTITIES (must be manually specified)
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pfield->CalculateCellCenteredField(pf->b, pf->bcc, pco,
                                                bis, bie, bjs, bje,
                                                pmb->ke+1, pmb->ke+NGHOST);
      }
      pmb->peos->PrimitiveToConserved(ph->w, pf->bcc, ph->u, pco,
                                      bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
      if (NSCALARS > 0) {
        pmb->peos->PassiveScalarPrimitiveToConserved(
            ps->r, ph->u, ps->s, pco, bis, bie, bjs, bje, pmb->ke+1, pmb->ke+NGHOST);
      }
    }
    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) &&
           (block_bcs[BoundaryFace::outer_x3] != BoundaryFlag::block)) {
      if (prad->rotate_phi == 1) {
        pradbvar->RotateHPi_OuterX3(time, dt, bis, bie, bjs, bje, pmb->ke, NGHOST);
      } else if (prad->rotate_phi == 2) {
        pradbvar->RotatePi_OuterX3(time, dt, bis, bie, bjs, bje, pmb->ke, NGHOST);
      }
    }
  }
  return;
}

//! \brief
//!
//! \note
//! - should "bvars_it" be fixed in this class member function? Or passed as argument?
void BoundaryValues::DispatchBoundaryFunctions(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh,
    AthenaArray<Real> &prim, FaceField &b, AthenaArray<Real> &ir,
    AthenaArray<Real> &u_cr,  BoundaryFace face,
    std::vector<BoundaryVariable *> bvars_subset) {
  if (block_bcs[face] ==  BoundaryFlag::user) {  // user-enrolled BCs
    pmy_mesh_->BoundaryFunction_[face](pmb, pco, prim, b, time, dt,
                                       il, iu, jl, ju, kl, ku, NGHOST);

    // user-defined  boundary for radiation
    if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)) {
      pmy_mesh_->RadBoundaryFunc_[face](pmb,pco,pmb->pnrrad,prim,b, ir,time,dt,
                                             il,iu,jl,ju,kl,ku,NGHOST);
    }

    if (CR_ENABLED) {
      pmy_mesh_->CRBoundaryFunc_[face](pmb,pco,pmb->pcr,prim, b, u_cr,time,dt,
                                              il,iu,jl,ju,kl,ku,NGHOST);
    }
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
  for (auto bvars_it = bvars_subset.begin(); bvars_it != bvars_subset.end();
       ++bvars_it) {
    switch (block_bcs[face]) {
      case BoundaryFlag::user: // handled above, outside loop over BoundaryVariable objs
        break;
      case BoundaryFlag::reflect:
        switch (face) {
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
        switch (face) {
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
      case BoundaryFlag::vacuum: // special boundary condition type for radiation
        switch (face) {
          case BoundaryFace::undef:
            ATHENA_ERROR(msg);
          case BoundaryFace::inner_x1:
            (*bvars_it)->VacuumInnerX1(time, dt, il, jl, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x1:
            (*bvars_it)->VacuumOuterX1(time, dt, iu, jl, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::inner_x2:
            (*bvars_it)->VacuumInnerX2(time, dt, il, iu, jl, kl, ku, NGHOST);
            break;
          case BoundaryFace::outer_x2:
            (*bvars_it)->VacuumOuterX2(time, dt, il, iu, ju, kl, ku, NGHOST);
            break;
          case BoundaryFace::inner_x3:
            (*bvars_it)->VacuumInnerX3(time, dt, il, iu, jl, ju, kl, NGHOST);
            break;
          case BoundaryFace::outer_x3:
            (*bvars_it)->VacuumOuterX3(time, dt, il, iu, jl, ju, ku, NGHOST);
            break;
        }
        break;
      case BoundaryFlag::polar_wedge:
        switch (face) {
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
                      << "MeshBlock boundary other than inner x2 or outer x2"
                      << std::endl;
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
    } // end switch (block_bcs[face])
  } // end loop over BoundaryVariable *
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ComputeShear(const Real time_fc, const Real time_int)
//! \brief Calculate the following quantities:
//! send_gid recv_gid send_lid recv_lid send_rank recv_rank,
//! send_size_hydro  recv_size_hydro: for MPI_Irecv
//! eps_,joverlap_: for update the conservative

void BoundaryValues::ComputeShear(const Real time_fc, const Real time_int) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  Mesh *pmesh = pmb->pmy_mesh;

  if (shearing_box == 1) {
    int nx2 = pmb->block_size.nx2;
    int js = pmb->js; int je = pmb->je;
    int jl = js-NGHOST; int ju = je+NGHOST;
    int level = pmb->loc.level - pmesh->root_level;

    // Update the amount of shear:
    Real x2size = pmesh->mesh_size.x2max - pmesh->mesh_size.x2min;
    Real dx = x2size/static_cast<Real>(nblx2*nx2);
    Real yshear, deltay;
    int joffset, Ngrids;

    // flux
    yshear = qomL_*time_fc;
    deltay = std::fmod(yshear, x2size);
    joffset = static_cast<int>(deltay/dx); // assumes uniform grid in azimuth
    Ngrids  = static_cast<int>(joffset/nx2);
    joverlap_flux_   = joffset - Ngrids*nx2;
    eps_flux_ = (std::fmod(deltay, dx))/dx;

    for (int upper=0; upper<2; upper++) {
      if (is_shear[upper]) {
        int *counts1 = sb_flux_data_[upper].send_count;
        int *counts2 = sb_flux_data_[upper].recv_count;
        SimpleNeighborBlock *nb1 = sb_flux_data_[upper].send_neighbor;
        SimpleNeighborBlock *nb2 = sb_flux_data_[upper].recv_neighbor;
        int *jmin1 = sb_flux_data_[upper].jmin_send;
        int *jmax1 = sb_flux_data_[upper].jmax_send;
        int *jmin2 = sb_flux_data_[upper].jmin_recv;
        int *jmax2 = sb_flux_data_[upper].jmax_recv;
        int jo     = (1-2*upper)*((1-upper)+joverlap_flux_);
        int Ng     = (1-2*upper)*Ngrids;

        for (int n=0; n<3; n++) {
          nb1[n].gid  = -1;  nb1[n].lid   = -1;
          nb1[n].rank = -1;  nb1[n].level = -1;
          nb2[n].gid  = -1;  nb2[n].lid   = -1;
          nb2[n].rank = -1;  nb2[n].level = -1;
          counts1[n]  = 0;   counts2[n]   = 0;
          jmin1[n]    = 0;   jmax1[n]     = 0;
          jmin2[n]    = 0;   jmax2[n]     = 0;
        }

        int jblock = 0;
        for (int j=0; j<nblx2; j++) {
          // find global index of current MeshBlock on the shearing boundary block list
          if (shbb_[upper][j].gid == pmb->gid) jblock = j;
        }
        int bshift;
        // send js+jo : je+jo
        // recv js-xgh:je+1+xgh
        // js+jo<=je+1+xgh-2*nx2
        //    jo<=xgh-nx2

        if (jo<=xgh_-nx2) {
          // case 1
          // send from this to 0: jb+(Ng-2), 1: jb+(Ng-1), 2: jb+(Ng  )
          // recv to this from 0: jb-(Ng-2), 1: jb-(Ng-1), 2: jb-(Ng  )
          bshift = -2;
        } else if (jo<=xgh_) {
          // case 2
          // send from this to 0: jb+(Ng-1), 1: jb+(Ng  ), 2: jb+(Ng+1)
          // recv to this from 0: jb-(Ng-1), 1: jb-(Ng  ), 2: jb-(Ng+1)
          bshift = -1;
        } else if (jo<=xgh_+nx2) {
          // case 3
          // send from this to 0: jb+(Ng  ), 1: jb+(Ng+1), 2: jb+(Ng+2)
          // recv to this from 0: jb-(Ng  ), 1: jb-(Ng+1), 2: jb-(Ng+2)
          bshift = 0;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in BoundaryValues::ComputeShear" << std::endl;
          ATHENA_ERROR(msg);
        }
        jmin1[0] = js+jo;
        jmax1[0] = je+1+xgh_+bshift*nx2;
        jmin1[1] = std::max(js+jo,js-xgh_+(bshift+1)*nx2);
        jmax1[1] = std::min(je+jo,je+xgh_+1+(bshift+1)*nx2);
        jmin1[2] = js-xgh_+(bshift+2)*nx2;
        jmax1[2] = je+jo;
        for (int n=0; n<3; n++) {
          counts1[n] = jmax1[n]-jmin1[n]+1;
          if (counts1[n]>0) {
            std::int64_t jtmp = (jblock+(Ng+bshift+n))%nblx2;
            if (jtmp <  0)     jtmp += nblx2;
            nb1[n].gid   = shbb_[upper][jtmp].gid;
            nb1[n].lid   = shbb_[upper][jtmp].lid;
            nb1[n].rank  = shbb_[upper][jtmp].rank;
            nb1[n].level = shbb_[upper][jtmp].level;
          }
        }
        jmin2[0] = js+jo-bshift*nx2;
        jmax2[0] = je+1+xgh_;
        jmin2[1] = std::max(js-xgh_,js+jo-(bshift+1)*nx2);
        jmax2[1] = std::min(je+1+xgh_,je+jo-(bshift+1)*nx2);
        jmin2[2] = js-xgh_;
        jmax2[2] = je+jo-(bshift+2)*nx2;
        for (int n=0; n<3; n++) {
          counts2[n] = jmax2[n]-jmin2[n]+1;
          if (counts2[n]>0) {
            std::int64_t jtmp = (jblock-(Ng+bshift+n))%nblx2;
            if (jtmp <  0)     jtmp += nblx2;
            nb2[n].gid   = shbb_[upper][jtmp].gid;
            nb2[n].lid   = shbb_[upper][jtmp].lid;
            nb2[n].rank  = shbb_[upper][jtmp].rank;
            nb2[n].level = shbb_[upper][jtmp].level;
          }
        }
      }
    } // end loop over inner, outer shearing boundaries

    // after integration
    yshear = qomL_*time_int;
    deltay = std::fmod(yshear, x2size);
    joffset = static_cast<int>(deltay/dx); // assumes uniform grid in azimuth
    Ngrids  = static_cast<int>(joffset/nx2);
    joverlap_   = joffset - Ngrids*nx2;
    eps_ = (std::fmod(deltay, dx))/dx;
    for (int upper=0; upper<2; upper++) {
      if (is_shear[upper]) {
        int *counts1 = sb_data_[upper].send_count;
        int *counts2 = sb_data_[upper].recv_count;
        SimpleNeighborBlock *nb1 = sb_data_[upper].send_neighbor;
        SimpleNeighborBlock *nb2 = sb_data_[upper].recv_neighbor;
        int *jmin1 = sb_data_[upper].jmin_send;
        int *jmax1 = sb_data_[upper].jmax_send;
        int *jmin2 = sb_data_[upper].jmin_recv;
        int *jmax2 = sb_data_[upper].jmax_recv;
        int jo     = (1-2*upper)*((1-upper)+joverlap_);
        int Ng     = (1-2*upper)*Ngrids;

        for (int n=0; n<4; n++) {
          nb1[n].gid  = -1;  nb1[n].lid   = -1;
          nb1[n].rank = -1;  nb1[n].level = -1;
          nb2[n].gid  = -1;  nb2[n].lid   = -1;
          nb2[n].rank = -1;  nb2[n].level = -1;
          counts1[n]  = 0;   counts2[n]   = 0;
          jmin1[n]    = 0;   jmax1[n]     = 0;
          jmin2[n]    = 0;   jmax2[n]     = 0;
        }

        int jblock = 0;
        for (int j=0; j<nblx2; j++) {
          // find global index of current MeshBlock on the shearing boundary block list
          if (shbb_[upper][j].gid == pmb->gid) jblock = j;
        }
        int bshift;
        if (jo<=xgh_+NGHOST-2*nx2) {
          // case 1
          // send from this to 0: jb+(Ng-3), 1: jb+(Ng-2), 2: jb+(Ng-1), 3: jb+(Ng)
          // recv to this from 0: jb-(Ng-3), 1: jb-(Ng-2), 2: jb-(Ng-1), 3: jb-(Ng)
          bshift = -3;
        } else if (jo<=xgh_+NGHOST-nx2) {
          // case 2
          // send from this to 0: jb+(Ng-2), 1: jb+(Ng-1), 2: jb+(Ng), 3: jb+(Ng+1)
          // recv to this from 0: jb-(Ng-2), 1: jb-(Ng-1), 2: jb-(Ng), 3: jb-(Ng+1)
          bshift = -2;
        } else if (jo<=xgh_+NGHOST) {
          // case 3
          // send from this to 0: jb+(Ng-1), 1: jb+(Ng), 2: jb+(Ng+1), 3: jb+(Ng+2)
          // recv to this from 0: jb-(Ng-1), 1: jb-(Ng), 2: jb-(Ng+1), 3: jb-(Ng+2)
          bshift = -1;
        } else if (jo<=xgh_+NGHOST+nx2) {
          // case 4
          // send from this to 0: jb+(Ng), 1: jb+(Ng+1), 2: jb+(Ng+2), 3: jb+(Ng+3)
          // recv to this from 0: jb-(Ng), 1: jb-(Ng+1), 2: jb-(Ng+2), 3: jb-(Ng+3)
          bshift = 0;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in BoundaryValues::ComputeShear" << std::endl;
          ATHENA_ERROR(msg);
        }
        jmin1[0] = js+jo;
        jmax1[0] = ju+1+xgh_+bshift*nx2;
        jmin1[1] = std::max(js+jo,jl-xgh_+(bshift+1)*nx2);
        jmax1[1] = std::min(je+jo,ju+xgh_+1+(bshift+1)*nx2);
        jmin1[2] = std::max(js+jo,jl-xgh_+(bshift+2)*nx2);
        jmax1[2] = std::min(je+jo,ju+xgh_+1+(bshift+2)*nx2);
        jmin1[3] = jl-xgh_+(bshift+3)*nx2;
        jmax1[3] = je+jo;
        for (int n=0; n<4; n++) {
          counts1[n] = jmax1[n]-jmin1[n]+1;
          if (counts1[n]>0) {
            std::int64_t jtmp = (jblock+(Ng+bshift+n))%nblx2;
            if (jtmp <  0)     jtmp += nblx2;
            nb1[n].gid   = shbb_[upper][jtmp].gid;
            nb1[n].lid   = shbb_[upper][jtmp].lid;
            nb1[n].rank  = shbb_[upper][jtmp].rank;
            nb1[n].level = shbb_[upper][jtmp].level;
          }
        }
        jmin2[0] = js+jo-bshift*nx2;
        jmax2[0] = ju+1+xgh_;
        jmin2[1] = std::max(jl-xgh_,js+jo-(bshift+1)*nx2);
        jmax2[1] = std::min(ju+1+xgh_,je+jo-(bshift+1)*nx2);
        jmin2[2] = std::max(jl-xgh_,js+jo-(bshift+2)*nx2);
        jmax2[2] = std::min(ju+1+xgh_,je+jo-(bshift+2)*nx2);
        jmin2[3] = jl-xgh_;
        jmax2[3] = je+jo-(bshift+3)*nx2;
        for (int n=0; n<4; n++) {
          counts2[n] = jmax2[n]-jmin2[n]+1;
          if (counts2[n]>0) {
            std::int64_t jtmp = (jblock-(Ng+bshift+n))%nblx2;
            if (jtmp <  0)     jtmp += nblx2;
            nb2[n].gid   = shbb_[upper][jtmp].gid;
            nb2[n].lid   = shbb_[upper][jtmp].lid;
            nb2[n].rank  = shbb_[upper][jtmp].rank;
            nb2[n].level = shbb_[upper][jtmp].level;
          }
        }
      }
    } // end loop over inner, outer shearing boundaries
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \brief Public function, to be called in MeshBlock ctor for keeping MPI tag bitfields
//! consistent across MeshBlocks, even if certain MeshBlocks only construct a subset of
//! physical variable classes

int BoundaryValues::AdvanceCounterPhysID(int num_phys) {
#ifdef MPI_PARALLEL
  //! \todo (felker):
  //! * add safety checks? input, output are positive, obey <= 31= MAX_NUM_PHYS
  int start_id = bvars_next_phys_id_;
  bvars_next_phys_id_ += num_phys;
  return start_id;
#else
  return 0;
#endif
}
