//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.cpp
//! \brief implementation of functions in class Gravity

// C headers

// C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../bvals/bvals_interfaces.hpp"
#include "../bvals/cc/bvals_cc.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"
#include "gravity.hpp"
#include "mg_gravity.hpp"

//! constructor, initializes data structures and parameters
// TODO(felker): change "MeshBlock *pmb" to reference member, set in initializer list
//----------------------------------------------------------------------------------------
//! \fn Gravity::Gravity(MeshBlock *pmb, ParameterInput *pin)
//! \brief Gravity constructor
Gravity::Gravity(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb), phi(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    coarse_phi(NHYDRO, pmb->ncc3, pmb->ncc2, pmb->ncc1,
              (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
               AthenaArray<Real>::DataStatus::empty)),
    empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
    four_pi_G(pmb->pmy_mesh->four_pi_G_),
    output_defect(false), fill_ghost(false),
    gbvar(pmb, &phi, &coarse_phi, empty_flux, false) {
  if (four_pi_G == 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Gravity::Gravity" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  if (SELF_GRAVITY_ENABLED == 2) {
    output_defect = pin->GetOrAddBoolean("gravity", "output_defect", false);
    if (output_defect)
      def.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);
    if (pmb->pmy_mesh->multilevel) {
      fill_ghost = pin->GetOrAddBoolean("gravity", "fill_ghost", false);
      if (fill_ghost) {
        int nf1 = pmb->block_size.nx2*pmb->block_size.nx3;
        int nf2 = pmb->block_size.nx1*pmb->block_size.nx3;
        int nf3 = pmb->block_size.nx1*pmb->block_size.nx2;
        fbuf_[inner_x1].NewAthenaArray(nf1);
        fbuf_[outer_x1].NewAthenaArray(nf1);
        fbuf_[inner_x2].NewAthenaArray(nf2);
        fbuf_[outer_x2].NewAthenaArray(nf2);
        fbuf_[inner_x3].NewAthenaArray(nf3);
        fbuf_[outer_x3].NewAthenaArray(nf3);
      }
    }
  }

  // using Gravity as an example of: containing full object members instead of pointer
  // memebers, construting BoundaryVariable composite obj (no default ctor) in Gravity
  // ctor initializer list, avoiding dynamically-managed memory and the need for a
  // user-provided dtor.

  // Enroll CellCenteredBoundaryVariable object
  gbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&gbvar);
  pmb->pbval->pgbvar = &gbvar;
}


//----------------------------------------------------------------------------------------
//! \fn Gravity::SaveFaceBoundaries()
//! \brief Save face boundary values for multigrid + mesh refinement
void Gravity::SaveFaceBoundaries() {
  int mylevel = pmy_block->pbval->nblevel[1][1][1];
  int is = pmy_block->is, ie = pmy_block->ie,
      js = pmy_block->js, je = pmy_block->je,
      ks = pmy_block->ks, ke = pmy_block->ke;

  if (pmy_block->pbval->nblevel[1][1][0] < mylevel) {
    int p = 0;
    BufferUtility::PackData(phi, fbuf_[inner_x1].data(), 0, 0,
                            is-1, is-1, js,   je,   ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[1][1][2] < mylevel) {
    int p = 0;
    BufferUtility::PackData(phi, fbuf_[outer_x1].data(), 0, 0,
                            is+1, is+1, js,   je,   ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[1][0][1] < mylevel) {
    int p = 0;
    BufferUtility::PackData(phi, fbuf_[inner_x2].data(), 0, 0,
                            is,   ie,   js-1, js-1, ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[1][2][1] < mylevel) {
    int p = 0;
    BufferUtility::PackData(phi, fbuf_[outer_x2].data(), 0, 0,
                            is,   ie,   je+1, je+1, ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[0][1][1] < mylevel) {
    int p = 0;
    BufferUtility::PackData(phi, fbuf_[inner_x3].data(), 0, 0,
                            is,   ie,   js,   je,   ks-1, ks-1, p);
  }
  if (pmy_block->pbval->nblevel[2][1][1] < mylevel) {
    int p = 0;
    BufferUtility::PackData(phi, fbuf_[outer_x3].data(), 0, 0,
                            is,   ie,   js,   je,   ke+1, ke+1, p);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn Gravity::RestoreFaceBoundaries()
//! \brief Restore face boundary values for multigrid + mesh refinement
void Gravity::RestoreFaceBoundaries() {
  int mylevel = pmy_block->pbval->nblevel[1][1][1];
  int p = 0;
  int is = pmy_block->is, ie = pmy_block->ie,
      js = pmy_block->js, je = pmy_block->je,
      ks = pmy_block->ks, ke = pmy_block->ke;

  if (pmy_block->pbval->nblevel[1][1][0] < mylevel) {
    int p = 0;
    BufferUtility::UnpackData(fbuf_[inner_x1].data(), phi, 0, 0,
                              is-1, is-1, js,   je,   ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[1][1][2] < mylevel) {
    int p = 0;
    BufferUtility::UnpackData(fbuf_[outer_x1].data(), phi, 0, 0,
                              is+1, is+1, js,   je,   ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[1][0][1] < mylevel) {
    int p = 0;
    BufferUtility::UnpackData(fbuf_[inner_x2].data(), phi, 0, 0,
                              is,   ie,   js-1, js-1, ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[1][2][1] < mylevel) {
    int p = 0;
    BufferUtility::UnpackData(fbuf_[outer_x2].data(), phi, 0, 0,
                              is,   ie,   je+1, je+1, ks,   ke,   p);
  }
  if (pmy_block->pbval->nblevel[0][1][1] < mylevel) {
    int p = 0;
    BufferUtility::UnpackData(fbuf_[inner_x3].data(), phi, 0, 0,
                              is,   ie,   js,   je,   ks-1, ks-1, p);
  }
  if (pmy_block->pbval->nblevel[2][1][1] < mylevel) {
    int p = 0;
    BufferUtility::UnpackData(fbuf_[outer_x3].data(), phi, 0, 0,
                              is,   ie,   js,   je,   ke+1, ke+1, p);
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn Gravity::ExpandPhysicalBoundaries()
//! \brief Expand physical boundary values to NGHOST = 2 and to edges/corners.
void Gravity::ExpandPhysicalBoundaries() {
  int is = pmy_block->is, ie = pmy_block->ie,
      js = pmy_block->js, je = pmy_block->je,
      ks = pmy_block->ks, ke = pmy_block->ke;

  // push face boundary values
  if (pmy_block->pbval->nblevel[1][1][0] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++)
        phi(k, j, is-2) = phi(k, j, is-1);
    }
  }
  if (pmy_block->pbval->nblevel[1][1][2] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++)
        phi(k, j, ie+2) = phi(k, j, ie+1);
    }
  }
  if (pmy_block->pbval->nblevel[1][0][1] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int i = is; i <= ie; i++)
        phi(k, js-2, i) = phi(k, js-1, i);
    }
  }
  if (pmy_block->pbval->nblevel[1][2][1] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int i = is; i <= ie; i++)
        phi(k, je+2, i) = phi(k, je+1, i);
    }
  }
  if (pmy_block->pbval->nblevel[0][1][1] < 0) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++)
        phi(ks-2, j, i) = phi(ks-1, j, i);
    }
  }
  if (pmy_block->pbval->nblevel[2][1][1] < 0) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++)
        phi(ke+2, j, i) = phi(ke+1, j, i);
    }
  }

  // fill edges
  if (pmy_block->pbval->nblevel[1][0][0] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(phi(k, js-1, is) + phi(k, js, is-1));
      phi(k, js-1, is-1) = p;
      phi(k, js-1, is-2) = p;
      phi(k, js-2, is-1) = p;
      phi(k, js-2, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[1][0][2] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(phi(k, js-1, ie) + phi(k, js, ie+1));
      phi(k, js-1, ie+1) = p;
      phi(k, js-1, ie+2) = p;
      phi(k, js-2, ie+1) = p;
      phi(k, js-2, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[1][2][0] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(phi(k, je+1, is) + phi(k, je, is-1));
      phi(k, je+1, is-1) = p;
      phi(k, je+1, is-2) = p;
      phi(k, je+2, is-1) = p;
      phi(k, je+2, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[1][2][2] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(phi(k, je+1, ie) + phi(k, je, ie+1));
      phi(k, je+1, ie+1) = p;
      phi(k, je+1, ie+2) = p;
      phi(k, je+2, ie+1) = p;
      phi(k, je+2, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][1][0] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(phi(ks-1, j, is) + phi(ks, j, is-1));
      phi(ks-1, j, is-1) = p;
      phi(ks-1, j, is-2) = p;
      phi(ks-2, j, is-1) = p;
      phi(ks-2, j, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][1][2] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(phi(ks-1, j, ie) + phi(ks, j, ie+1));
      phi(ks-1, j, ie+1) = p;
      phi(ks-1, j, ie+2) = p;
      phi(ks-2, j, ie+1) = p;
      phi(ks-2, j, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][1][0] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(phi(ke+1, j, is) + phi(ke, j, is-1));
      phi(ke+1, j, is-1) = p;
      phi(ke+1, j, is-2) = p;
      phi(ke+2, j, is-1) = p;
      phi(ke+2, j, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][1][2] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(phi(ke+1, j, ie) + phi(ke, j, ie+1));
      phi(ke+1, j, ie+1) = p;
      phi(ke+1, j, ie+2) = p;
      phi(ke+2, j, ie+1) = p;
      phi(ke+2, j, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][0][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(phi(ks-1, js, i) + phi(ks, js-1, i));
      phi(ks-1, js-1, i) = p;
      phi(ks-1, js-2, i) = p;
      phi(ks-2, js-1, i) = p;
      phi(ks-2, js-2, i) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][2][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(phi(ks-1, je, i) + phi(ks, je+1, i));
      phi(ks-1, je+1, i) = p;
      phi(ks-1, je+2, i) = p;
      phi(ks-2, je+1, i) = p;
      phi(ks-2, je+2, i) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][0][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(phi(ke+1, js, i) + phi(ke, js-1, i));
      phi(ke+1, js-1, i) = p;
      phi(ke+1, js-2, i) = p;
      phi(ke+2, js-1, i) = p;
      phi(ke+2, js-2, i) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][2][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(phi(ke+1, je, i) + phi(ke, je+1, i));
      phi(ke+1, je+1, i) = p;
      phi(ke+1, je+2, i) = p;
      phi(ke+2, je+1, i) = p;
      phi(ke+2, je+2, i) = p;
    }
  }

  // fill corners
  if (pmy_block->pbval->nblevel[0][0][0] < 0) {
    Real p = (phi(ks, js-1, is-1) + phi(ks-1, js, is-1) + phi(ks-1, js-1, is))/3.0;
    phi(ks-1, js-1, is-1) = p;
    phi(ks-1, js-1, is-2) = p;
    phi(ks-1, js-2, is-1) = p;
    phi(ks-1, js-2, is-2) = p;
    phi(ks-2, js-1, is-1) = p;
    phi(ks-2, js-1, is-2) = p;
    phi(ks-2, js-2, is-1) = p;
    phi(ks-2, js-2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[0][0][2] < 0) {
    Real p = (phi(ks, js-1, ie+1) + phi(ks-1, js, ie+1) + phi(ks-1, js-1, ie))/3.0;
    phi(ks-1, js-1, ie+1) = p;
    phi(ks-1, js-1, ie+2) = p;
    phi(ks-1, js-2, ie+1) = p;
    phi(ks-1, js-2, ie+2) = p;
    phi(ks-2, js-1, ie+1) = p;
    phi(ks-2, js-1, ie+2) = p;
    phi(ks-2, js-2, ie+1) = p;
    phi(ks-2, js-2, ie+2) = p;
  }
  if (pmy_block->pbval->nblevel[0][2][0] < 0) {
    Real p = (phi(ks, je+1, is-1) + phi(ks-1, je, is-1) + phi(ks-1, je+1, is))/3.0;
    phi(ks-1, je+1, is-1) = p;
    phi(ks-1, je+1, is-2) = p;
    phi(ks-1, je+2, is-1) = p;
    phi(ks-1, je+2, is-2) = p;
    phi(ks-2, je+1, is-1) = p;
    phi(ks-2, je+1, is-2) = p;
    phi(ks-2, je+2, is-1) = p;
    phi(ks-2, je+2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[2][0][0] < 0) {
    Real p = (phi(ke, js-1, is-1) + phi(ke+1, js, is-1) + phi(ke+1, js-1, is))/3.0;
    phi(ke+1, js-1, is-1) = p;
    phi(ke+1, js-1, is-2) = p;
    phi(ke+1, js-2, is-1) = p;
    phi(ke+1, js-2, is-2) = p;
    phi(ke+2, js-1, is-1) = p;
    phi(ke+2, js-1, is-2) = p;
    phi(ke+2, js-2, is-1) = p;
    phi(ke+2, js-2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[0][2][2] < 0) {
    Real p = (phi(ks, je+1, ie+1) + phi(ks-1, je, ie+1) + phi(ks-1, je+1, ie))/3.0;
    phi(ks-1, je+1, ie+1) = p;
    phi(ks-1, je+1, ie+2) = p;
    phi(ks-1, je+2, ie+1) = p;
    phi(ks-1, je+2, ie+2) = p;
    phi(ks-2, je+1, ie+1) = p;
    phi(ks-2, je+1, ie+2) = p;
    phi(ks-2, je+2, ie+1) = p;
    phi(ks-2, je+2, ie+2) = p;
  }
  if (pmy_block->pbval->nblevel[2][0][2] < 0) {
    Real p = (phi(ke, js-1, ie+1) + phi(ke+1, js, ie+1) + phi(ke+1, js-1, ie))/3.0;
    phi(ke+1, js-1, ie+1) = p;
    phi(ke+1, js-1, ie+2) = p;
    phi(ke+1, js-2, ie+1) = p;
    phi(ke+1, js-2, ie+2) = p;
    phi(ke+2, js-1, ie+1) = p;
    phi(ke+2, js-1, ie+2) = p;
    phi(ke+2, js-2, ie+1) = p;
    phi(ke+2, js-2, ie+2) = p;
  }
  if (pmy_block->pbval->nblevel[2][2][0] < 0) {
    Real p = (phi(ke, je+1, is-1) + phi(ke+1, je, is-1) + phi(ke+1, je+1, is))/3.0;
    phi(ke+1, je+1, is-1) = p;
    phi(ke+1, je+1, is-2) = p;
    phi(ke+1, je+2, is-1) = p;
    phi(ke+1, je+2, is-2) = p;
    phi(ke+2, je+1, is-1) = p;
    phi(ke+2, je+1, is-2) = p;
    phi(ke+2, je+2, is-1) = p;
    phi(ke+2, je+2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[2][2][2] < 0) {
    Real p = (phi(ke, je+1, ie+1) + phi(ke+1, je, ie+1) + phi(ke+1, je+1, ie))/3.0;
    phi(ke+1, je+1, ie+1) = p;
    phi(ke+1, je+1, ie+2) = p;
    phi(ke+1, je+2, ie+1) = p;
    phi(ke+1, je+2, ie+2) = p;
    phi(ke+2, je+1, ie+1) = p;
    phi(ke+2, je+1, ie+2) = p;
    phi(ke+2, je+2, ie+1) = p;
    phi(ke+2, je+2, ie+2) = p;
  }

  return;
}

