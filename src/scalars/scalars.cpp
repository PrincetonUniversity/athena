//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file scalars.cpp
//  \brief implementation of functions in class PassiveScalars

// C headers

// C++ headers
#include <algorithm>
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "scalars.hpp"

// constructor, initializes data structures and parameters

PassiveScalars::PassiveScalars(MeshBlock *pmb, ParameterInput *pin)  :
    s(NSCALARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    s_flux{ {NSCALARS, pmb->ncells3, pmb->ncells2, pmb->ncells1+1},
            {NSCALARS, pmb->ncells3, pmb->ncells2+1, pmb->ncells1,
             (pmb->pmy_mesh->f2_ ? AthenaArray<Real>::DataStatus::allocated :
              AthenaArray<Real>::DataStatus::empty)},
            {NSCALARS, pmb->ncells3+1, pmb->ncells2, pmb->ncells1,
             (pmb->pmy_mesh->f3_ ? AthenaArray<Real>::DataStatus::allocated :
              AthenaArray<Real>::DataStatus::empty)}
    },
    coarse_s_(NSCALARS, pmb->ncc3, pmb->ncc2, pmb->ncc1,
              (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
               AthenaArray<Real>::DataStatus::empty)),
    sbvar(pmb, &s, &coarse_s_, s_flux),
    pmy_block(pmb) {
  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, ncells3 = pmb->ncells3;
  Mesh *pm = pmy_block->pmy_mesh;

  pmb->RegisterMeshBlockData(s);

  // Allocate optional passive scalar variable memory registers for time-integrator
  if (pmb->precon->xorder == 4) {
    // fourth-order cell-centered approximations
    s_cc.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  }

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel == true) { 
    pmy_block->pmr->AddToRefinement(&s, &coarse_s_);
  }

  // enroll CellCenteredBoundaryVariable object
  sbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&sbvar);
  pmb->pbval->bvars_main_int.push_back(&sbvar);

  // Allocate memory for scratch arrays
  dxw_.NewAthenaArray(ncells1);
  wl_.NewAthenaArray(NSCALARS, ncells1);
  wr_.NewAthenaArray(NSCALARS, ncells1);
  wlb_.NewAthenaArray(NSCALARS, ncells1);
  x1face_area_.NewAthenaArray(ncells1+1);
  if (pm->f2_) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if (pm->f3_) {
    x3face_area_.NewAthenaArray(ncells1);
    x3face_area_p1_.NewAthenaArray(ncells1);
  }
  cell_volume_.NewAthenaArray(ncells1);

  // fourth-order 4D scratch arrays
  wl3d_.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  wr3d_.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  scr1_nkji_.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  laplacian_l_fc_.NewAthenaArray(ncells1);
  laplacian_r_fc_.NewAthenaArray(ncells1);
}
