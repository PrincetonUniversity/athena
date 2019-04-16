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
    sbvar(pmb, &s, &coarse_s_, s_flux),
    pmy_block(pmb) {
  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, ncells3 = pmb->ncells3;
  Mesh *pm = pmy_block->pmy_mesh;

  // Allocate optional passive scalar variable memory registers for time-integrator
  if (pmb->precon->xorder == 4) {
    // fourth-order cell-centered approximations
    s_cc.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  }

  s_flux[X1DIR].NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1+1);
  if (pm->f2_)
    s_flux[X2DIR].NewAthenaArray(NSCALARS, ncells3, ncells2+1, ncells1);
  if (pm->f3_)
    s_flux[X3DIR].NewAthenaArray(NSCALARS, ncells3+1, ncells2, ncells1);

  // allocate prolongation buffers
  if (pm->multilevel == true) {
    int ncc1 = pmb->block_size.nx1/2 + 2*NGHOST;
    int ncc2 = 1;
    if (pm->f2_) ncc2 = pmb->block_size.nx2/2 + 2*NGHOST;
    int ncc3 = 1;
    if (pm->f3_) ncc3 = pmb->block_size.nx3/2 + 2*NGHOST;
    coarse_s_.NewAthenaArray(NSCALARS, ncc3, ncc2, ncc1);
    // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
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

  // fourth-order hydro
  // 4D scratch arrays
  wl3d_.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  wr3d_.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  scr1_nkji_.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  laplacian_l_fc_.NewAthenaArray(ncells1);
  laplacian_r_fc_.NewAthenaArray(ncells1);
}
