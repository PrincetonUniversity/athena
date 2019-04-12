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
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "scalars.hpp"

// constructor, initializes data structures and parameters

PassiveScalars::PassiveScalars(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for primitive/conserved variables
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  s.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  pmb->RegisterMeshBlockData(s);

  // fourth-order cell-centered approximations
  s_cc.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);

  s_flux[X1DIR].NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1+1);
  if (pmy_block->block_size.nx2 > 1)
    s_flux[X2DIR].NewAthenaArray(NSCALARS, ncells3, ncells2+1, ncells1);
  if (pmy_block->block_size.nx3 > 1)
    s_flux[X3DIR].NewAthenaArray(NSCALARS, ncells3+1, ncells2, ncells1);

  // allocate prolongation buffers
  if (pmy_block->pmy_mesh->multilevel == true) {
    int ncc1 = pmb->block_size.nx1/2 + 2*NGHOST;
    int ncc2 = 1;
    if (pmb->block_size.nx2 > 1) ncc2 = pmb->block_size.nx2/2 + 2*NGHOST;
    int ncc3 = 1;
    if (pmb->block_size.nx3 > 1) ncc3 = pmb->block_size.nx3/2 + 2*NGHOST;
    coarse_s_.NewAthenaArray(NSCALARS, ncc3, ncc2, ncc1);
    // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
    pmy_block->pmr->AddToRefinement(&s, &coarse_s_);
  }

  // TODO(KGF): change to RAII
  // create object to interface with BoundaryValues
  psbval  = new CellCenteredBoundaryVariable(pmy_block, &s, &coarse_s_, s_flux);
  psbval->bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(psbval);

  pmb->pbval->bvars_main_int.push_back(psbval);

  // Allocate memory for scratch arrays
  dxw_.NewAthenaArray(ncells1);
  wl_.NewAthenaArray(NSCALARS, ncells1);
  wr_.NewAthenaArray(NSCALARS, ncells1);
  wlb_.NewAthenaArray(NSCALARS, ncells1);
  x1face_area_.NewAthenaArray(ncells1+1);
  if (pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if (pmy_block->block_size.nx3 > 1) {
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

// destructor

PassiveScalars::~PassiveScalars() {
  s.DeleteAthenaArray();

  // fourth-order integrator
  s_cc.DeleteAthenaArray();

  s_flux[X1DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) s_flux[X2DIR].DeleteAthenaArray();
  if (pmy_block->block_size.nx3 > 1) s_flux[X3DIR].DeleteAthenaArray();

  dxw_.DeleteAthenaArray();
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
  wlb_.DeleteAthenaArray();
  x1face_area_.DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if (pmy_block->block_size.nx3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();

  // fourth-order integrator
  // 4D scratch arrays
  wl3d_.DeleteAthenaArray();
  wr3d_.DeleteAthenaArray();
  scr1_nkji_.DeleteAthenaArray();
  laplacian_l_fc_.DeleteAthenaArray();
  laplacian_r_fc_.DeleteAthenaArray();

  delete psbval;
}
