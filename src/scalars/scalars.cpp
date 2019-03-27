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
#include "./scalars.hpp"

// constructor, initializes data structures and parameters

PassiveScalars::PassiveScalars(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for primitive/conserved variables
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  // Allocate memory registers for primitive/conserved variables for time-integrator
  s.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);

  // fourth-order hydro cell-centered approximations
  s_cc.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);

  s_flux[X1DIR].NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1+1);
  if (pmy_block->block_size.nx2 > 1)
    s_flux[X2DIR].NewAthenaArray(NSCALARS, ncells3, ncells2+1, ncells1);
  if (pmy_block->block_size.nx3 > 1)
    s_flux[X3DIR].NewAthenaArray(NSCALARS, ncells3+1, ncells2, ncells1);

  // KGF: change to RAII
  psbval  = new CellCenteredBoundaryVariable(pmy_block, &s, s_flux);
  psbval->bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(psbval);

  pmb->pbval->bvars_main_int.push_back(psbval);

  // Allocate memory for scratch arrays
  dxw_.NewAthenaArray(ncells1);
  wl_.NewAthenaArray((NWAVE), ncells1);
  wr_.NewAthenaArray((NWAVE), ncells1);
  wlb_.NewAthenaArray((NWAVE), ncells1);
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
  wl3d_.NewAthenaArray((NWAVE), ncells3, ncells2, ncells1);
  wr3d_.NewAthenaArray((NWAVE), ncells3, ncells2, ncells1);
  scr1_nkji_.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  laplacian_l_fc_.NewAthenaArray(ncells1);
  laplacian_r_fc_.NewAthenaArray(ncells1);
}

// destructor

PassiveScalars::~PassiveScalars() {
  s.DeleteAthenaArray();

  // fourth-orderintegrator
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

  // remove this class's BoundaryVariable from vector of pointers in BoundaryValues
  // KGF: assumes that this destructor is called BEFORE ~BoundaryValues())
  //  pmy_block->pbval->bvars.erase(pmy_block->pbval->bvars.begin() + phbval->bvar_index);
  // need to update the rest of the vector entries "bvar_index"!
  delete psbval;
}
