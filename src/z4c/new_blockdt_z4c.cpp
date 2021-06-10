//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_blockdt_z4c.cpp
//  \brief computes timestep using CFL condition on a MeshBlock

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

#include "z4c.hpp"

//! \fn Real Z4c::NewBlockTimeStep(void)
//  \brief calculate the minimum timestep within a MeshBlock
Real Z4c::NewBlockTimeStep(void) {
  MeshBlock * pmb = pmy_block;

  Real min_dt = (FLT_MAX);

  // Just slice dt1_, dt2_, dt3_ directly
  AthenaArray<Real> &dt1 = dt1_, &dt2 = dt2_, &dt3 = dt3_;  // (x1 slices)

  for (int k = pmb->ks; k <= pmb->kve; ++k) {
    for (int j = pmb->js; j <= pmb->jve; ++j) {

      // Note: these should be uniform for Cart. and equiv. to dxn
      pmb->pcoord->CenterWidth1(k, j, pmb->is, pmb->ive, dt1);
      pmb->pcoord->CenterWidth2(k, j, pmb->is, pmb->ive, dt2);
      pmb->pcoord->CenterWidth3(k, j, pmb->is, pmb->ive, dt3);

      for (int i = pmb->is; i <= pmb->ive; ++i) {
        Real & dt_1 = dt1(i);
        min_dt = std::min(min_dt, dt_1);
      }
      if (pmb->block_size.nx2 > 1) {
        for (int i = pmb->is; i <= pmb->ive; ++i) {
          Real & dt_2 = dt2(i);
          min_dt = std::min(min_dt, dt_2);
        }
      }
      if (pmb->block_size.nx3 > 1) {
        for (int i = pmb->is; i <= pmb->ive; ++i) {
          Real & dt_3 = dt3(i);
          min_dt = std::min(min_dt, dt_3);
        }
      }
    }
  }

  min_dt *= pmb->pmy_mesh->cfl_number;

  pmb->new_block_dt_ = min_dt;
  return min_dt;
}
