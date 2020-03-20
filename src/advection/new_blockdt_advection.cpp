//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_blockdt_advection.cpp
//  \brief computes timestep using CFL condition on a MeshBlock

// C/C++ headers
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

#include "advection.hpp"

//! \fn Real Advection::NewBlockTimeStep(void)
//  \brief calculate the minimum timestep within a MeshBlock
Real Advection::NewBlockTimeStep(void) {
  MeshBlock * pmb = pmy_block;
  int tid = 0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  Real real_max = std::numeric_limits<Real>::max();
  Real min_dt = real_max;

  // AthenaArray<Real> dt1, dt2, dt3;
  // dt1.InitWithShallowCopy(dt1_);
  // dt2.InitWithShallowCopy(dt2_);
  // dt3.InitWithShallowCopy(dt3_);

  // Just slice dt1_, dt2_, dt3_ directly
  AthenaArray<Real> &dt1 = dt1_, &dt2 = dt2_, &dt3 = dt3_;  // (x1 slices)

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

      pmb->pcoord->CenterWidth1(k, j, is, ie, dt1);
      pmb->pcoord->CenterWidth2(k, j, is, ie, dt2);
      pmb->pcoord->CenterWidth3(k, j, is, ie, dt3);

      for (int i = is; i <= ie; ++i) {
        dt1(i) /= std::abs(cx1);
        dt2(i) /= std::abs(cx2);
        dt3(i) /= std::abs(cx3);
      }
      for (int i = is; i <= ie; ++i) {
        Real & dt_1 = dt1(i);
        min_dt = std::min(min_dt, dt_1);
      }
      if (pmb->block_size.nx2 > 1) {
        for (int i = is; i <= ie; ++i) {
          Real & dt_2 = dt2(i);
          min_dt = std::min(min_dt, dt_2);
        }
      }
      if (pmb->block_size.nx3 > 1) {
        for (int i = is; i <= ie; ++i) {
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
