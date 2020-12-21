//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_hydro.cpp
//  \brief functions that perform flux correction for hydrodynamic variables

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // std::memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../globals.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../orbital_advection/orbital_advection.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "../bvals_cc.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn int HydroBoundaryVariable::LoadFluxBoundaryBufferSameLevel(Real *buf,
//!                                                  const NeighborBlock& nb)
//! \brief Set surface hydro flux buffers for sending to a block on the same level

int HydroBoundaryVariable::LoadFluxBoundaryBufferSameLevel(Real *buf,
                                                       const NeighborBlock& nb) {
  MeshBlock *pmb=pmy_block_;
  Real qomL = pbval_->qomL_;
  int p = 0;
  if (pbval_->shearing_box == 1 && nb.shear
      && (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1)) {
    int i;
    int sign;
    if (nb.fid == BoundaryFace::inner_x1) {
      i = pmb->is;
      sign = -1;
    } else {
      i = pmb->ie + 1;
      sign =  1;
    }
    if(pmb->porb->orbital_advection_defined) {
      for (int nn=nl_; nn<=nu_; nn++) {
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int j=pmb->js; j<=pmb->je; j++) {
            buf[p++] = x1flux(nn,k,j,i);
          }
        }
      }
    } else {
      for (int nn=nl_; nn<=nu_; nn++) {
        if(nn==IVY) {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++) {
              buf[p++] = x1flux(nn,k,j,i)+sign*qomL*x1flux(IDN,k,j,i);
            }
          }
        } else if (nn==IEN) {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++) {
              buf[p++] = x1flux(nn,k,j,i)
                           +0.5*SQR(qomL)*x1flux(IDN,k,j,i)
                           +sign*qomL*x1flux(IVY,k,j,i);
            }
          }
        } else {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++) {
              buf[p++] = x1flux(nn,k,j,i);
            }
          }
        }
      }
    }
  }
  return p;
}
