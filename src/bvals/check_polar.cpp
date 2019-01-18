//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file check_polar.cpp
//  \brief check compatibilty of user-selected 'polar' and/or 'polar_wedge' boundary flags
//         with compile-time and run-time solver configurations

// C headers

// C++ headers
#include <iomanip>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::CheckPolarBoundaries(void)
//  \brief Check for any compatibility issues if polar-type boundary flags are selected.
//  Called after setting 6x boundary functions in MeshBlock's BoundaryValues() constructor

void BoundaryValues::CheckPolarBoundaries() {
  // Check that spherical-like coordinates were specified at compile time (this check
  // must be manually updated when a new coordinate system is added to coordinates/)
  if (((std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0)
       && (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") != 0)
       && (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") != 0)
       // no safety-checks or restrictions on user-defined metrics
       && (std::strcmp(COORDINATE_SYSTEM, "gr_user") != 0))) {
    std::stringstream msg;
    msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
        << "The use of 'polar' or 'polar_wedge' boundary flags is limited \n"
        << "to spherical-like coordinate systems. \n"
        << "Specified COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }
  // Check that AMR is disabled (SMR is ok)
  if (pmy_mesh_->multilevel == true) {
    if (pmy_mesh_->adaptive == true) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "The use of AMR with any 'polar' or 'polar_wedge' boundary \n"
          << "flags is currently unsupported" << std::endl;
      ATHENA_ERROR(msg);
    } else {
      // TODO(kfelker): SMR check: all blocks along a pole are at the same refinement lvl
    }
  }
  // 3D spherical-like coordinates with at least one 'polar' or 'polar_wedge' x2 boundary
  if (pmy_block_->block_size.nx3 > 1) {
    // Mesh must extend from 0.0 to 2*PI (exactly) along the azimuthal x3 dimension
    if (pmy_mesh_->mesh_size.x3max != static_cast<Real>(TWO_PI)
        || pmy_mesh_->mesh_size.x3min != static_cast<Real>(0.0)) {
      // PI, TWO_PI macro constants are specified to 16 digits after the decimal pt.
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "3D spherical-like coordinates with a x2 'polar' or 'polar_wedge' \n"
          << "boundary requires that the Mesh spans the full [0, 2*pi] azimuthal \n"
          << "range. Current x2 boundary selections are: \n"
          << "ix2_bc=" << block_bcs[INNER_X2] << "\n"
          << "ox2_bc=" << block_bcs[OUTER_X2] << "\n"
          << "range. Current x3 domain limits are: \n" << std::scientific
          << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
          << "x3min=" << pmy_mesh_->mesh_size.x3min << "\n"
          << "x3max=" << pmy_mesh_->mesh_size.x3max << std::endl;
      ATHENA_ERROR(msg);
    }
    // Also, lower and upper x3 boundaries (of Mesh) must both be periodic
    if (pmy_mesh_->mesh_bcs[INNER_X3] != PERIODIC_BNDRY
        || pmy_mesh_->mesh_bcs[OUTER_X3] != PERIODIC_BNDRY) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "3D spherical-like coordinates with a x2 'polar' or 'polar_wedge' \n"
          << "boundary requires that both x3 boundary conditions are periodic \n"
          << "Current x3 boundary selections are: \n"
          << "ix3_bc=" << block_bcs[INNER_X3] << "\n"
          << "ox3_bc=" << block_bcs[OUTER_X3] << std::endl;
      ATHENA_ERROR(msg);
    }
    // 'polar' vs. 'polar_wedge' selection in 3D: use the latter if and only if the Mesh
    // doesn't span exactly x2min=0.0,x2max=pi for lower,upper boundaries, respectively
    if ((pmy_mesh_->mesh_size.x2min == static_cast<Real>(0.0)
         && block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE)
        || (pmy_mesh_->mesh_size.x2min != static_cast<Real>(0.0)
            && block_bcs[INNER_X2] == POLAR_BNDRY)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "3D spherical-like coordinates with x2min=" << pmy_mesh_->mesh_size.x2min
          << "\nwith lower boundary condition ix2_bc=" << block_bcs[INNER_X2] << "\n"
          << "Use 'polar' boundary if and only if x2min=" << std::scientific
          << std::setprecision(std::numeric_limits<Real>::max_digits10 -1) << 0.0
          << "\nUse 'polar_wedge' if x2min>" << 0.0 << std::endl;
      ATHENA_ERROR(msg);
    }
    if ((pmy_mesh_->mesh_size.x2max == static_cast<Real>(PI)
         && block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE)
        || (pmy_mesh_->mesh_size.x2max != static_cast<Real>(PI)
            && block_bcs[OUTER_X2] == POLAR_BNDRY)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "3D spherical-like coordinates with x2max=" << pmy_mesh_->mesh_size.x2max
          << "\nwith lower boundary condition ox2_bc=" << block_bcs[OUTER_X2] << "\n"
          << "Use 'polar' boundary if and only if x2max=" << std::scientific
          << std::setprecision(std::numeric_limits<Real>::max_digits10 -1) << PI
          << "\nUse 'polar_wedge' if x2max <" << PI << std::endl;
      ATHENA_ERROR(msg);
    }
    // single MeshBlock or even number of MeshBlocks wrapping around poles in 3D
    if (block_bcs[INNER_X2] == POLAR_BNDRY || block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
      if (pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Number of MeshBlocks around the pole must be 1 or even." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    if (block_bcs[OUTER_X2] == POLAR_BNDRY || block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
      if (pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Number of MeshBlocks around the pole must be 1 or even." << std::endl;
        ATHENA_ERROR(msg);
      }
    }
  } else { // 2D or 1D:
    // 'polar' boundary flag can only be used in 3D
    if (block_bcs[INNER_X2] == POLAR_BNDRY || block_bcs[OUTER_X2] == POLAR_BNDRY) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Use 'polar_wedge', not 'polar' boundary flag for 2D spherical-like \n"
          << "coordinate boundaries. Current x2 boundary selections are: \n"
          << "ix2_bc=" << block_bcs[INNER_X2] << "\n"
          << "ox2_bc=" << block_bcs[OUTER_X2] << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  return;
}
