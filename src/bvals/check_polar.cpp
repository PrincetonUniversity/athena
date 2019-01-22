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
    }
  }

  // Mesh must extend from 0.0 to PI (exactly) along the polar x2 dimension
  if ((block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE || block_bcs[INNER_X2] == POLAR_BNDRY)
      && pmy_mesh_->mesh_size.x2min != static_cast<Real>(0.0)) {
    // PI, TWO_PI macro constants are specified to 16 digits after the decimal pt.
    std::stringstream msg;
    msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
        << "2D or 3D spherical-like coordinates with a lower x2 boundary flag \n"
        << "'polar' or 'polar_wedge' requires that the Mesh extends to exactly \n"
        << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
        << "x2min=" << std::scientific << 0.0 << "\n"
        << "Current x2 boundary selections are: \n"
        << "ix2_bc=" << GetBoundaryString(block_bcs[INNER_X2]) << "\n"
        << "ox2_bc=" << GetBoundaryString(block_bcs[OUTER_X2]) << "\n"
        << "Current x2 domain limits are: \n"
        << "x2min=" << pmy_mesh_->mesh_size.x2min << "\n"
        << "x2max=" << pmy_mesh_->mesh_size.x2max << std::endl;
    ATHENA_ERROR(msg);
  }
  if ((block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE || block_bcs[OUTER_X2] == POLAR_BNDRY)
      && pmy_mesh_->mesh_size.x2max != static_cast<Real>(PI)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
        << "2D or 3D spherical-like coordinates with an upper x2 boundary flag \n"
        << "'polar' or 'polar_wedge' requires that the Mesh extends to exactly \n"
        << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
        << "x2max=" << std::scientific << PI << "\n"
        << "Current x2 boundary selections are: \n"
        << "ix2_bc=" << GetBoundaryString(block_bcs[INNER_X2]) << "\n"
        << "ox2_bc=" << GetBoundaryString(block_bcs[OUTER_X2]) << "\n"
        << "Current x2 domain limits are: \n"
        << "x2min=" << pmy_mesh_->mesh_size.x2min << "\n"
        << "x2max=" << pmy_mesh_->mesh_size.x2max << std::endl;
    ATHENA_ERROR(msg);
  }

  // 3D spherical-like coordinates with at least one 'polar' or 'polar_wedge' x2 boundary
  if (pmy_block_->block_size.nx3 > 1) {
    // 'polar' vs. 'polar_wedge' selection in 3D: use the latter if and only if the Mesh
    // doesn't span exactly x3min=0.0,x3max=2*pi for lower,upper boundaries, respectively
    if ((pmy_mesh_->mesh_size.x3min != static_cast<Real>(0.0)
         || pmy_mesh_->mesh_size.x3max != static_cast<Real>(TWO_PI))
        && (block_bcs[INNER_X2] == POLAR_BNDRY || block_bcs[OUTER_X2] == POLAR_BNDRY)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "3D spherical-like coordinates with current x3 domain limits:\n"
          << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
          << "x3min=" << std::scientific << pmy_mesh_->mesh_size.x3min << "\n"
          << "x3max=" << pmy_mesh_->mesh_size.x3max << "\n"
          << "Use 'polar' boundary flag(s) if and only if: \n"
          << "x3min=" << 0.0 << "\nAND\n"
          << "x3max=" << TWO_PI << "\n"
          << "Otherwise, use 'polar_wedge' flag(s)" << std::endl;
      ATHENA_ERROR(msg);
    }
    if (block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE
        || block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
      if ((pmy_mesh_->mesh_size.x3min == static_cast<Real>(0.0)
           && pmy_mesh_->mesh_size.x3max == static_cast<Real>(TWO_PI))) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "3D spherical-like coordinates with current x3 domain limits:\n"
            << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
            << "x3min=" << std::scientific << pmy_mesh_->mesh_size.x3min << "\n"
            << "x3max=" << pmy_mesh_->mesh_size.x3max << "\n"
            << "Use 'polar_wedge' boundary flag(s) if and only if: \n"
            << "x3min >" << 0.0 << "\nOR\n"
            << "x3max <" << TWO_PI << "\n"
            << "Otherwise, use 'polar' flag(s)" << std::endl;
        ATHENA_ERROR(msg);
      }
      // in 3D, 'polar_wedge' only makes sense if each MeshBlock spans dx3=pi/n, n=int
      Real block_dx3 = pmy_block_->block_size.x3max - pmy_block_->block_size.x3min;
      if (std::fmod(PI, block_dx3) > std::numeric_limits<Real>::epsilon()) {
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "3D spherical-like coordinates with current x3 domain limits:\n"
            << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
            << "x3min=" << std::scientific << pmy_mesh_->mesh_size.x3min << "\n"
            << "x3max=" << pmy_mesh_->mesh_size.x3max << "\n"
            << "Use 'polar_wedge' boundary flag(s) if and only if each MeshBlock\n"
            << "spans dx3=PI/n, for some integer n. Current MeshBlock azimuthal length:\n"
            << "dx3=" << block_dx3 << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    // Azimuthal x3 boundaries (of Mesh) must both be periodic if 'polar' or 'polar_wedge'
    // flags are used even once
    if (pmy_mesh_->mesh_bcs[INNER_X3] != PERIODIC_BNDRY
        || pmy_mesh_->mesh_bcs[OUTER_X3] != PERIODIC_BNDRY) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "3D spherical-like coordinates with at least one x2 'polar' \n"
          << "boundary requires that both x3 boundary conditions are periodic \n"
          << "Current x3 boundary selections are: \n"
          << "ix3_bc=" << GetBoundaryString(block_bcs[INNER_X3]) << "\n"
          << "ox3_bc=" << GetBoundaryString(block_bcs[OUTER_X3]) << std::endl;
      ATHENA_ERROR(msg);
    }

    // even number of cells in the azimuthal direction for the root grid of the Mesh
    if (pmy_mesh_->mesh_size.nx3%2!=0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "The use of 'polar' or 'polar_wedge' boundary flags in 3D\n"
          << "requires a Mesh with an even number of cells in the\n"
          << "azimuthal direction." << std::endl;
      ATHENA_ERROR(msg);
    }
    // single MeshBlock or even number of MeshBlocks wrapping around poles in 3D
    if (pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "The use of 'polar' or 'polar_wedge' boundary flags in 3D\n"
          << "requires that the number of MeshBlocks around the pole\n"
          << "must be 1 or even." << std::endl;
      ATHENA_ERROR(msg);
    }
  } else { // 2D or 1D:
    // 'polar' boundary flag can only be used in 3D
    if (block_bcs[INNER_X2] == POLAR_BNDRY || block_bcs[OUTER_X2] == POLAR_BNDRY) {
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Use 'polar_wedge', not 'polar' boundary flag for 2D spherical-like \n"
          << "coordinate boundaries. Current x2 boundary selections are: \n"
          << "ix2_bc=" << GetBoundaryString(block_bcs[INNER_X2]) << "\n"
          << "ox2_bc=" << GetBoundaryString(block_bcs[OUTER_X2]) << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  return;
}
