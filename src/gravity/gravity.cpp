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
#include "../bvals/cc/bvals_cc.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "gravity.hpp"
#include "mg_gravity.hpp"

//! constructor, initializes data structures and parameters

// TODO(felker): change "MeshBlock *pmb" to reference member, set in initializer list
Gravity::Gravity(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb), phi(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
    four_pi_G(pmb->pmy_mesh->four_pi_G_),
    gbvar(pmb, &phi, nullptr, empty_flux) {
  if (four_pi_G == 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Gravity::Gravity" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // using Gravity as an example of: containing full object members instead of pointer
  // memebers, construting BoundaryVariaable composite obj (no default ctor) in Gravity
  // ctor initializer list, avoiding dynamically-managed memory and the need for a
  // user-provided dtor.

  // Enroll CellCenteredBoundaryVariable object
  gbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&gbvar);
}
