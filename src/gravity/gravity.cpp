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
#include <utility>
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
    output_source(false),
    gbvar(pmb, &phi, nullptr, empty_flux),
    accumulated_src_(pmb->ncells3, pmb->ncells2, pmb->ncells1) {
  if (four_pi_G == 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Gravity::Gravity" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  output_source = pin->GetOrAddBoolean("gravity", "output_source", false);

  // using Gravity as an example of: containing full object members instead of pointer
  // memebers, construting BoundaryVariaable composite obj (no default ctor) in Gravity
  // ctor initializer list, avoiding dynamically-managed memory and the need for a
  // user-provided dtor.

  // Enroll CellCenteredBoundaryVariable object
  gbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&gbvar);
}

void Gravity::EnrollSource(AthenaArray<Real> &arr, int idx) {
  // make sure no duplication
  UnenrollSource(arr, idx);
  enrolled_src_.emplace_back(&arr, idx);
}

void Gravity::UnenrollSource(AthenaArray<Real> &arr, int idx) {
  for (auto it = enrolled_src_.begin(); it != enrolled_src_.end(); it++) {
    if (it->first == &arr && it->second == idx) {
      // no duplication
      enrolled_src_.erase(it);
      return;
    }
  }
}

void Gravity::ComputeSource() {
  if (enrolled_src_.size() == 1) {
    // only one array enrolled, point src to the array
    auto parr = enrolled_src_.front().first;
    int idx = enrolled_src_.front().second;
    src.InitWithShallowSlice(*parr, 4, idx, 1);
  } else {
    // zero or more than one array enrolled, sum everything to buffer and point src to it
    accumulated_src_.ZeroClear();
    int sz = accumulated_src_.GetSize();
    AthenaArray<Real> arr;    // shallow copy of each enrolled array
    for (auto s : enrolled_src_) {
      auto parr = s.first;
      int idx = s.second;
      arr.InitWithShallowSlice(*parr, 4, idx, 1);
#pragma omp simd
      for (int i = 0; i < sz; i++)
        accumulated_src_(i) += arr(i);
    }
    src.InitWithShallowSlice(accumulated_src_, 4, 0, 1);
  }
}
