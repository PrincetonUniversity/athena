//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_scalars.cpp
//  \brief implements functions in EquationOfState class for passive scalars

// C headers

// C++ headers
#include <cmath>   // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "eos.hpp"

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPassiveScalarFloors(AthenaArray<Real> &prim, int i)
// \brief Apply species concentration floor to cell-averaged passive scalars or
// reconstructed L/R cell interface states (if PPM is used, e.g.) along:
// (NSCALARS x) x1 slices

void EquationOfState::ApplyPassiveScalarFloors(AthenaArray<Real> &r, int i) {
  // TODO(felker): process user-input "hydro/sfloor" in each EquationOfState ctor
  // 8x .cpp files + more in general/. Is there a better way to avoid code duplication?

  // currently, assumes same floor is applied to all NSCALARS species
  // TODO(felker): generalize this to allow separate floors per species
  for (int n=0; n<NSCALARS; ++n) {
    Real& r_n  = r(n,i);
    // apply (prim) density floor
    r_n = (r_n > scalar_floor_) ?  r_n : scalar_floor_;
  }
  return;
}
