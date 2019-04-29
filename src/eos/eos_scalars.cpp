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
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "eos.hpp"


//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPassiveScalarsFloor(AthenaArray<Real> &prim, int i)
// \brief Apply species concentration floor to cell-averaged passive scalars or
// reconstructed L/R cell interface states (if PPM is used, e.g.)
// along (NSCALARS x) x1 slices

void EquationOfState::ApplyPassiveScalarsFloor(AthenaArray<Real> &s, int i) {
  for (int n=0; n<NSCALARS; ++n) {
    Real& s_n  = s(n,i);
    // apply (prim) density floor
    s_n = (s_n > scalar_floor_) ?  s_n : scalar_floor_;
  }
  return;
}
