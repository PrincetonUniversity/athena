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
// \!fn void EquationOfState::PassiveScalarConservedToPrimitive(AthenaArray<Real> &s,
//           const AthenaArray<Real> &w, const AthenaArray<Real> &r_old,
//           AthenaArray<Real> &r, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku)
// \brief Converts conserved into primitive passive scalar variables

void EquationOfState::PassiveScalarConservedToPrimitive(
    AthenaArray<Real> &s, const AthenaArray<Real> &w, const AthenaArray<Real> &r_old,
    AthenaArray<Real> &r,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  for (int n=0; n<NSCALARS; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          const Real &w_d  = w(IDN,k,j,i);
          const Real di = 1.0/w_d;

          //for (int n=0; n<NSCALARS; ++n) {
          Real& s_n  = s(n,k,j,i);
          Real& r_n  = r(n,k,j,i);
          // apply passive scalars floor to conserved variable first, then transform:
          // (multi-D fluxes may have caused it to drop below floor)
          s_n = (s_n < scalar_floor_*w_d) ?  scalar_floor_*w_d : s_n;
          r_n = s_n*di;
          // TODO(felker): continue to monitor the acceptability of this absolute 0. floor
          // (may create very large global conservation violations, e.g. the first few
          // cycles of the slotted cylinder test)

          //r_n = (r_n < scalar_floor_) ? scalar_floor_ : r_n;
          //s_n = r_n*w_d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PassiveScalarPrimitiveToConserved(const AthenaArray<Real> &r
//           const AthenaArray<Real> &w, AthenaArray<Real> &s, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PassiveScalarPrimitiveToConserved(
    const AthenaArray<Real> &r, const AthenaArray<Real> &w,
    AthenaArray<Real> &s, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  for (int n=0; n<NSCALARS; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          const Real &w_d  = w(IDN,k,j,i);
          //for (int n=0; n<NSCALARS; ++n) {
          Real& s_n  = s(n,k,j,i);
          const Real& r_n  = r(n,k,j,i);
          s_n = r_n*w_d;
        }
      }
    }
  }
  return;
}

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
