//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_scalars.cpp
//! \brief implements functions in EquationOfState class for passive scalars

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
#include "../scalars/scalars.hpp"
#include "eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::PassiveScalarConservedToPrimitive(AthenaArray<Real> &s,
//!           const AthenaArray<Real> &u, const AthenaArray<Real> &r_old,
//!           AthenaArray<Real> &r, Coordinates *pco,
//!           int il, int iu, int jl, int ju, int kl, int ku)
//! \brief Converts conserved into primitive passive scalar variables

void EquationOfState::PassiveScalarConservedToPrimitive(
    AthenaArray<Real> &s, const AthenaArray<Real> &u, const AthenaArray<Real> &r_old,
    AthenaArray<Real> &r,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  for (int n=0; n<NSCALARS; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          const Real &d  = u(IDN,k,j,i);

          //for (int n=0; n<NSCALARS; ++n) {
          Real& s_n  = s(n,k,j,i);
          Real& r_n  = r(n,k,j,i);
          // apply passive scalars floor to conserved variable first, then transform:
          // (multi-D fluxes may have caused it to drop below floor)
          s_n = (s_n < scalar_floor_ * d) ?  scalar_floor_ * d : s_n;
          r_n = s_n/d;
          // TODO(felker): continue to monitor the acceptability of this absolute 0. floor
          // (may create very large global conservation violations, e.g. the first few
          // cycles of the slotted cylinder test)

          //r_n = (r_n < scalar_floor_) ? scalar_floor_ : r_n;
          //s_n = r_n * d;
        }
      }
    }
  }
  return;
}

// TODO(felker): a ton of overlap with ConservedToPrimitiveCellAverage in eos.hpp.
// AthenaArray<Real> targets + 2x function calls (pointwise EOS and flooring)
void EquationOfState::PassiveScalarConservedToPrimitiveCellAverage(
    AthenaArray<Real> &s, const AthenaArray<Real> &r_old, AthenaArray<Real> &r,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  MeshBlock *pmb = pmy_block_;
  Hydro *ph = pmb->phydro;
  PassiveScalars *ps = pmb->pscalars;

  int nl = 0;
  int nu = NSCALARS - 1;
  Real h = pco->dx1f(il);
  Real C = (h*h)/24.0;

  // Fourth-order accurate approx to cell-centered Hydro conserved and primitive variables
  AthenaArray<Real> &u_cc = ph->u_cc;
  // Passive scalrs
  AthenaArray<Real> &s_cc = ps->s_cc, &r_cc = ps->r_cc;
  // Laplacians of cell-averaged conserved and 2nd order accurate primitive variables
  // (reuse Hydro scratch arrays)
  AthenaArray<Real> &laplacian_cc = ph->scr1_nkji_;  // ps->scr1_nkji_; // private

  // Compute and store Laplacian of cell-averaged conserved variables
  pco->Laplacian(s, laplacian_cc, il, iu, jl, ju, kl, ku, nl, nu);

  // Compute fourth-order approximation to cell-centered conserved variables
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // We do not actually need to store all cell-centered conserved variables,
          // but the ConservedToPrimitive() implementation operates on 4D arrays
          s_cc(n,k,j,i) = s(n,k,j,i) - C*laplacian_cc(n,k,j,i);
        }
      }
    }
  }

  // Compute Laplacian of 2nd-order approximation to cell-averaged primitive variables
  pco->Laplacian(r, laplacian_cc, il, iu, jl, ju, kl, ku, nl, nu);

  // Convert cell-centered conserved values to cell-centered primitive values
  PassiveScalarConservedToPrimitive(s_cc, u_cc, r_cc, r_cc, pco, il, iu,
                                    jl, ju, kl, ku);

  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // Compute fourth-order approximation to cell-averaged primitive variables
          r(n,k,j,i) = r_cc(n,k,j,i) + C*laplacian_cc(n,k,j,i);
          // Reapply primitive variable floors to cell-average of (prim) dimensionless
          // concentration WITHOUT correcting cell-average of (cons) scalar mass
          ApplyPassiveScalarFloors(r, n, k, j, i);
          //ApplyPassiveScalarPrimitiveConservedFloors(s, w, r, n, k, j, i);
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::PassiveScalarPrimitiveToConserved(const AthenaArray<Real> &r
//!           const AthenaArray<Real> &u, AthenaArray<Real> &s, Coordinates *pco,
//!           int il, int iu, int jl, int ju, int kl, int ku);
//! \brief Converts primitive variables into conservative variables

void EquationOfState::PassiveScalarPrimitiveToConserved(
    const AthenaArray<Real> &r, const AthenaArray<Real> &u,
    AthenaArray<Real> &s, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  for (int n=0; n<NSCALARS; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          const Real &d  = u(IDN,k,j,i);
          //for (int n=0; n<NSCALARS; ++n) {
          Real& s_n  = s(n,k,j,i);
          const Real& r_n  = r(n,k,j,i);
          s_n = r_n * d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ApplyPassiveScalarFloors(AthenaArray<Real> &prim, int n,
//!                                                     int k, int j, int i)
//! \brief Apply species concentration floor to cell-averaged passive scalars or
//! reconstructed L/R cell interface states (if PPM is used, e.g.) along:
//! (NSCALARS x) x1 slices

void EquationOfState::ApplyPassiveScalarFloors(AthenaArray<Real> &r, int n, int k, int j,
                                               int i) {
  // TODO(felker): process user-input "hydro/sfloor" in each EquationOfState ctor
  // 8x .cpp files + more in general/. Is there a better way to avoid code duplication?

  // currently, assumes same floor is applied to all NSCALARS species
  // TODO(felker): generalize this to allow separate floors per species
  for (int n=0; n<NSCALARS; ++n) {
    Real& r_n  = r(n,i);
    // apply (prim) dimensionless concentration floor WITHOUT adjusting passive scalar
    // mass (conserved), unlike in floor in standard EOS
    r_n = (r_n > scalar_floor_) ?  r_n : scalar_floor_;
  }
  return;
}

// currently unused. previously, only called in above 4th order routine:
// PassiveScalarConservedToPrimitiveCellAverage()

void EquationOfState::ApplyPassiveScalarPrimitiveConservedFloors(
    AthenaArray<Real> &s, const AthenaArray<Real> &w, AthenaArray<Real> &r,
    int n, int k, int j, int i) {
  const Real &w_d  = w(IDN,k,j,i);
  const Real di = 1.0/w_d;
  Real& s_n  = s(n,k,j,i);
  Real& r_n  = r(n,k,j,i);

  s_n = (s_n < scalar_floor_*w_d) ?  scalar_floor_*w_d : s_n;

  // this next line, when applied indiscriminately, erases the accuracy gains performed in
  // the 4th order stencils, since <r> != <s>*<1/di>, in general
  r_n = s_n*di;
  // however, if r_n is riding the variable floor, it probably should be applied so that
  // s_n = rho*r_n is consistent (more concerned with conservation than order of accuracy
  // when quantities are floored)
  return;
}
