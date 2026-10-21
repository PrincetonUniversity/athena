//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file source_time_test.cpp
//! \brief Problem generator for testing the time integration of source terms.
//!
//! The gas is uniform and at rest, so the flux divergence vanishes and the evolution is
//! driven only by source terms, whose time integral is known exactly. Two source terms
//! can be used:
//!   - constant acceleration (hydro/grav_acc1 etc., built in): the momentum is
//!     rho*g*t exactly, and
//!   - a user-defined energy source (enrolled if problem/src_amp is not zero) that
//!     adds src_amp*cos(src_omega*t) per unit time to the total energy of every cell, so
//!     that the energy changes by src_amp*sin(src_omega*t)/src_omega. It depends on the
//!     time at which the time integrator evaluates the source terms.

// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
Real src_amp, src_omega;

//========================================================================================
//! \fn void CosineEnergySource(MeshBlock *pmb, const Real time, const Real dt,
//!                             const AthenaArray<Real> &prim,
//!                             const AthenaArray<Real> &prim_scalar,
//!                             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
//!                             AthenaArray<Real> &cons_scalar)
//! \brief Adds dt*src_amp*cos(src_omega*time) to the total energy of every real cell
//========================================================================================

void CosineEnergySource(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &prim_scalar,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                        AthenaArray<Real> &cons_scalar) {
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        cons(IEN,k,j,i) += dt*src_amp*std::cos(src_omega*time);
      }
    }
  }
}
}  // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Enrolls the user-defined energy source if problem/src_amp is not zero
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  src_amp = pin->GetReal("problem", "src_amp");
  src_omega = pin->GetReal("problem", "src_omega");
  if (src_amp != 0.0) EnrollUserExplicitSourceFunction(CosineEnergySource);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Uniform gas at rest, with density 1 and pressure 1
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real gm1 = peos->GetGamma() - 1.0;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = 1.0/gm1;
      }
    }
  }
  return;
}
