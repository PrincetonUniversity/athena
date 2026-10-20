//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dt_growth_test.cpp
//! \brief Problem generator for testing the growth of the time step, in particular after
//! a restart.
//!
//! A user-defined time step limits the step to problem/dt_small until
//! problem/dt_jump_time, and to problem/dt_large afterwards. If dt_large is much larger
//! than the CFL time step, the time step is then limited only by the factor of two by
//! which it may grow from one cycle to the next. The initial condition is a uniform gas
//! with a smooth Gaussian bump (problem/amp, problem/width) in density and pressure.

// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
Real dt_small, dt_large, dt_jump_time;

//========================================================================================
//! \fn Real JumpingTimeStep(MeshBlock *pmb)
//! \brief User time step: dt_small before dt_jump_time and dt_large afterwards
//========================================================================================

Real JumpingTimeStep(MeshBlock *pmb) {
  return (pmb->pmy_mesh->time < dt_jump_time) ? dt_small : dt_large;
}
}  // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Reads the parameters of the user time step and enrolls it
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  dt_small = pin->GetReal("problem", "dt_small");
  dt_large = pin->GetReal("problem", "dt_large");
  dt_jump_time = pin->GetReal("problem", "dt_jump_time");
  EnrollUserTimeStepFunction(JumpingTimeStep);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Uniform gas at rest with a Gaussian bump in density and pressure
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real amp = pin->GetReal("problem", "amp");
  const Real width = pin->GetReal("problem", "width");
  const Real gm1 = peos->GetGamma() - 1.0;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        const Real bump = amp*std::exp(-SQR(pcoord->x1v(i)/width));
        phydro->u(IDN,k,j,i) = 1.0 + bump;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = (1.0 + bump)/gm1;
      }
    }
  }
  return;
}
