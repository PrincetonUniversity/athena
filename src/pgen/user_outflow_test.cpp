//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file user_outflow_test.cpp
//! \brief Problem generator for testing user-defined boundary functions against the
//! built-in outflow boundary condition.
//!
//! With problem/user_bc=1 the zero-gradient (outflow) boundary condition in x1 is
//! implemented by user-defined boundary functions that copy the first active cell into
//! the ghost cells of the primitive array they are given. Run with
//! mesh/ix1_bc=ox1_bc=user, the result must be identical to the built-in outflow
//! condition (mesh/ix1_bc=ox1_bc=outflow). The initial condition is a uniform gas at
//! rest with a smooth Gaussian bump (problem/amp, problem/width) in density and pressure
//! that reaches the boundaries during the run.

// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
//========================================================================================
//! \fn void UserOutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//!                             FaceField &b, Real time, Real dt,
//!                             int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief Zero-gradient boundary condition on the inner x1 face for the hydro variables
//========================================================================================

void UserOutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                        FaceField &b, Real time, Real dt,
                        int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,il-i) = prim(n,k,j,il);
        }
      }
    }
  }
}

//========================================================================================
//! \fn void UserOutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//!                             FaceField &b, Real time, Real dt,
//!                             int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief Zero-gradient boundary condition on the outer x1 face for the hydro variables
//========================================================================================

void UserOutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                        FaceField &b, Real time, Real dt,
                        int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,iu+i) = prim(n,k,j,iu);
        }
      }
    }
  }
}
}  // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Enrolls the user-defined outflow boundary functions if problem/user_bc=1
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (pin->GetInteger("problem", "user_bc") == 1) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, UserOutflowInnerX1);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, UserOutflowOuterX1);
  }
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
