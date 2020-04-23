//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file add_z4c_rhs.cpp
//  \brief adds the Z4c equation RHS to the state vector

// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void Z4c::AddZ4cRHS
//  \brief Adds RHS to weighted average of variables from
//  previous step(s) of time integrator algorithm
//
//  This function operates on the interior points of the MeshBlock

void Z4c::AddZ4cRHS(AthenaArray<Real> & rhs, Real const wght, AthenaArray<Real> &u_out) {
  for(int n = 0; n < N_Z4c; ++n) {
    ILOOP3(k,j,i) {
      u_out(n,k,j,i) += wght*(pmy_block->pmy_mesh->dt)*rhs(n,k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Z4c::WeightedAveW
//  \brief Compute weighted average of state vector U in time integrator step
//
//  This function operates on the interior points of the MeshBlock

void Z4c::WeightedAve(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
                      AthenaArray<Real> &u_in2, const Real wght[3]) {
  // consider every possible simplified form of weighted sum operator:
  // U = a*U + b*U1 + c*U2
  // if c=0, c=b=0, or c=b=a=0 (in that order) to avoid extra FMA operations

  // u_in2 may be an unallocated AthenaArray if using a 2S time integrator
  if (wght[2] != 0.0) {
    for (int n=0; n<N_Z4c; ++n) {
      ILOOP3(k,j,i) {
        u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i);
      }
    }
  }
  else if (wght[1] != 0.0) {
    for (int n=0; n<N_Z4c; ++n) {
      ILOOP3(k,j,i) {
        u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i);
      }
    }
  }
  else if (wght[0] != 0.0) {
    for (int n=0; n<N_Z4c; ++n) {
      ILOOP3(k,j,i) {
        u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i);
      }
    }
  }
  else {
    u_out.ZeroClear();
  }
  return;
}

