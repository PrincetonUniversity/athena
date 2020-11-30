//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
//======================================================================================
//! \file rotatingsystem.cpp
//  \brief Adds coriolis force and centrifugal force
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"

// this class header
#include "hydro_srcterms.hpp"

//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::RotatingSystemSourceTerms
//              (const Real dt, const AthenaArray<Real> *flux,
//               const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
//  \brief source terms for the rotating system

void HydroSourceTerms::RotatingSystemSourceTerms
                 (const Real dt, const AthenaArray<Real> *flux,
                  const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  if(std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    // dM1/dt = 2 \rho \Omega vp +\rho r \Omega^2
    // dM2/dt = -2 \rho \Omega vr
    // dE/dt  = \rho r \Omega^2 vr
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den  = prim(IDN,k,j,i);
          Real mom1 = den*prim(IVX,k,j,i);
          Real rv   = pmb->pcoord->x1v(i);
          Real src  = SQR(rv*Omega_0_)*pmb->pcoord->coord_src1_i_(i); // (rOmega)^2/r
          cons(IM1,k,j,i) += dt*(2.0*Omega_0_*(den*prim(IVY,k,j,i))+den*src);
          cons(IM2,k,j,i) -= 2.0*dt*Omega_0_*mom1;
          if(NON_BAROTROPIC_EOS) {
            // This is consistent with the pointmass.
            cons(IEN,k,j,i) +=
              dt*0.5*src*(flux[X1DIR](IDN,k,j,i)+flux[X1DIR](IDN,k,j,i+1));
          }
        }
      }
    }
  } else if(std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    // dM1/dt = 2 \rho vc vp / r
    //          +\rho (vc)^2 / r
    // dM2/dt = 2 \rho vp cot(\theta) vc / r
    //          + \rho cot(\theta) (vc)^2 /r
    // dM3/dt = -\rho vr (2 vc)/r
    //          -\rho vt (2 vc) cot(\theta) /r
    // dE/dt  = \rho vr (vc)^2/r
    //          + \rho vt (vc)^2 cot(\theta) /r
    // vc     = r sin(\theta)\Omega
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        Real cv1 = pmb->pcoord->coord_src1_j_(j); // cot(theta)
        Real cv3 = pmb->pcoord->coord_src3_j_(j); // cot(\theta)
        Real sv  = std::sin(pmb->pcoord->x2v(j)); // sin(\theta)
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den  = prim(IDN,k,j,i);
          Real mom1 = den*prim(IVX,k,j,i);
          Real mom2 = den*prim(IVY,k,j,i);
          Real mom3 = den*prim(IVZ,k,j,i);
          Real rv   = pmb->pcoord->x1v(i);
          Real vc   = rv*sv*Omega_0_;
          Real ri   = pmb->pcoord->coord_src1_i_(i); // 1/r
          Real src  = SQR(vc); // vc^2
          Real force = ri*(2.0*vc*mom3+src*den);
          cons(IM1,k,j,i) += dt*force;
          cons(IM2,k,j,i) += dt*force*cv1;
          cons(IM3,k,j,i) -= 2.0*dt*ri*vc*(mom1+cv3*mom2);

          if(NON_BAROTROPIC_EOS) {
            // This is consistent with the pointmass.
            Real rho_v1 = 0.5*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i));
            Real rho_v2 = 0.5*(flux[X2DIR](IDN,k,j+1,i)+flux[X2DIR](IDN,k,j,i));
            cons(IEN,k,j,i) += dt*src*(rho_v1+cv1*rho_v2)/rv;
          }
        }
      }
    }
  }
  return;
}
