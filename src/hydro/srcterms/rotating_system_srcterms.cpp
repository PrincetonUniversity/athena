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
//  \brief source terms for rotating system

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
          Real den    = prim(IDN,k,j,i);               // \rho
          Real rv     = pmb->pcoord->x1v(i);           // r
          Real vc     = rv*Omega_0_;                   //
          Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                        + 0.5*den*prim(IVX,k,j,i);

          Real src_i = dt*rv*SQR(Omega_0_);
          cons(IM1,k,j,i) += den*(dt*2.0*Omega_0_*prim(IVY,k,j,i)
                                  +src_i);
          cons(IM2,k,j,i) += -dt*2.0*Omega_0_*rho_v1;
          if(NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) += src_i*rho_v1;
        }
      }
    }
  } else if(std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    // dM1/dt = 2 \rho (r sin(\theta) \Omega) vp /r
    //          +\rho (r sin(\theta) \Omega)^2/r
    // dM2/dt = 2 \rho vp cot(\theta) r sin(\theta) \Omega /r
    //          + \rho cot(\theta) (r sint(\theta)\Omega)^2 /r
    // dM3/dt = -\rho vr (2 r sin(\theta)\Omega)/r
    //          -\rho vt (2 r sin(\theta)\Omega) cot(\theta) /r
    // dE/dt  = \rho vr (r sint(\theta)\Omega)^2/r
    //          + \rho vt (r sint(\theta)\Omega)^2 cot(\theta) /r
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den = prim(IDN,k,j,i);
          Real vp  = prim(IVZ,k,j,i);
          Real rv  = pmb->pcoord->x1v(i);
          Real ri  = pmb->pcoord->coord_src1_i_(i);// 1/r
          Real vc  = rv*std::sin(pmb->pcoord->x2v(j))*Omega_0_;

          Real cv1    = pmb->pcoord->coord_src1_j_(j);
          Real cv3    = pmb->pcoord->coord_src3_j_(j);
          Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                        + 0.5*den*prim(IVX,k,j,i);
          Real rho_v2 = 0.25*(flux[X2DIR](IDN,k,j+1,i)+flux[X2DIR](IDN,k,j,i))
                        + 0.5*den*prim(IVY,k,j,i);

          Real src_i1 = dt*2.0*vc*vp*ri;
          Real src_i2 = dt*SQR(vc)*ri;
          cons(IM1,k,j,i) += den*(src_i1+src_i2);
          cons(IM2,k,j,i) += den*(src_i1+src_i2)*cv1;
          cons(IM3,k,j,i) += -2.0*dt*vc*ri*(rho_v1+rho_v2*cv3);
          if(NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) += src_i2*(rho_v1+rho_v2*cv1);
        }
      }
    }
  }
  return;
}
