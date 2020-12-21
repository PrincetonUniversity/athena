//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file rotating_system_srcterms.cpp
//! \brief Adds coriolis force and centrifugal force
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
//!             (const Real dt, const AthenaArray<Real> *flux,
//!              const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
//! \brief source terms for the rotating system

void HydroSourceTerms::RotatingSystemSourceTerms
                 (const Real dt, const AthenaArray<Real> *flux,
                  const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  if(std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    // dM1/dt = 2 \rho vc vp /r +\rho vc^2/r
    // dM2/dt = -2 \rho vc vr /r
    // dE/dt  = \rho vc^2 vr /r
    // vc     = r \Omega
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den  = prim(IDN,k,j,i);
          Real mom1 = den*prim(IVX,k,j,i);
          Real ri   = pmb->pcoord->coord_src1_i_(i);
          Real rv   = pmb->pcoord->x1v(i);
          Real vc   = rv*Omega_0_;
          Real src  = SQR(vc); // (rOmega)^2
          Real flux_c = 0.5*(flux[X1DIR](IDN,k,j,i)+flux[X1DIR](IDN,k,j,i+1));
          cons(IM1,k,j,i) += dt*ri*(2.0*vc*(den*prim(IVY,k,j,i))+den*src);
          cons(IM2,k,j,i) -= dt*ri*vc*(mom1 + flux_c);
          if(NON_BAROTROPIC_EOS) {
            // This is consistent with the pointmass.
            cons(IEN,k,j,i) += dt*ri*src*flux_c;
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
          Real rv   = pmb->pcoord->x1v(i);
          Real ri   = pmb->pcoord->coord_src1_i_(i); // 1/r
          Real vc   = rv*sv*Omega_0_;
          Real src  = SQR(vc); // vc^2
          Real force = den*ri*(2.0*vc*prim(IVZ,k,j,i)+src);
          Real flux_xc = 0.5*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i));
          Real flux_yc = 0.5*(flux[X2DIR](IDN,k,j+1,i)+flux[X2DIR](IDN,k,j,i));
          cons(IM1,k,j,i) += dt*force;
          cons(IM2,k,j,i) += dt*force*cv1;
          cons(IM3,k,j,i) -= dt*ri*vc*(den*prim(IVX,k,j,i)+flux_xc
                                       +cv3*(den*prim(IVY,k,j,i)+flux_yc));
          if(NON_BAROTROPIC_EOS) {
            // This is consistent with the pointmass.
            cons(IEN,k,j,i) += dt*src*(flux_xc+cv1*flux_yc)/rv;
          }
        }
      }
    }
  }
  return;
}
