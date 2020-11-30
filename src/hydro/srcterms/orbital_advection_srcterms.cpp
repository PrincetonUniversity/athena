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
//! \file orbital_advection_srcterms.cpp
//  \brief Add source terms for orbital advection
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../hydro.hpp"

// this class header
#include "hydro_srcterms.hpp"

//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::OrbitalAdvectionSourceTerms
//                          (const Real dt, const AthenaArray<Real> *flux,
//                           const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
//  \brief Add source terms for orbital advection
void HydroSourceTerms::OrbitalAdvectionSourceTerms
                       (const Real dt, const AthenaArray<Real> *flux,
                        const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Field *pf = pmb->pfield;
  OrbitalAdvection *porb = pmb->porb;
  AthenaArray<Real> &vKc = porb->vKc;
  AthenaArray<Real> &dvKc1 = porb->dvKc1;
  AthenaArray<Real> &dvKc2 = porb->dvKc2;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) { // default orbital velocity
      // dM1/dt = 2\Omega \rho vy
      // dM2/dt = -\rho vx (2-q)\Omega
      // dE/dt  = q\Omega0 (\rho vy vx -b1 b2)
      //           + 2\rho\Omega vx (vk+q\Omega x)
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den  = prim(IDN,k,j,i);
            Real mom1 = den*prim(IVX,k,j,i);
            Real vy   = prim(IVY,k,j,i);
            cons(IM1,k,j,i) +=dt*Omega_0_*2.0*(den*vy);
            cons(IM2,k,j,i) -=dt*Omega_0_*(2.0-qshear_)*mom1;
            if (NON_BAROTROPIC_EOS) {
              Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i)+flux[X1DIR](IDN,k,j,i+1))
                            +0.5*mom1;
              Real temp = -rho_v1*vy;
              if (MAGNETIC_FIELDS_ENABLED) {
                temp += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
              }
              cons(IEN,k,j,i) += -dt*temp*qshear_*Omega_0_;
            }
          }
        }
      }
    } else { // user-defined orbital velocity
      // dM1/dt = 2\Omega \rho vy + 2\Omega \rho vc
      // dM2/dt = -\rho vx (2\Omega+dvk/dx) - \rho vz dvk/dz
      // dE/dt  = (b1 b2 -\rho vy vx) dvk/dx
      //           + (b2 b3 - \rho vy vz) dvk/dz
      //           + 2\rho\Omega vx vc
      // vc     = vk+q\Omega x
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den   = prim(IDN,k,j,i);
            Real vk    = vKc(k,i);
            Real dvk_i = dvKc1(k,i);
            Real vc    = vk+qshear_*Omega_0_*pmb->pcoord->x1v(i);
            Real mom1  = den*prim(IVX,k,j,i);
            Real vy    = prim(IVY,k,j,i);
            Real src   = 2.0*Omega_0_*vc;
            // 2D components
            cons(IM1,k,j,i) += dt*(Omega_0_*2.0*(den*vy)+src*den);
            cons(IM2,k,j,i) -= dt*(2.0*Omega_0_+dvk_i)*mom1;
            if (NON_BAROTROPIC_EOS) {
              Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                            +0.5*mom1;
              Real temp = -rho_v1*vy;
              if (MAGNETIC_FIELDS_ENABLED) {
                temp += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
              }
              cons(IEN,k,j,i) += dt*(temp*dvk_i+src*rho_v1);
            }
            // 3D components
            if (pmb->block_size.nx3 > 1) {
              Real dvk_k  = dvKc2(k,i);
              if (dvk_k != 0.0) {
                Real mom3 = den*prim(IVZ,k,j,i);
                cons(IM2,k,j,i) -= dt*dvk_k*mom3;
                if (NON_BAROTROPIC_EOS) {
                  Real rho_v3 = 0.25*(flux[X3DIR](IDN,k+1,j,i)+flux[X3DIR](IDN,k,j,i))
                                +0.5*mom3;
                  Real temp = -rho_v3*vy;
                  if (MAGNETIC_FIELDS_ENABLED) {
                    cons(IEN,k,j,i) += pf->bcc(IB2,k,j,i)*pf->bcc(IB3,k,j,i);
                  }
                  cons(IEN,k,j,i) += dt*temp*dvk_k;
                }
              }
            }
          }
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) { // default orbital velocity
      if (pmb->block_size.nx3 > 1) { //3D
        // dM1/dt = 2\rho vp vc/r
        // dM2/dt = -\rho vr (2 vc /r+(dvk/dr-vk/r))-\rho vz dvk/dz
        // dM3/dt = -\rho GM z /d^3
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        //          +(dvk/dz)(b2 b3 - \rho vp vz) -GM z /d^3 \rho vz
        // vc     = vk + r\Omega
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real den    = prim(IDN,k,j,i);
              Real rv     = pmb->pcoord->x1v(i);
              Real ri     = pmb->pcoord->coord_src1_i_(i); // 1/r
              Real zv     = pmb->pcoord->x3v(k);
              Real dv     = std::sqrt(SQR(rv)+SQR(zv));    // d = sqrt(r^2+z^2)
              Real vc     = std::sqrt(gm_/dv)*rv/dv;       // vk + r\Omega
              Real rdok_i = -1.5*vc*rv/SQR(dv);            // dvk/dr-vk/r
              Real dvk_k  = -1.5*vc*zv/SQR(dv);            // dvk/dz
              Real mom1   = den*prim(IVX,k,j,i);
              Real vp     = prim(IVY,k,j,i);
              Real mom3   = den*prim(IVZ,k,j,i);
              cons(IM1,k,j,i) += 2.0*dt*vc*ri*(den*vp);
              cons(IM2,k,j,i) -= dt*((2.0*vc*ri+rdok_i)*mom1+dvk_k*mom3);
              Real src_k = -gm_*zv/(dv*dv*dv);
              cons(IM3,k,j,i) += dt*src_k*den;
              if (NON_BAROTROPIC_EOS) {
                Real rho_v3 = 0.5*(flux[X3DIR](IDN,k+1,j,i)+flux[X3DIR](IDN,k,j,i));
                cons(IEN,k,j,i) += dt*src_k*rho_v3;

                Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                              +0.5*mom1;
                rho_v3 = 0.5*rho_v3+0.5*mom3;
                Real temp1 = -rho_v1*vp;
                Real temp2 = -rho_v3*vp;
                if (MAGNETIC_FIELDS_ENABLED) {
                  temp1 += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
                  temp2 += pf->bcc(IB2,k,j,i)*pf->bcc(IB3,k,j,i);
                }
                cons(IEN,k,j,i) += dt*temp1*rdok_i+dt*temp2*dvk_k;
              }
            }
          }
        }
      } else { // 2D
        // dM1/dt = 2\rho vp vc/r
        // dM2/dt = -\rho vr (2 vc /r+(dvk/dr-vk/r))
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        // vc     = vk+r\Omega
        int ks = pmb->ks;
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den    = prim(IDN,ks,j,i);
            Real rv     = pmb->pcoord->x1v(i);
            Real ri     = pmb->pcoord->coord_src1_i_(i);
            Real vc     = std::sqrt(gm_/rv);
            Real rdok_i = -1.5*vc/rv;
            Real mom1   = den*prim(IVX,ks,j,i);
            Real vp     = prim(IVY,ks,j,i);
            cons(IM1,ks,j,i) += 2.0*dt*vc*ri*(den*vp);
            cons(IM2,ks,j,i) -= dt*(2.0*vc*ri+rdok_i)*mom1;
            if (NON_BAROTROPIC_EOS) {
              Real rho_v1 = 0.25*(flux[X1DIR](IDN,ks,j,i+1)+flux[X1DIR](IDN,ks,j,i))
                            +0.5*mom1;
              Real temp = -rho_v1*vp;
              if (MAGNETIC_FIELDS_ENABLED) {
                temp += pf->bcc(IB1,ks,j,i)*pf->bcc(IB2,ks,j,i);
              }
              cons(IEN,ks,j,i) += dt*temp*rdok_i;
            }
          }
        }
      }
    } else { // user-defined orbital velocity
      if (pmb->block_size.nx3== 1) { // 2D
        // dM1/dt = 2\rho vp vc/r
        //          + \rho (vc^2/r - gm/r^2)
        // dM2/dt = -\rho vr (2 vc /r+(dvk/dr-vk/r))
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        //          +\rho vr (vc^2/r - gm/r^2)
        // vc     = vk+r\Omega
        int ks = pmb->ks;
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den    = prim(IDN,ks,j,i);
            Real rv     = pmb->pcoord->x1v(i);
            Real ri     = pmb->pcoord->coord_src1_i_(i);
            Real vk     = vKc(ks,i);
            Real vc     = vk + Omega_0_*rv;
            Real rdok_i = dvKc1(ks,i)-vk*ri;
            Real mom1   = den*prim(IVX,ks,j,i);
            Real vp     = prim(IVY,ks,j,i);
            Real src    = (SQR(vc)-gm_/rv)*ri; // vc^2/r-gm/r^2
            cons(IM1,ks,j,i) += dt*(2.0*vc*ri*(den*vp)+src*den);
            cons(IM2,ks,j,i) -= dt*(2.0*vc*ri+rdok_i)*mom1;
            if (NON_BAROTROPIC_EOS) {
              // This is consistent with the pointmass.
              Real rho_v1 = 0.5*(flux[X1DIR](IDN,ks,j,i+1)+flux[X1DIR](IDN,ks,j,i));
              cons(IEN,ks,j,i) += dt*src*rho_v1;

              rho_v1 = 0.5*rho_v1+0.5*mom1;
              Real temp = -rho_v1*vp;
              if (MAGNETIC_FIELDS_ENABLED) {
                temp += pf->bcc(IB1,ks,j,i)*pf->bcc(IB2,ks,j,i);
              }
              cons(IEN,ks,j,i) += dt*temp*rdok_i;
            }
          }
        }
      } else { // 3D
        std::cout << "[OribtalAdvectionSourceTerm]: not compatible to 3D with "
                  << "user-defined vK in cylindrical coordinates ." << std::endl;
        return;
      }
    }
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    // dM1/dt = 2\rho vp vc/r
    //          + \rho (vc^2/r - gm/r^2)
    // dM2/dt = 2\rho vp vc/r cot(\theta)
    //          + \rho/r vc^2 cot(\theta)
    // dM3/dt = -\rho vr (2 vc /r + rdok_i)
    //          -\rho vt cot(\theta)(2 vc /r + sdvk_j)
    // dE/dt  = rdok_i(b1b3-\rho vr vp)
    //          +sdvk_j (b2 b3 - \rho vt vp)
    //          + \rho vr (vc^2/r - gm/r^2)
    //          + \rho vt/r cot(\theta) vc^2
    // vc     = vk + r sin(\theta) \Omega
    // rdok_i = dvk/dr-vk/r
    // sdvk_j = (dvk/dtheta - vk cot(\theta))/r
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        Real cv1    = pmb->pcoord->coord_src1_j_(j);   // cot(\theta)
        Real cv3    = pmb->pcoord->coord_src3_j_(j);   // cot(\theta)
        Real sv     = std::sin(pmb->pcoord->x2v(j));   // sin(\theta)
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den    = prim(IDN,k,j,i);
          Real rv     = pmb->pcoord->x1v(i);
          Real ri     = pmb->pcoord->coord_src1_i_(i); // 1/r
          Real vk     = vKc(j,i);                      // vk
          Real vc     = vk + Omega_0_*rv*sv;
          Real rdok_i = dvKc1(j,i)-vk*ri;              // dvk/dr-vk/r
          Real sdvk_j = (dvKc2(j,i)-cv3*vk)*ri;        // (dvk/dtheta - vk cot(\theta))/r

          Real mom1   = den*prim(IVX,k,j,i);
          Real mom2   = den*prim(IVY,k,j,i);
          Real vp     = prim(IVZ,k,j,i);
          Real src    = SQR(vc);
          Real grav   = gm_/rv;
          cons(IM1,k,j,i) += dt*ri*(2.0*vc*(den*vp)+(src-grav)*den);
          cons(IM2,k,j,i) += dt*ri*cv1*(2.0*vc*(den*vp)+src*den);
          cons(IM3,k,j,i) -= dt*((2.0*vc*ri+rdok_i)*mom1+(2.0*vc*ri+sdvk_j)*mom2*cv3);
          //energy
          if (NON_BAROTROPIC_EOS) {
            // This is consistent with the pointmass.
            Real rho_v1 = 0.5*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i));
            Real rho_v2 = 0.5*(flux[X2DIR](IDN,k,j+1,i)+flux[X2DIR](IDN,k,j,i));
            cons(IEN,k,j,i) += dt*((src-grav)*rho_v1+src*cv1*rho_v2)/rv;

            rho_v1 = 0.5*rho_v1+0.5*mom1;
            rho_v2 = 0.5*rho_v2+0.5*mom2;
            Real temp1 = -rho_v1*vp;
            Real temp2 = -rho_v2*vp;
            if (MAGNETIC_FIELDS_ENABLED) {
             temp1 += pf->bcc(IB1,k,j,i)*pf->bcc(IB3,k,j,i);
             temp2 += pf->bcc(IB2,k,j,i)*pf->bcc(IB3,k,j,i);
            }
            cons(IEN,k,j,i) += dt*(temp1*rdok_i+temp2*sdvk_j);
          }
        }
      }
    }
  }
  return;
}
