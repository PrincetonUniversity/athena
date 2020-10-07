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
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) { // prepared orbital velocity
      // dM1/dt = 2\Omega \rho vy
      // dM2/dt = -\rho vx (2-q)\Omega
      // dE/dt  = q\Omega0 (\rho vy vx -b1 b2)
      //           + 2\rho\Omega vx (vk+q\Omega x)
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den    = prim(IDN,k,j,i);
            Real vy     = prim(IVY,k,j,i);
            Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                          + 0.5*den*prim(IVX,k,j,i);
            cons(IM1,k,j,i) += 2.0*dt*den*Omega_0_*vy;
            cons(IM2,k,j,i) += (qshear_-2.0)*dt*Omega_0_*rho_v1;
            if (NON_BAROTROPIC_EOS) {
              Real temp = -rho_v1*vy;
              if (MAGNETIC_FIELDS_ENABLED) {
                temp += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
              }
              cons(IEN,k,j,i) += -dt*temp*qshear_*Omega_0_;
            }
          }
        }
      }
    } else { // user define orbital velocity
      // dM1/dt = 2\Omega \rho vy + 2\Omega \rho(vk+q\Omega x)
      // dM2/dt = -\rho vx (2\Omega+dvk/dx) - \rho vz dvk/dz
      // dE/dt  = (b1 b2 -\rho vy vx) dvk/dx
      //           + (b2 b3 - \rho vy vz dvk/dz)
      //           + 2\rho\Omega vx (vk+q\Omega x)
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den   = prim(IDN,k,j,i);
            Real vk    = vKc(k,i);
            Real dvk_i = dvKc1(k,i);
            Real vc    = vk+qshear_*Omega_0_*pmb->pcoord->x1v(i);
            Real vy     = prim(IVY,k,j,i);
            Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                          + 0.5*den*prim(IVX,k,j,i);
            // 2D components
            Real src_i = 2.0*dt*Omega_0_*vc;
            cons(IM1,k,j,i) += 2.0*dt*den*Omega_0_*vy+den*src_i;
            cons(IM2,k,j,i) += -dt*rho_v1
                                *(2.0*Omega_0_ + dvk_i);
            if (NON_BAROTROPIC_EOS) {
              Real temp = -rho_v1*vy;
              if (MAGNETIC_FIELDS_ENABLED) {
                temp += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
              }
              cons(IEN,k,j,i) += dt*temp*dvk_i
                                 +src_i*rho_v1;
            }
            // 3D components
            if (pmb->block_size.nx3 > 1) {
              Real rho_v3 = 0.5*(flux[X3DIR](IDN,k+1,j,i)+flux[X3DIR](IDN,k,j,i));
              Real dvk_k  = dvKc2(k,i);
              if (dvk_k != 0.0) {
                Real src_j  = -dt*rho_v3*dvk_k;
                cons(IM2,k,j,i) += src_j;
                if (NON_BAROTROPIC_EOS) {
                  cons(IEN,k,j,i) += src_j*vy;
                  if (MAGNETIC_FIELDS_ENABLED) {
                    cons(IEN,k,j,i) += dt*dvk_k*pf->bcc(IB2,k,j,i)*pf->bcc(IB3,k,j,i);
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) { // prepared orbital velocity
      if (pmb->block_size.nx3 > 1) { //3D
        // dM1/dt = 2\rho vp (vk + r\Omega)/r
        //          +\rho/r((vk+r\Omega)^2-gm r^2/d^3)
        // dM2/dt = -\rho vr (vk/r+2\Omega+dvk/dr)-\rho vz dvk/dz
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        //          +(dvk/dz)(b2 b3 - \rho vp vz)
        //          +\rho vr/r ((vk+r\Omega)^2-gm/r)
        // w/o pointmass
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real den    = prim(IDN,k,j,i);               // \rho
              Real rv     = pmb->pcoord->x1v(i);           // r
              Real ri     = pmb->pcoord->coord_src1_i_(i); // 1/r
              Real zv     = pmb->pcoord->x3v(k);           // z
              Real dv     = std::sqrt(SQR(rv)+SQR(zv));    // d
              Real vc     = std::sqrt(gm_/dv)*rv/dv;       // vk + r\Omega
              Real vk     = vc - Omega_0_*rv;              // vk
              Real rdok_i = -1.5*vc*rv/SQR(dv);            // r d(vk/r)/dr
              Real dvk_k  = -1.5*vc*zv/SQR(dv);            // dvk/dz
              Real vp     = prim(IVY,k,j,i);               // vp
              Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                            + 0.5*den*prim(IVX,k,j,i);
              Real rho_v3 = 0.25*(flux[X3DIR](IDN,k+1,j,i)+flux[X3DIR](IDN,k,j,i))
                            + 0.5*den*prim(IVZ,k,j,i);

              cons(IM1,k,j,i) += 2.0*den*dt*ri*vc*vp;
              cons(IM2,k,j,i) += -dt*
                                  ((2.0*vc*ri+rdok_i)*rho_v1
                                  +dvk_k*rho_v3);
              Real src_k = -dt*gm_*zv/(dv*dv*dv);
              cons(IM3,k,j,i) += src_k*den;
              if (NON_BAROTROPIC_EOS) {
                Real temp1 = -rho_v1*prim(IVY,k,j,i);
                Real temp2 = -rho_v3*prim(IVY,k,j,i);
                if (MAGNETIC_FIELDS_ENABLED) {
                  temp1 += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
                  temp2 += pf->bcc(IB2,k,j,i)*pf->bcc(IB3,k,j,i);
                }
                cons(IEN,k,j,i) += dt*(temp1*rdok_i+temp2*dvk_k)
                                   +src_k*rho_v3;
              }
            }
          }
        }
      } else { // 2D
        // dM1/dt = 2\rho vp (vk + r\Omega)/r
        // dM2/dt = -\rho vr (vk/r+2\Omega+dvk/dr)
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        // w/o pointmass
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real den    = prim(IDN,k,j,i);
              Real rv     = pmb->pcoord->x1v(i);
              Real ri     = pmb->pcoord->coord_src1_i_(i);
              Real vc     = std::sqrt(gm_/rv);
              Real vk     = vc - Omega_0_*rv;
              Real rdok_i = -1.5*vc/rv;
              Real vp     = prim(IVY,k,j,i);
              Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                            + 0.5*den*prim(IVX,k,j,i);

              cons(IM1,k,j,i) += 2.0*dt*den*ri*vc*vp;
              cons(IM2,k,j,i) += -(2.0*vc*ri+rdok_i)*dt*rho_v1;
              if (NON_BAROTROPIC_EOS) {
                Real temp = -rho_v1*prim(IVY,k,j,i);
                if (MAGNETIC_FIELDS_ENABLED) {
                  temp += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
                }
                cons(IEN,k,j,i) += dt*temp*rdok_i;
              }
            }
          }
        }
      }
    } else { // user define orbital velocity
      if (pmb->block_size.nx3== 1) { // 2D
        // dM1/dt = 2\rho vp (vk + r\Omega)/r
        //          + \rho/r(vk+r\Omega)^2
        // dM2/dt = -\rho vr (vk/r+2\Omega+dvk/dr)
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        //          +\rho vr/r (vk+r\Omega)^2
        // w pointmass
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real den    = prim(IDN,k,j,i);
              Real rv     = pmb->pcoord->x1v(i);
              Real ri     = pmb->pcoord->coord_src1_i_(i);
              Real vk     = vKc(k,i);
              Real vc     = vk + Omega_0_*rv;
              Real rdok_i = dvKc1(k,i)-vk*ri;
              Real vp     = prim(IVY,k,j,i);
              Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                            + 0.5*den*prim(IVX,k,j,i);

              Real src_i = dt*ri*SQR(vc);
              cons(IM1,k,j,i) += 2.0*dt*den*ri*vc*vp + den*src_i;
              cons(IM2,k,j,i) += -(2.0*vc*ri+rdok_i)*dt*rho_v1;
              if (NON_BAROTROPIC_EOS) {
                Real temp = -rho_v1*prim(IVY,k,j,i);
                if (MAGNETIC_FIELDS_ENABLED) {
                  temp += pf->bcc(IB1,k,j,i)*pf->bcc(IB2,k,j,i);
                }
                cons(IEN,k,j,i) += dt*temp*rdok_i + src_i*rho_v1;
              }
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
    // dM1/dt = 2\rho vp (vk + r sin(\theta) \Omega)/r
    //          + \rho/r((vk+r sin(\theta ) \Omega)^2 - gm/r)
    // dM2/dt = 2\rho vp (vk + r sin(\theta) \Omega)/r cot(\theta)
    //          + \rho/r((vk+r sin(\theta ) \Omega)^2 cot(\theta)
    // dM3/dt = -\rho vr ((vk+2r sin(\theta)\Omega)/r+dvk/dr)
    //          -\rho vt (cot(\theta)(vk+2r sin(\theta)\Omega)+dvk/dtheta)/r
    // dE/dt  = (dvk/dr-vk/r)(b1b3-\rho vr vp)
    //          +(dvk/dtheta - vk cot(\theta))/r (b2 b3 - \rho vt vp)
    //          + \rho vr/r ((vk+r sin(\theta) \Omega)^2 - gm/r)
    //          + \rho vt/r cot(\theta)(vk+r sin(\theta) \Omega)^2
    // w/o pointmass
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        Real cv1    = pmb->pcoord->coord_src1_j_(j);  // cot(\theta)
        Real cv3    = pmb->pcoord->coord_src3_j_(j);  // cot(\theta)
        Real sv     = std::sin(pmb->pcoord->x2v(j));  // sin(\theta)
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den    = prim(IDN,k,j,i);              // \rho
          Real rv     = pmb->pcoord->x1v(i);          // r
          Real ri     = pmb->pcoord->coord_src1_i_(i);// 1/r

          Real vk     = vKc(j,i);                     // vk
          Real vc     = vk + Omega_0_*rv*sv;          //
          Real rdok_i = dvKc1(j,i)-vk*ri;             // (dvk/dr-vk/r)
          Real sdvk_j = dvKc2(j,i)-cv3*vk;            // (dvk/dtheta - vk cot(\theta))

          Real vp     = prim(IVZ,k,j,i);
          Real rho_v1 = 0.25*(flux[X1DIR](IDN,k,j,i+1)+flux[X1DIR](IDN,k,j,i))
                        + 0.5*den*prim(IVX,k,j,i);
          Real rho_v2 = 0.25*(flux[X2DIR](IDN,k,j+1,i)+flux[X2DIR](IDN,k,j,i))
                        + 0.5*den*prim(IVY,k,j,i);

          cons(IM1,k,j,i) += dt*ri*den*(2.0*vp*vc+SQR(vc)-gm_/rv);
          cons(IM2,k,j,i) += dt*ri*cv1*den*(2.0*vp+vc)*vc;
          cons(IM3,k,j,i) += -dt*((2.0*vc*ri+rdok_i)*rho_v1
                                 +(2.0*vc*cv3+sdvk_j)*ri*rho_v2);
          //energy
          if (NON_BAROTROPIC_EOS) {
            Real temp1 = -vp*rho_v1;
            Real temp2 = -vp*rho_v2;
            if (MAGNETIC_FIELDS_ENABLED) {
             temp1 += pf->bcc(IB1,k,j,i)*pf->bcc(IB3,k,j,i);
             temp2 += pf->bcc(IB2,k,j,i)*pf->bcc(IB3,k,j,i);
            }
            cons(IEN,k,j,i) += dt*(temp1*rdok_i+
                                 ri*(temp2*sdvk_j
                                     +rho_v1*(SQR(vc)-gm_/rv)
                                     +rho_v2*cv1*SQR(vc)));
          }
        }
      }
    }
  }
  return;
}
