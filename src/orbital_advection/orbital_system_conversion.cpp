//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file orbital_system_conversion.cpp
//! \brief functions for orbital system conversion

// C/C++ headers
#include <cstring>    // memcpy
#include <iostream>   // cout, endl

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "orbital_advection.hpp"


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::ConvertOrbitalSystem(const AthenaArray<Real> &w0,
//!                            const AthenaArray<Real> &u0, const OrbitalTransform trans)
//! \brief convert the orbital system to the normal system

void OrbitalAdvection::ConvertOrbitalSystem(const AthenaArray<Real> &w0,
                       const AthenaArray<Real> &u0, const OrbitalTransform trans) {
  int flag = static_cast<int>(trans);
  if ((orbital_system_conversion_done&flag) > 0) {
    int il = pmb_->is-(NGHOST); int jl = pmb_->js; int kl = pmb_->ks;
    int iu = pmb_->ie+(NGHOST); int ju = pmb_->je; int ku = pmb_->ke;
    if (nc2>1) {
      jl -= NGHOST;
      ju += NGHOST;
    }
    if (nc3>1) {
      kl -= NGHOST;
      ku += NGHOST;
    }
    // prim
    if((orbital_system_conversion_done&flag)%2==1) {
      if(orbital_direction == 1) {
        for(int k=kl; k<=ku; k++) {
          for(int j=jl; j<=ju; j++) {
#pragma omp simd
            for(int i=il; i<=iu; i++) {
              w_orb(IDN,k,j,i) = w0(IDN,k,j,i);
              w_orb(IVX,k,j,i) = w0(IVX,k,j,i);
              w_orb(IVY,k,j,i) = w0(IVY,k,j,i) + vKc(k,i);
              w_orb(IVZ,k,j,i) = w0(IVZ,k,j,i);
              if (NON_BAROTROPIC_EOS)
                w_orb(IPR,k,j,i) = w0(IPR,k,j,i);
            }
          }
        }
      } else if(orbital_direction ==2) {
        for(int k=kl; k<=ku; k++) {
          for(int j=jl; j<=ju; j++) {
#pragma omp simd
            for(int i=il; i<=iu; i++) {
              w_orb(IDN,k,j,i) = w0(IDN,k,j,i);
              w_orb(IVX,k,j,i) = w0(IVX,k,j,i);
              w_orb(IVY,k,j,i) = w0(IVY,k,j,i);
              w_orb(IVZ,k,j,i) = w0(IVZ,k,j,i)+ vKc(j,i);
              if (NON_BAROTROPIC_EOS)
                w_orb(IPR,k,j,i) = w0(IPR,k,j,i);
            }
          }
        }
      }
      orbital_system_conversion_done -= 1;
    }
    // cons
    if((orbital_system_conversion_done&flag)>=2) {
      if(orbital_direction == 1) {
        for(int k=kl; k<=ku; k++) {
          for(int j=jl; j<=ju; j++) {
#pragma omp simd
            for(int i=il; i<=iu; i++) {
              Real den = u0(IDN,k,j,i);
              Real m2  = u0(IM2,k,j,i);
              Real mk  = den*vKc(k,i);
              u_orb(IDN,k,j,i) = den;
              u_orb(IM1,k,j,i) = u0(IM1,k,j,i);
              u_orb(IM2,k,j,i) = m2 + mk;
              u_orb(IM3,k,j,i) = u0(IM3,k,j,i);
              if (NON_BAROTROPIC_EOS)
                u_orb(IEN,k,j,i) = u0(IEN,k,j,i)+mk*(m2+0.5*mk)/den;
            }
          }
        }
      } else if(orbital_direction ==2) {
        for(int k=kl; k<=ku; k++) {
          for(int j=jl; j<=ju; j++) {
#pragma omp simd
            for(int i=il; i<=iu; i++) {
              Real den = u0(IDN,k,j,i);
              Real m3  = u0(IM3,k,j,i);
              Real mk  = den*vKc(j,i);
              u_orb(IDN,k,j,i) = den;
              u_orb(IM1,k,j,i) = u0(IM1,k,j,i);
              u_orb(IM2,k,j,i) = u0(IM2,k,j,i);
              u_orb(IM3,k,j,i) = m3+ mk;
              if (NON_BAROTROPIC_EOS)
                u_orb(IEN,k,j,i) = u0(IEN,k,j,i)+mk*(m3+0.5*mk)/den;
            }
          }
        }
      }
      orbital_system_conversion_done -= 2;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::ResetOrbitalSystemConversionFlag()
//! \brief orbital_system_conversion_done flag for HydroDiffusion() and Outputs

void OrbitalAdvection::ResetOrbitalSystemConversionFlag() {
  orbital_system_conversion_done = 3;
  return;
}
