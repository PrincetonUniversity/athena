//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief functions for orbital system output

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
//! \fn void OrbitalAdvection::SetOrbitalSystemOutput(const AthenaArray<Real> &src,
//                                                    bool cons_flag)
//  \brief calculate profiles including orbital velocity
void OrbitalAdvection::SetOrbitalSystemOutput(const AthenaArray<Real> &src) {
  if(!orbital_system_output_done) {
    int il = pmb_->is-(NGHOST); int jl = pmb_->js-(NGHOST); int kl = pmb_->ks;
    int iu = pmb_->ie+(NGHOST); int ju = pmb_->je+(NGHOST); int ku = pmb_->ke;
    if (nc3>1) {
      kl -= NGHOST;
      ku += NGHOST;
    }
    if(orbital_direction == 1) {
      for(int k=kl; k<=ku; k++) {
        for(int j=jl; j<=ju; j++) {
#pragma omp simd
          for(int i=il; i<=iu; i++) {
            w_orb(IDN,k,j,i) = src(IDN,k,j,i);
            w_orb(IVX,k,j,i) = src(IVX,k,j,i);
            w_orb(IVY,k,j,i) = src(IVY,k,j,i) + vKc(k,i);
            w_orb(IVZ,k,j,i) = src(IVZ,k,j,i);
            if (NON_BAROTROPIC_EOS)
              w_orb(IPR,k,j,i) = src(IPR,k,j,i);
          }
        }
      }
    } else if(orbital_direction ==2) {
      for(int k=kl; k<=ku; k++) {
        for(int j=jl; j<=ju; j++) {
#pragma omp simd
          for(int i=il; i<=iu; i++) {
            w_orb(IDN,k,j,i) = src(IDN,k,j,i);
            w_orb(IVX,k,j,i) = src(IVX,k,j,i);
            w_orb(IVY,k,j,i) = src(IVY,k,j,i);
            w_orb(IVZ,k,j,i) = src(IVZ,k,j,i)+ vKc(j,i);
            if (NON_BAROTROPIC_EOS)
              w_orb(IPR,k,j,i) = src(IPR,k,j,i);
          }
        }
      }
    }
    pmb_->peos->PrimitiveToConserved(w_orb, pf_->bcc, u_orb,
                                     pco_, il, iu, jl, ju, kl, ku);
    orbital_system_output_done = true;
  }
  return;
}
