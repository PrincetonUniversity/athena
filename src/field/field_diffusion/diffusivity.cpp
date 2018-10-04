//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//=======================================================================================
//  \brief functions to calculate diffusion fluxes


// C++ headers
#include <iostream>   // cout, endl

// Athena headers
#include "../../athena.hpp"          // Real
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../../coordinates/coordinates.hpp" // Coordinates
#include "../../hydro/hydro.hpp"     // Fluid
#include "../field.hpp"              // Field
#include "field_diffusion.hpp"
#include "../../mesh/mesh.hpp"       // MeshBlock
#include "../../parameter_input.hpp" // ParameterInput

//======================================================================================
//! \file diffusivity.cpp
//  \brief implements functions that are related to non-ideal MHD diffusivities
//  Include the following functions:
//
//    ConstDiffusivity          - magnetic diffusivity with constant proportional
//                                 coefficients (default choice)
//    CalcCurrent               - calculate current density
//
//    OhmicEMF                  - calculate EMF from Ohmic resistivity
//
//    AmbipolarEMF              - calculate EMF from ambipolar diffusion
//
//    PoyntingFlux              - calculate the Poynting flux due to dissipation
//=======================================================================================

//--------------------------------------------------------------------------------------
// Magnetic diffusivity from constant coefficients

void ConstDiffusivity(FieldDiffusion *pfdif, MeshBlock *pmb, const AthenaArray<Real> &w,
     const AthenaArray<Real> &bmag, const int is, const int ie, const int js,
     const int je, const int ks, const int ke) {
  if (pfdif->eta_ohm > 0.0) { // Ohmic resistivity is turned on
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
#pragma omp simd
        for(int i=is; i<=ie; i++)
          pfdif->etaB(I_O, k,j,i) = pfdif->eta_ohm;
      }
    }
  }

  if (pfdif->eta_hall != 0.0) { // Hall diffusivity is turned on
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
#pragma omp simd
        for(int i=is; i<=ie; i++)
          pfdif->etaB(I_H, k,j,i) = pfdif->eta_hall*bmag(k,j,i)/w(IDN,k,j,i);
      }
    }
  }

  if (pfdif->eta_ad > 0.0) { // ambipolar diffusivity is turned on
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
#pragma omp simd
        for(int i=is; i<=ie; i++)
          pfdif->etaB(I_A, k,j,i) = pfdif->eta_ad*SQR(bmag(k,j,i));
      }
    }
  }

  return;
}


//-------------------------------------------------------------------------------------
// Calculate current density

void FieldDiffusion::CalcCurrent(FaceField &b) {
  MeshBlock *pmb = pmy_block;
  Coordinates *pco = pmb->pcoord;

  int ext2=0,ext3=0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  if (pmb->block_size.nx2 > 1) ext2 = 1;
  if (pmb->block_size.nx3 > 1) ext3 = 1;

  AthenaArray<Real> J1,J2,J3,b1i,b2i,b3i;
  J1.InitWithShallowCopy(jedge_.x1e);
  J2.InitWithShallowCopy(jedge_.x2e);
  J3.InitWithShallowCopy(jedge_.x3e);
  b1i.InitWithShallowCopy(b.x1f);
  b2i.InitWithShallowCopy(b.x2f);
  b3i.InitWithShallowCopy(b.x3f);

  AthenaArray<Real> area,len,len_m1;
  area.InitWithShallowCopy(face_area_);
  len.InitWithShallowCopy(edge_length_);
  len_m1.InitWithShallowCopy(edge_length_m1_);

  for (int k=ks-ext3; k<=ke+2*ext3; ++k) {
    for (int j=js-2*ext2; j<=je+2*ext2; ++j) {
      pco->VolCenterFace2Area(k-ext3,j,is-2,ie+1,area);
      pco->VolCenter3Length(k-ext3,j,is-2,ie+2,len);
#pragma omp simd
      for (int i=is-1; i<=ie+2; ++i) {
        J2(k,j,i) = -(len(i)*b3i(k,j,i) - len(i-1)*b3i(k,j,i-1))/area(i-1);
      }
    }
  }

  for (int k=ks-2*ext3; k<=ke+2*ext3; ++k) {
    for (int j=js-ext2; j<=je+2*ext2; ++j) {
      pco->VolCenterFace3Area(k,j-ext2,is-2,ie+1,area);
      pco->VolCenter2Length(k,j-ext2,is-2,ie+2,len);
#pragma omp simd
      for (int i=is-1; i<=ie+2; ++i) {
        J3(k,j,i) = (len(i)*b2i(k,j,i) - len(i-1)*b2i(k,j,i-1))/area(i-1);
      }
    }
  }

  if (pmb->block_size.nx2 > 1) {
    for (int k=ks-ext3; k<=ke+2*ext3; ++k) {
      for (int j=js-1; j<=je+2; ++j) {
        pco->VolCenterFace1Area(k-ext3,j-1,is-2,ie+2,area);
        pco->VolCenter3Length(k-ext3,j-1,is-2,ie+2,len_m1);
        pco->VolCenter3Length(k-ext3,j  ,is-2,ie+2,len);
#pragma omp simd
        for (int i=is-2; i<=ie+2; ++i) {
          J1(k,j,i) = (len(i)*b3i(k,j,i) - len_m1(i)*b3i(k,j-1,i))/area(i);
        }
      }
    }

    for (int k=ks-2*ext3; k<=ke+2*ext3; ++k) {
      for (int j=js-1; j<=je+2; ++j) {
        pco->VolCenterFace3Area(k,j-1,is-2,ie+1,area);
        pco->VolCenter1Length(k,j-1,is-2,ie+1,len_m1);
        pco->VolCenter1Length(k,j  ,is-2,ie+1,len);
#pragma omp simd
        for (int i=is-1; i<=ie+2; ++i) {
          J3(k,j,i) -= (len(i-1)*b1i(k,j,i) - len_m1(i-1)*b1i(k,j-1,i))/area(i-1);
        }
      }
    }
  }

  if (pmb->block_size.nx3 > 1) {
    for (int k=ks-1; k<=ke+2; ++k) {
      for (int j=js-1; j<=je+2; ++j) {
        pco->VolCenterFace1Area(k-1,j-1,is-2,ie+2,area);
        pco->VolCenter2Length(k-1,j-1,is-2,ie+2,len_m1);
        pco->VolCenter2Length(k  ,j-1,is-2,ie+2,len);
#pragma omp simd
        for (int i=is-2; i<=ie+2; ++i) {
          J1(k,j,i) -= (len(i)*b2i(k,j,i) - len_m1(i)*b2i(k-1,j,i))/area(i);
        }
      }
    }

    for (int k=ks-1; k<=ke+2; ++k) {
      for (int j=js-2; j<=je+2; ++j) {
        pco->VolCenterFace2Area(k-1,j,is-2,ie+1,area);
        pco->VolCenter1Length(k-1,j,is-2,ie+1,len_m1);
        pco->VolCenter1Length(k  ,j,is-2,ie+1,len);
#pragma omp simd
        for (int i=is-1; i<=ie+2; ++i) {
          J2(k,j,i) += (len(i-1)*b1i(k,j,i) - len_m1(i-1)*b1i(k-1,j,i))/area(i-1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
// EMF from Ohmic resistivity

void FieldDiffusion::OhmicEMF(const FaceField &b, const AthenaArray<Real> &bc,
                              EdgeField &e) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1,e2,e3,J1,J2,J3;
  e1.InitWithShallowCopy(e.x1e);
  e2.InitWithShallowCopy(e.x2e);
  e3.InitWithShallowCopy(e.x3e);
  J1.InitWithShallowCopy(jedge_.x1e);
  J2.InitWithShallowCopy(jedge_.x2e);
  J3.InitWithShallowCopy(jedge_.x3e);

// 1D update:
  if (pmb->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      Real eta_O = 0.5*(etaB(I_O,ks,js,i-1)+etaB(I_O,ks,js,i));

      e2(ks  ,js  ,i) += eta_O * J2(ks,js,i);
      e2(ke+1,js  ,i)  = e2(ks,js,i);
      e3(ks  ,js  ,i) += eta_O * J3(ks,js,i);
      e3(ks  ,je+1,i)  = e3(ks,js,i);
    }
    return;
  }

// 2D update:
  if (pmb->block_size.nx3 == 1) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        Real eta_O = 0.5*(etaB(I_O,ks,j,i)+etaB(I_O,ks,j-1,i));
        e1(ks  ,j,i) += eta_O * J1(ks,j,i);
        e1(ke+1,j,i)  = e1(ks,j,i);

        eta_O = 0.5*(etaB(I_O,ks,j,i)+etaB(I_O,ks,j,i-1));
        e2(ks  ,j,i) += eta_O * J2(ks,j,i);
        e2(ke+1,j,i)  = e2(ks,j,i);

        eta_O = 0.25*(etaB(I_O,ks,j  ,i)+etaB(I_O,ks,j  ,i-1)
                     +etaB(I_O,ks,j-1,i)+etaB(I_O,ks,j-1,i-1));
        e3(ks,  j,i) += eta_O * J3(ks,j,i);
      }
    }
    return;
  }

// 3D update:
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        Real eta_O = 0.25*(etaB(I_O,k  ,j,i)+etaB(I_O,k  ,j-1,i)
                          +etaB(I_O,k-1,j,i)+etaB(I_O,k-1,j-1,i));
        e1(k,j,i) += eta_O * J1(k,j,i);

        eta_O = 0.25*(etaB(I_O,k  ,j,i)+etaB(I_O,k  ,j,i-1)
                     +etaB(I_O,k-1,j,i)+etaB(I_O,k-1,j,i-1));
        e2(k,j,i) += eta_O * J2(k,j,i);

        eta_O = 0.25*(etaB(I_O,k,j  ,i)+etaB(I_O,k,j  ,i-1)
                     +etaB(I_O,k,j-1,i)+etaB(I_O,k,j-1,i-1));
        e3(k,j,i) += eta_O * J3(k,j,i);
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
// EMF from ambipolar diffusion

void FieldDiffusion::AmbipolarEMF(const FaceField &b, const AthenaArray<Real> &bc,
                                  EdgeField &e) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1,e2,e3,J1,J2,J3;
  e1.InitWithShallowCopy(e.x1e);
  e2.InitWithShallowCopy(e.x2e);
  e3.InitWithShallowCopy(e.x3e);
  J1.InitWithShallowCopy(jedge_.x1e);
  J2.InitWithShallowCopy(jedge_.x2e);
  J3.InitWithShallowCopy(jedge_.x3e);

// 1D update:
  if (pmb->block_size.nx2 == 1) {
    for (int i=is; i<=ie+1; ++i) {
      Real eta_A = 0.5*(etaB(I_A,ks,js,i-1)+etaB(I_A,ks,js,i));

      Real intBx = b.x1f(ks,js,i);
      Real intBy = 0.5*(b.x2f(ks,js,i) + b.x2f(ks,js,i-1));
      Real intBz = 0.5*(b.x3f(ks,js,i) + b.x3f(ks,js,i-1));

      Real Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
      Real JdotB = J2(ks,js,i)*intBy + J3(ks,js,i)*intBz;

      e2(ks  ,js  ,i) += eta_A * (J2(ks,js,i) - JdotB*intBy/Bsq);
      e3(ks  ,js  ,i) += eta_A * (J3(ks,js,i) - JdotB*intBz/Bsq);

      e2(ke+1,js  ,i)  = e2(ks,js,i);
      e3(ks  ,je+1,i)  = e3(ks,js,i);
    }
    return;
  }

// 2D update:
  if (pmb->block_size.nx3 == 1) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        // emf.x
        Real eta_A = 0.5*(etaB(I_A,ks,j,i)+etaB(I_A,ks,j-1,i));

        Real intJx = J1(ks,j,i);
        Real intJy = 0.25*(J2(ks,j,  i) + J2(ks,j,  i+1)
                          +J2(ks,j-1,i) + J2(ks,j-1,i+1));
        Real intJz = 0.5 *(J3(ks,j,  i)   + J3(ks,j,i+1));

        Real intBx = 0.5*(bc(IB1,ks,j,i)+bc(IB1,ks,j-1,i));
        Real intBy = b.x2f(ks,j,i);
        Real intBz = 0.5*(b.x3f(ks,j,i)+b.x3f(ks,j-1,i));

        Real Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
        Real JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        e1(ks  ,j,i) += eta_A * (J1(ks,j,i) - JdotB*intBx/Bsq);
        e1(ke+1,j,i)  = e1(ks,j,i);

        // emf.y
        eta_A = 0.5*(etaB(I_A,ks,j,i)+etaB(I_A,ks,j,i-1));

        intJx = 0.25*(J1(ks,j,i  ) + J1(ks,j+1,i  )
                     +J1(ks,j,i-1) + J1(ks,j+1,i-1));
        intJy = J2(ks,j,i);
        intJz = 0.5 *(J3(ks,j,i  ) + J3(ks,j+1,i));

        intBx = b.x1f(ks,j,i);
        intBy = 0.5*(bc(IB2,ks,j,i)+bc(IB2,ks,j,i-1));
        intBz = 0.5*(b.x3f(  ks,j,i)+b.x3f(  ks,j,i-1));

        Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        e2(ks  ,j,i) += eta_A * (J2(ks,j,i) - JdotB*intBy/Bsq);
        e2(ke+1,j,i)  = e2(ks,j,i);

        // emf.z
        eta_A = 0.25*(etaB(I_A,ks,j  ,i)+etaB(I_A,ks,j  ,i-1)
                     +etaB(I_A,ks,j-1,i)+etaB(I_A,ks,j-1,i-1));

        intJx = 0.5*(J1(ks,j,i) + J1(ks,j,i-1));
        intJy = 0.5*(J2(ks,j,i) + J2(ks,j-1,i));
        intJz = J3(ks,j,i);

        intBx = 0.5*(b.x1f(ks,j,i) + b.x1f(ks,j-1,i));
        intBy = 0.5*(b.x2f(ks,j,i) + b.x2f(ks,j,i+1));
        intBz = 0.25*(b.x3f(ks,j  ,i) + b.x3f(ks,j  ,i-1)
                     +b.x3f(ks,j-1,i) + b.x3f(ks,j-1,i-1));

        Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        e3(ks,  j,i) += eta_A * (J3(ks,j,i) - JdotB*intBz/Bsq);
        e3(ke+1,j,i) = e3(ks,j,i);
      }
    }
    return;
  }

// 3D update:
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        // emf.x
        Real eta_A = 0.25*(etaB(I_A,k  ,j,i)+etaB(I_A,k  ,j-1,i)
                          +etaB(I_A,k-1,j,i)+etaB(I_A,k-1,j-1,i));

        Real intJx = J1(k,j,i);
        Real intJy = 0.25*(J2(k,j,  i) + J2(k,j,  i+1)
                          +J2(k,j-1,i) + J2(k,j-1,i+1));
        Real intJz = 0.25*(J3(k  ,j,i)   + J3(k  ,j,i+1)
                          +J3(k-1,j,i)   + J3(k-1,j,i+1));

        Real intBx = 0.25*(bc(IB1,k  ,j,i)+bc(IB1,k  ,j-1,i)
                          +bc(IB1,k-1,j,i)+bc(IB1,k-1,j-1,i));
        Real intBy = 0.5*(b.x2f(k,j,i) + b.x2f(k-1,j,i));
        Real intBz = 0.5*(b.x3f(k,j,i) + b.x3f(k,j-1,i));

        Real Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
        Real JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        e1(k,j,i) += eta_A * (J1(k,j,i) - JdotB*intBx/Bsq);

        // emf.y
        eta_A = 0.25*(etaB(I_A,k  ,j,i)+etaB(I_A,k  ,j,i-1)
                     +etaB(I_A,k-1,j,i)+etaB(I_A,k-1,j,i-1));

        intJx = 0.25*(J1(k,j,  i) + J1(k,j,  i-1)
                     +J1(k,j+1,i) + J1(k,j+1,i-1));
        intJy = J2(k,j,i);
        intJz = 0.25*(J3(k  ,j,i)   + J3(k  ,j+1,i)
                     +J3(k-1,j,i)   + J3(k-1,j+1,i));

        intBx = 0.5*(b.x1f(k,j,i) + b.x1f(k-1,j,i));
        intBy = 0.25*(bc(IB2,k  ,j,i)+bc(IB2,k  ,j,i-1)
                     +bc(IB2,k-1,j,i)+bc(IB2,k-1,j,i-1));
        intBz = 0.5*(b.x3f(k,j,i) + b.x3f(k,j,i-1));

        Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        e2(k,j,i) += eta_A * (J2(k,j,i) - JdotB*intBy/Bsq);

        // emf.z
        eta_A = 0.25*(etaB(I_A,k,j  ,i)+etaB(I_A,k,j  ,i-1)
                     +etaB(I_A,k,j-1,i)+etaB(I_A,k,j-1,i-1));

        intJx = 0.25*(J1(k  ,j,i) + J1(k  ,j,i-1)
                     +J1(k+1,j,i) + J1(k+1,j,i-1));
        intJy = 0.25*(J2(k  ,j,i) + J2(k  ,j-1,i)
                     +J2(k+1,j,i) + J2(k+1,j-1,i));
        intJz = J3(k,j,i);

        intBx = 0.5*(b.x1f(k,j,i) + b.x1f(k,j-1,i));
        intBy = 0.5*(b.x2f(k,j,i) + b.x2f(k,j,i-1));
        intBz = 0.25*(bc(IB3,k,j  ,i)+bc(IB3,k,j  ,i-1)
                     +bc(IB3,k,j-1,i)+bc(IB3,k,j-1,i-1));

        Bsq   = SQR(intBx)+SQR(intBy)+SQR(intBz)+TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        e3(k,j,i) += eta_A * (J3(k,j,i) - JdotB*intBz/Bsq);
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//Poynting flux from non-ideal MHD

void FieldDiffusion::PoyntingFlux(EdgeField &e, const AthenaArray<Real> &bc) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1,e2,e3,f1,f2,f3;
  e1.InitWithShallowCopy(e.x1e);
  e2.InitWithShallowCopy(e.x2e);
  e3.InitWithShallowCopy(e.x3e);
  f1.InitWithShallowCopy(pflux.x1f);
  f2.InitWithShallowCopy(pflux.x2f);
  f3.InitWithShallowCopy(pflux.x3f);

// 1D update:
  if (pmb->block_size.nx2 == 1) {
#pragma omp simd
    for (int i=is; i<=ie+1; ++i) {
      f1(ks,js,i) = -0.5*(bc(IB2,ks,js,i)+bc(IB2,ks,js,i-1))*e3(ks,js,i)
                    +0.5*(bc(IB3,ks,js,i)+bc(IB3,ks,js,i-1))*e2(ks,js,i);
    }
    return;
  }

// 2D update:
  if (pmb->block_size.nx3 == 1) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        f1(ks,j,i) = -0.25*(bc(IB2,ks,j,i)+bc(IB2,ks,j,i-1))
                          *(e3(ks,j,i) + e3(ks,j+1,i))
                      +0.5*(bc(IB3,ks,j,i)+bc(IB3,ks,j,i-1))*e2(ks,j,i);

        f2(ks,j,i) = -0.5*(bc(IB3,ks,j,i)+bc(IB3,ks,j-1,i))*e1(ks,j,i)
                    +0.25*(bc(IB1,ks,j,i)+bc(IB1,ks,j-1,i))
                         *(e3(ks,j,i) + e3(ks,j,i+1));
      }
    }
    return;
  }

// 3D update:
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        f1(k,j,i) = -0.25*(bc(IB2,k,j,i)+bc(IB2,k,j,i-1))
                         *(e3(k,j,i) + e3(k,j+1,i))
                    +0.25*(bc(IB3,k,j,i)+bc(IB3,k,j,i-1))
                         *(e2(k,j,i) + e2(k+1,j,i));

        f2(k,j,i) = -0.25*(bc(IB3,k,j,i)+bc(IB3,k,j-1,i))
                         *(e1(k,j,i) + e1(k+1,j,i))
                    +0.25*(bc(IB1,k,j,i)+bc(IB1,k,j-1,i))
                         *(e3(k,j,i) + e3(k,j,i+1));

        f3(k,j,i) = -0.25*(bc(IB1,k,j,i)+bc(IB1,k-1,j,i))
                         *(e2(k,j,i) + e2(k,j,i+1))
                    +0.25*(bc(IB2,k,j,i)+bc(IB2,k-1,j,i))
                         *(e1(k,j,i) + e1(k,j+1,i));
      }
    }
  }
  return;
}
