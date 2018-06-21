//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena++ headers
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

//---------------------------------------------------------------------------------------
// Calculate isotropic thermal conduction

void HydroDiffusion::ThermalFlux_iso(const AthenaArray<Real> &prim,
              const AthenaArray<Real> &cons, AthenaArray<Real> *cndflx) {
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real kappaf, denf, dTdx, dTdy, dTdz;

  // i-direction
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) {
      if(pmb_->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        dTdx = (prim(IPR,k,j,i)/prim(IDN,k,j,i) - prim(IPR,k,j,i-1)/
                prim(IDN,k,j,i-1))/pco_->dx1v(i-1);
        x1flux(k,j,i) -= kappaf*denf*dTdx;
      }
    }
  }

  // j-direction
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx3 == 1) // 2D
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j-1,i));
          dTdy = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k,j-1,i)/
                    prim(IDN,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1);
          x2flux(k,j,i) -= kappaf*denf*dTdy;
        }
      }
    }
  } // zero flux for 1D

  // k-direction
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) // 2D or 3D
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D
      il=is-1, iu=ie+1;
  }
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k-1,j,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k-1,j,i));
          dTdz = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k-1,j,i)/
                   prim(IDN,k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
          x3flux(k,j,i) -= kappaf*denf*dTdz;
        }
      }
    }
  } // zero flux for 1D/2D

  return;
}


//---------------------------------------------------------------------------------------
// Calculate anisotropic thermal conduction

void HydroDiffusion::ThermalFlux_aniso(const AthenaArray<Real> &p,
                 const AthenaArray<Real> &c, AthenaArray<Real> *flx) {
  return;
}




//----------------------------------------------------------------------------------------
// constant viscosity

void ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, int is, int ie, int js, int je, int ks, int ke) {
  if (phdif->kappa_iso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ISO,k,j,i) = phdif->kappa_iso;
      }
    }
  }
  if (phdif->kappa_aniso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ANI,k,j,i) = phdif->kappa_aniso;
      }
    }
  }
  return;
}
