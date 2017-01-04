//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena++ headers
#include "diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::Viscosity
//  \brief Adds viscous stress as flux to the hydro flux

void HydroDiffusion::Viscosity(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &cons, AthenaArray<Real> *diflx)
{
  AthenaArray<Real> &x1flux=diflx[X1DIR];
  AthenaArray<Real> &x2flux=diflx[X2DIR];
  AthenaArray<Real> &x3flux=diflx[X3DIR];
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real denf;
  Real nuiso1 = nuiso_;
  Real nuiso2 = -(2./3.);


// step-1. calculate the flux due to viscous stress
  Divv(prim,divv_);

// step-2. calculate the flux across each face
// i-direction
  for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je; ++j){
      // compute fluxes
      FaceXdx(k,j,is,ie+1,prim,fx_);
      FaceXdy(k,j,is,ie+1,prim,fy_);
      FaceXdz(k,j,is,ie+1,prim,fz_);
      // store fluxes
      for (int i=is; i<=ie+1; ++i){
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        x1flux(IM1,k,j,i) = -denf*nuiso1*(fx_(i)+nuiso2*(divv_(k,j,i)+divv_(k,j,i-1)));
        x1flux(IM2,k,j,i) = -denf*nuiso1*fy_(i);
        x1flux(IM3,k,j,i) = -denf*nuiso1*fz_(i);
        if(NON_BAROTROPIC_EOS)
          x1flux(IEN,k,j,i) = 0.5*((prim(IM1,k,j,i-1)+prim(IM1,k,j,i))*x1flux(IM1,k,j,i) +
                                (prim(IM2,k,j,i-1)+prim(IM2,k,j,i))*x1flux(IM2,k,j,i) +
                                (prim(IM3,k,j,i-1)+prim(IM3,k,j,i))*x1flux(IM3,k,j,i));
      }
  }}

// j-direction
  for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je+1; ++j){
      // compute fluxes
      FaceYdx(k,j,is,ie,prim,fx_);
      FaceYdy(k,j,is,ie,prim,fy_);
      FaceYdz(k,j,is,ie,prim,fz_);
      // store fluxes
      for(int i=is; i<=ie; i++) {
        denf = 0.5*(prim(IDN,k,j+1,i)+prim(IDN,k,j,i));
        x2flux(IM1,k,j,i) = -denf*nuiso1*fx_(i);
        x2flux(IM2,k,j,i) = -denf*nuiso1*
                         (fy_(i)+nuiso2*(divv_(k,j+1,i)+divv_(k,j,i)));
        x2flux(IM3,k,j,i) = -denf*nuiso1*fz_(i);
        if (NON_BAROTROPIC_EOS)
          x2flux(IEN,k,j,i) = 0.5*((prim(IM1,k,j,i)+prim(IM1,k,j+1,i))*x2flux(IM1,k,j,i) +
                                (prim(IM2,k,j,i)+prim(IM2,k,j+1,i))*x2flux(IM2,k,j,i) +
                                (prim(IM3,k,j,i)+prim(IM3,k,j+1,i))*x2flux(IM3,k,j,i));
      }
  }}


// k-direction
  for (int k=ks; k<=ke+1; ++k){
    for (int j=js; j<=je; ++j){
      // compute fluxes
      FaceZdx(k,j,is,ie,prim,fx_);
      FaceZdy(k,j,is,ie,prim,fy_);
      FaceZdz(k,j,is,ie,prim,fz_);
      // store fluxes
      for(int i=is; i<=ie; i++) {
        denf = 0.5*(prim(IDN,k+1,j,i)+prim(IDN,k,j,i));
        x3flux(IM1,k,j,i) = -denf*nuiso1*fx_(i);
        x3flux(IM2,k,j,i) = -denf*nuiso1*fy_(i);
        x3flux(IM3,k,j,i) = -denf*nuiso1*(fz_(i)+nuiso2*(divv_(k+1,j,i)+divv_(k,j,i)));
        if (NON_BAROTROPIC_EOS)
          x3flux(IEN,k,j,i) = 0.5*((prim(IM1,k,j,i)+prim(IM1,k+1,j,i))*x3flux(IM1,k,j,i) +
                                  (prim(IM2,k,j,i)+prim(IM2,k+1,j,i))*x3flux(IM2,k,j,i) +
                                  (prim(IM3,k,j,i)+prim(IM3,k+1,j,i))*x3flux(IM3,k,j,i));
      }
  }}

// step-2. add the flux to the total flux solved via Rieman solver
// e.g., the flux across x1-face:
//   x1flux(IM1,k,j,i)+=sigma_11;
//   x1flux(IM2,k,j,i)+=sigma_21;
//   x1flux(IM3,k,j,i)+=sigma_31;
//   x1flux(IEN,k,j,i)+=sigma_11*v1+sigma_21*v2+sigma_31*v3
// where sigma_ij is the stress tensor of viscosity


  return;
}

//-------------------------------------------------------------------------------------
// Calculate Divv
void HydroDiffusion::Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv)
{

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  int il = is-1; int iu = ie+1;
  int jl, ju, kl, ku;
  Real area_p1, area;
  Real vel_p1, vel;

  if(pmb_->block_size.nx2 == 1) // 1D
    jl=js, ju=je, kl=ks, ku=ke;
  else if(pmb_->block_size.nx3 == 1) // 2D
    jl=js-1, ju=je+1, kl=ks, ku=ke;
  else // 3D
    jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;

  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      // calculate x1-flux divergence
      pmb_->pcoord->Face1Area(k,j,il,iu+1,x1area_);
      for (int i=il; i<=iu; ++i){
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(prim(IM1,k,j,i+1)+prim(IM1,k,j,i));
        vel     = 0.5*(prim(IM1,k,j,i)+prim(IM1,k,j,i-1));
        divv(k,j,i) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      if (pmb_->block_size.nx2 > 1) {
        pmb_->pcoord->Face2Area(k,j  ,il,iu,x2area_);
        pmb_->pcoord->Face2Area(k,j+1,il,iu,x2area_p1_);
        for (int i=il; i<=iu; ++i){
          area_p1 = x2area_p1_(i);
          area    = x2area_(i);
          vel_p1  = 0.5*(prim(IM2,k,j+1,i)+prim(IM2,k,j,i));
          vel     = 0.5*(prim(IM2,k,j,i)+prim(IM2,k,j-1,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      if (pmb_->block_size.nx3 > 1) {
        pmb_->pcoord->Face3Area(k  ,j,il,iu,x3area_);
        pmb_->pcoord->Face3Area(k+1,j,il,iu,x3area_p1_);
        for (int i=il; i<=iu; ++i){
          area_p1 = x3area_p1_(i);
          area    = x3area_(i);
          vel_p1  = 0.5*(prim(IM3,k+1,j,i)+prim(IM3,k,j,i));
          vel     = 0.5*(prim(IM3,k,j,i)+prim(IM3,k-1,j,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      pmb_->pcoord->CellVolume(k,j,il,iu,vol_);
      for (int i=il; i<=iu; ++i){
        divv(k,j,i) = divv(k,j,i)/vol_(i);
      }
    }
  }

  return;
}

// v_{x1;x1}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  for (int i=il; i<=iu; ++i){
    len(i) = 2.0*(prim(IM1,k,j,i)-prim(IM1,k,j,i-1))/pco_->dx1v(i-1);
  }
  return;
}

// v_{x2;x1}+v_{x1;x2}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
  for (int i=il; i<=iu; ++i){
    len(i) = pco_->h2f(i)*(prim(IM2,k,j,i)/pco_->h2v(i)
                    -prim(IM2,k,j,i-1)/pco_->h2v(i-1))/pco_->dx1v(i-1)
             +0.5*(prim(IM1,k,j+1,i)+prim(IM1,k,j+1,i-1)-prim(IM1,k,j-1,i)-prim(IM1,k,j-1,i-1))
              /pco_->h2f(i)/(pco_->dx2v(j-1)+pco_->dx2v(j));
  }
  return;
}

// v_{x3;x1}+v_{x1;x3}  covariant derivative at x1 interface
void HydroDiffusion::FaceXdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pco_->h31f(i)*(prim(IM3,k,j,i)/pco_->h31v(i)-prim(IM3,k,j,i-1)/pco_->h31v(i-1))/pco_->dx1v(i-1)
             +0.5*(prim(IM1,k+1,j,i)+prim(IM1,k+1,j,i-1)-prim(IM1,k-1,j,i)-prim(IM1,k-1,j,i-1))
             /pco_->h31f(i)/pco_->h32v(j)/(pco_->dx3v(k-1)+pco_->dx3v(k));
  }
  return;
}

// v_{x1;x2}+v_{x2;x1}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1)
             +pco_->h2v(i)*0.5*((prim(IM2,k,j,i+1)+prim(IM2,k,j-1,i+1))/pco_->h2v(i+1)
                         -(prim(IM2,k,j,i-1)+prim(IM2,k,j-1,i-1))/pco_->h2v(i-1))
             /(pco_->dx1v(i-1)+pco_->dx1v(i));
  }
  return;
}

// v_{x2;x2}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = 2.0*(prim(IM2,k,j,i)-prim(IM2,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1)
              +(prim(IM1,k,j,i)+prim(IM1,k,j-1,i))/pco_->h2v(i)*pco_->dh2vd1(i);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void HydroDiffusion::FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pco_->h32f(j)*(prim(IM3,k,j,i)/pco_->h32f(j)-prim(IM3,k,j-1,i)/pco_->h32f(j-1))/pco_->h2v(i)/pco_->dx2v(j-1)
             +0.5*(prim(IM2,k+1,j,i)+prim(IM2,k+1,j-1,i)-prim(IM2,k-1,j,i)-prim(IM2,k-1,j-1,i))/pco_->h31v(i)/pco_->h32f(j)
             /(pco_->dx3v(k-1)+pco_->dx3v(k));
  }
  return;
}

// v_{x1;x3}+v_{x3;x1}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k-1,j,i))/pco_->dx3v(k-1)
             +0.5*pco_->h31v(i)*((prim(IM3,k,j,i+1)+prim(IM3,k-1,j,i+1))/pco_->h31v(i+1)-(prim(IM3,k,j,i-1)+prim(IM3,k-1,j,i-1))/pco_->h31v(i-1))
             /(pco_->dx1v(i-1)+pco_->dx1v(i));
  }
  return;
}

// v_{x2;x3}+v_{x3;x2}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k-1,j,i))/pco_->h31v(i)/pco_->h32v(j)/pco_->dx3v(k-1)
             +0.5*pco_->h32v(j)*((prim(IM3,k,j+1,i)+prim(IM3,k-1,j+1,i))/pco_->h32v(j+1)-(prim(IM3,k,j-1,i)+prim(IM3,k-1,j-1,i))/pco_->h32v(j-1))
             /pco_->h2v(i)/(pco_->dx2v(j-1)+pco_->dx2v(j));
  }
  return;
}

// v_{x3;x3}  covariant derivative at x3 interface
void HydroDiffusion::FaceZdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = 2.0*(prim(IM3,k,j,i)-prim(IM3,k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j) + (prim(IM1,k,j,i)+prim(IM1,k-1,j,i))*pco_->dh31vd1(i)/pco_->h31v(i) +
      (prim(IM2,k,j,i)+prim(IM2,k-1,j,i))*pco_->dh32vd2(j)/pco_->h32v(j);
  }
  return;
}


