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
//
//======================================================================================
//! \file viscosity.cpp
//  \brief implements functions that compute viscosity terms
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../parameter_input.hpp"

// this class header
#include "viscosity.hpp"

// viscosity constructor 

Viscosity::Viscosity(Hydro *pf, ParameterInput *pin)
{
  pmy_hydro_ = pf;
  nuiso_ = pin->GetOrAddReal("problem","nuiso",0.0);
  int ncells3 =1;

// Allocate memory for scratch vectors
  int ncells1 = pf->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pf->pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pf->pmy_block->block_size.nx3 > 1) ncells3 = pf->pmy_block->block_size.nx3 + 2*(NGHOST);

  area.NewAthenaArray(ncells1);
  area_p1.NewAthenaArray(ncells1);
  vol.NewAthenaArray(ncells1);
  visflx_.NewAthenaArray((NWAVE),ncells1);
  jvisflx_j_.NewAthenaArray((NWAVE),ncells1);
  kvisflx_k_.NewAthenaArray((NWAVE),ncells2,ncells1);
  dx_.NewAthenaArray(ncells1);
  dy_.NewAthenaArray(ncells1);
  dz_.NewAthenaArray(ncells1);
  divv_.NewAthenaArray(ncells3,ncells2,ncells1); 
}
 
// destructor

Viscosity::~Viscosity()
{
  area.DeleteAthenaArray();
  area_p1.DeleteAthenaArray();
  vol.DeleteAthenaArray();
  visflx_.DeleteAthenaArray();
  jvisflx_j_.DeleteAthenaArray();
  kvisflx_k_.DeleteAthenaArray();
  dx_.DeleteAthenaArray();
  dy_.DeleteAthenaArray();
  dz_.DeleteAthenaArray();
  divv_.DeleteAthenaArray();
}

Real Viscosity::nuiso1(int n,int k, int j, int i)
{
   return (nuiso_);
}

// --------------------------------------------------------------------------
//! \fn Viscosity::nuiso2()
// \brief dilatational coefficient of viscosity or the second coefficient of viscosity

Real Viscosity::cnuiso2(int k, int j, int i)
{
   return (-2./3.);
}

void Viscosity::ViscosityTerms(const Real dt,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real s11=1., s22=1., s33=1., s12=1., s13=1., s23=1.;
  Real denf, nu;

  if(nuiso_ == 0.0) return;

// precalculate divv
  pmb->pcoord->Divv(prim, divv_);

  for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je; ++j){

      pmb->pcoord->CellVolume(k,j,is,ie,vol);

// i-direction
      pmb->pcoord->Face1Area(k,j,is,ie+1,area);

      pmb->pcoord->FaceXdx(k, j, is, ie+1, prim, dx_);
      pmb->pcoord->FaceXdy(k, j, is, ie+1, prim, dy_);
      pmb->pcoord->FaceXdz(k, j, is, ie+1, prim, dz_);

#pragma simd
      for (int i=is; i<=ie+1; ++i){
        nu = nuiso1(IM1, k,j,i);
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        visflx_(IM1,i) = denf*nu*(2.*dx_(i)+0.5*cnuiso2(k,j,i)*(divv_(i)+divv_(i-1)))*s11;
        visflx_(IM2,i) = denf*nu*dy_(i)*s12;
        visflx_(IM3,i) = denf*nu*dz_(i)*s13;
	if (NON_BAROTROPIC_EOS) 
          visflx_(IEN,i) = 0.5*((prim(IM1,k,j,i-1)+prim(IM1,k,j,i))*visflx_(IM1,i) +
				(prim(IM2,k,j,i-1)+prim(IM2,k,j,i))*visflx_(IM2,i) +
				(prim(IM3,k,j,i-1)+prim(IM3,k,j,i))*visflx_(IM3,i)); 
      }

      // update cons quantities
      for (int n=1; n<NHYDRO; ++n){
#pragma simd
        for (int i=is; i<=ie; ++i){
          cons(n,k,j,i) += dt*(area(i+1)*visflx_(n,i+1) - area(i)*visflx_(n,i))/vol(i);
        }
      }

      // add coordinate (geometric) source terms
      pmb->pcoord->VisSrcTermsX1(k,j,dt,visflx_,prim,cons);


// j-direction
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,area_p1);

	// For the first time, also need to calculate jvisflx_j_
        if (j==js) {
          pmb->pcoord->FaceYdx(k, j, is, ie, prim, dx_);
          pmb->pcoord->FaceYdy(k, j, is, ie, prim, dy_);
          pmb->pcoord->FaceYdz(k, j, is, ie, prim, dz_);

          for (int i=is; i<=ie; ++i){
	    nu = nuiso1(IM2,k,j,i);
	    denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j-1,i));
            jvisflx_j_(IM1,i) = denf*nu*dx_(i)*s12;
            jvisflx_j_(IM2,i) = denf*nu*
                                (2.*dy_(i)+0.5*cnuiso2(k,j,i)*(divv_(k,j,i)+divv_(k,j-1,i)))*s22;
            jvisflx_j_(IM3,i) = denf*nu*dz_(i)*s23;
	    if (NON_BAROTROPIC_EOS) 
              jvisflx_j_(IEN,i) = 0.5*((prim(IM1,k,j-1,i)+prim(IM1,k,j,i))*jvisflx_j_(IM1,i) +
                                       (prim(IM2,k,j-1,i)+prim(IM2,k,j,i))*jvisflx_j_(IM2,i) +
                                       (prim(IM3,k,j-1,i)+prim(IM3,k,j,i))*jvisflx_j_(IM3,i));
          }
        }
	// calculate visflx_
        pmb->pcoord->FaceYdx(k, j+1, is, ie, prim, dx_);
        pmb->pcoord->FaceYdy(k, j+1, is, ie, prim, dy_);
        pmb->pcoord->FaceYdz(k, j+1, is, ie, prim, dz_);

        for (int i=is; i<=ie; ++i){
          nu = nuiso1(IM2,k,j+1,i);
          denf = 0.5*(prim(IDN,k,j+1,i)+prim(IDN,k,j,i));
          visflx_(IM1,i) = denf*nu*dx_(i)*s12;
          visflx_(IM2,i) = denf*nu*
                           (2.*dy_(i)+0.5*cnuiso2(k,j+1,i)*(divv_(k,j+1,i)+divv_(k,j,i)))*s22;
          visflx_(IM3,i) = denf*nu*dz_(i)*s23;
	  if (NON_BAROTROPIC_EOS)
            visflx_(IEN,i) = 0.5*((prim(IM1,k,j,i)+prim(IM1,k,j+1,i))*visflx_(IM1,i) +
                                  (prim(IM2,k,j,i)+prim(IM2,k,j+1,i))*visflx_(IM2,i) +
                                  (prim(IM3,k,j,i)+prim(IM3,k,j+1,i))*visflx_(IM3,i));
        }
    
	// update cons quantities
        for (int n=1; n<NHYDRO; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            cons(n,k,j,i) += dt*(area_p1(i)*visflx_(n,i) - area(i)*jvisflx_j_(n,i))/vol(i);
          }
        }

	// add coordinate (geometric) source terms
        pmb->pcoord->VisSrcTermsX2(k,j,dt,jvisflx_j_,visflx_,prim,cons);          

	// store viscous flux at j+1 to next j
        if(j<je)
          jvisflx_j_ = visflx_;
      } //j-direction

// k-direction
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,area_p1);

	// For the first time, also need to calculate kvisflx_k_
        if (k==ks){
          pmb->pcoord->FaceZdx(k, j, is, ie, prim, dx_);
          pmb->pcoord->FaceZdy(k, j, is, ie, prim, dy_);
          pmb->pcoord->FaceZdz(k, j, is, ie, prim, dz_);

          for (int i=is; i<=ie; ++i){
	    nu = nuiso1(IM3,k,j,i);
	    denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k-1,j,i));
            kvisflx_k_(IM1,j,i) = denf*nu*dx_(i)*s13;
            kvisflx_k_(IM2,j,i) = denf*nu*dy_(i)*s23;
            kvisflx_k_(IM3,j,i) = denf*nu*
                                   (2.*dz_(i)+0.5*cnuiso2(k,j,i)*(divv_(k,j,i)+divv_(k-1,j,i)))*s33; 
            if (NON_BAROTROPIC_EOS)
              kvisflx_k_(IEN,j,i) = 0.5*((prim(IM1,k-1,j,i)+prim(IM1,k,j,i))*kvisflx_k_(IM1,j,i) +
                                         (prim(IM2,k-1,j,i)+prim(IM2,k,j,i))*kvisflx_k_(IM2,j,i) +
                                         (prim(IM3,k-1,j,i)+prim(IM3,k,j,i))*kvisflx_k_(IM3,j,i));
          }
        }
	// calculate visflx_
        pmb->pcoord->FaceZdx(k+1, j, is, ie, prim, dx_);
        pmb->pcoord->FaceZdy(k+1, j, is, ie, prim, dy_);
        pmb->pcoord->FaceZdz(k+1, j, is, ie, prim, dz_);

        for (int i=is; i<=ie; ++i){
	  nu = nuiso1(IM3,k+1,j,i);
          denf = 0.5*(prim(IDN,k+1,j,i)+prim(IDN,k,j,i));
          visflx_(IM1,i) = denf*nu*dx_(i)*s13;
          visflx_(IM2,i) = denf*nu*dy_(i)*s23;
          visflx_(IM3,i) = denf*nu*
                           (2.*dz_(i)+0.5*cnuiso2(k+1,j,i)*(divv_(k+1,j,i)+divv_(k,j,i)))*s33;                       
	  if (NON_BAROTROPIC_EOS)
            visflx_(IEN,i) = 0.5*((prim(IM1,k,j,i)+prim(IM1,k+1,j,i))*visflx_(IM1,i) +
                                  (prim(IM2,k,j,i)+prim(IM2,k+1,j,i))*visflx_(IM2,i) +
                                  (prim(IM3,k,j,i)+prim(IM3,k+1,j,i))*visflx_(IM3,i));
        }
    
	// update cons quantities
        for (int n=1; n<NHYDRO; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            cons(n,k,j,i) += dt*(area_p1(i)*visflx_(n,i) - area(i)*kvisflx_k_(n,j,i))/vol(i);
          }
        }

	// add coordinate (geometric) source terms
        pmb->pcoord->VisSrcTermsX3(k,j,dt,kvisflx_k_,visflx_,prim,cons);

	// store viscous flux at k+1 to next 
        if(k<ke){
          for (int n=1; n<NHYDRO; ++n){
#pragma simd
            for (int i=is; i<=ie; ++i){
              kvisflx_k_(n,j,i) = visflx_(n,i);
            }
          }
        }
      } //k-direction
    } //j-loop
  } //k-loop
  return;
}

Real Viscosity::VisDt(Real len, int k, int j, int i)
{
  Real Visdt;
  if(pmy_hydro_->pmy_block->block_size.nx3>1){
    Visdt = SQR(len)/6.0/nuiso1(IDN,k,j,i);
  } else {
    if(pmy_hydro_->pmy_block->block_size.nx2>1)
      Visdt = SQR(len)/8.0/nuiso1(IDN,k,j,i);
    else
      Visdt = SQR(len)/4.0/nuiso1(IDN,k,j,i);
  }
  return Visdt;
}
