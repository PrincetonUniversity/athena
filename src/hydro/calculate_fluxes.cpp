//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//  \brief Calculate hydro/MHD fluxes

// C/C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"   // reapply floors to face-centered reconstructed states
#include "../field/field.hpp"
#include "../gravity/gravity.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, int order) {
  MeshBlock *pmb=pmy_block;
  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  AthenaArray<Real> b1, b2, b3, w_x1f, w_x2f, w_x3f, e2x1, e3x1, e1x2, e3x2, e1x3, e2x3;
  if (MAGNETIC_FIELDS_ENABLED) {
    b1.InitWithShallowCopy(b.x1f);
    b2.InitWithShallowCopy(b.x2f);
    b3.InitWithShallowCopy(b.x3f);
    w_x1f.InitWithShallowCopy(pmb->pfield->wght.x1f);
    w_x2f.InitWithShallowCopy(pmb->pfield->wght.x2f);
    w_x3f.InitWithShallowCopy(pmb->pfield->wght.x3f);
    e2x1.InitWithShallowCopy(pmb->pfield->e2_x1f);
    e3x1.InitWithShallowCopy(pmb->pfield->e3_x1f);
    e1x2.InitWithShallowCopy(pmb->pfield->e1_x2f);
    e3x2.InitWithShallowCopy(pmb->pfield->e3_x2f);
    e1x3.InitWithShallowCopy(pmb->pfield->e1_x3f);
    e2x3.InitWithShallowCopy(pmb->pfield->e2_x3f);
  }

  AthenaArray<Real> wl, wr, dxw;
  wl.InitWithShallowCopy(wl_);
  wr.InitWithShallowCopy(wr_);
  dxw.InitWithShallowCopy(dxw_);

  // fourth-order hydro quantities:
  // face-centered reconstructed primitive variables, fluxes, and their Laplacians:
  AthenaArray<Real> wl_fc, wr_fc, flux_fc, laplacian_l_fc, laplacian_r_fc;
  wl_fc.InitWithShallowCopy(wl_fc_);
  wr_fc.InitWithShallowCopy(wr_fc_);
  flux_fc.InitWithShallowCopy(flux_fc_);
  laplacian_l_fc.InitWithShallowCopy(scr1_nkji_);
  laplacian_r_fc.InitWithShallowCopy(scr2_nkji_);

//----------------------------------------------------------------------------------------
// i-direction

  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }

  // reconstruct L/R states
  if (order == 1) {
    pmb->precon->DonorCellX1(pmb, kl, ku, jl, ju, is, ie+1, w, bcc, wl, wr);
  } else if (order == 2) {
    pmb->precon->PiecewiseLinearX1(pmb, kl, ku, jl, ju, is, ie+1, w, bcc, wl, wr);
  } else {
    pmb->precon->PiecewiseParabolicX1(pmb, kl, ku, jl, ju, is, ie+1, w, bcc, wl, wr);
  }

  // compute fluxes, store directly into 3D arrays
  // x1flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
  // x1flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
  RiemannSolver(kl, ku, jl, ju, is, ie+1, IVX, b1, wl, wr, x1flux, e3x1, e2x1);

  // begin fourth-order hydro:
  // Compute Laplacian of primitive Riemann states on x1 faces
  pmb->pcoord->LaplacianX1(wl, laplacian_l_fc, is, ie+1, jl, ju, kl, ku, 0, NHYDRO-1);
  pmb->pcoord->LaplacianX1(wr, laplacian_r_fc, is, ie+1, jl, ju, kl, ku, 0, NHYDRO-1);

  // TODO(kfelker): assuming uniform mesh with dx1f=dx2f=dx3f, so this should factor out
  // TODO(kfelker): also, this may need to be dx1v, since Laplacian is cell-centered
  Real h = pmb->pcoord->dx1f(is);  // pco->dx1f(i); inside loop
  Real C = (h*h)/24.0;

  // Approximate x1 face-centered primitive Riemann states
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw);
        for (int i=is; i<=ie+1; ++i) {
          wl_fc_(n,k,j,i) = wl(n,k,j,i) - C*laplacian_l_fc(n,k,j,i);
          wr_fc_(n,k,j,i) = wr(n,k,j,i) - C*laplacian_r_fc(n,k,j,i);
          // reapply primitive variable floors to face-centered L/R Riemann states
          // TODO(kfelker): only needs to be called 1x for all NHYDRO
          pmb->peos->ApplyPrimitiveFloors(wl_fc_, k, j, i);
          pmb->peos->ApplyPrimitiveFloors(wr_fc_, k, j, i);
        }
      }
    }
  }

  // Compute x1 interface fluxes from face-centered primitive variables
  // TODO(kfelker): check that e3x1,e2x1 arguments added in late 2017 work here
  RiemannSolver(kl, ku, jl, ju, is, ie+1, IVX, b1, wl_fc_, wr_fc_, flux_fc, e3x1, e2x1);

  // Compute Laplacian of second-order accurate face-averaged flux on x1 faces
  pmb->pcoord->LaplacianX1(x1flux, laplacian_l_fc, is, ie+1, jl, ju, kl, ku, 0, NHYDRO-1);

  // Correct face-averaged fluxes (Guzik eq. 10)
  for(int n=0; n<NHYDRO; n++) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // Use 1-cell width ghost buffer
        if(k>=ks && k<=ke && j>=js && j<=je) {
          // Currently using x1 vectorized face spacing. Load larger buffer? Use dx1v?
          // Originally, pmb->pcoord->dx1v(i)
          pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw);
          for(int i=is; i<=ie+1; i++) {
            x1flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_l_fc(n,k,j,i);
          }
        }
      }
    }
  }
  // end x1 fourth-order hydro

  // compute weights for GS07 CT algorithm
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw);
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x1flux(IDN,k,j,i)
                      / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
        Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
        w_x1f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
      }
    }}
  }

//----------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

    // set the loop limits
    il=is, iu=ie, kl=ks, ku=ke;
    if (MAGNETIC_FIELDS_ENABLED) {
      if (pmb->block_size.nx3 == 1) // 2D
        il=is-1, iu=ie+1, kl=ks, ku=ke;
      else // 3D
        il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
    }

    // reconstruct L/R states at j
    if (order == 1) {
      pmb->precon->DonorCellX2(pmb,kl,ku,js,je+1,il,iu,w,bcc,wl,wr);
    } else if (order == 2) {
      pmb->precon->PiecewiseLinearX2(pmb,kl,ku,js,je+1,il,iu,w,bcc,wl,wr);
    } else {
      pmb->precon->PiecewiseParabolicX2(pmb,kl,ku,js,je+1,il,iu,w,bcc,wl,wr);
    }

    // compute fluxes, store directly into 3D arrays
    // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
    // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
    RiemannSolver(kl,ku,js,je+1,il,iu,IVY,b2,wl,wr,x2flux,e1x2,e3x2);

    // compute weights for GS07 CT algorithm
    if (MAGNETIC_FIELDS_ENABLED) {
      for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
        pmb->pcoord->CenterWidth2(k,j,il,iu,dxw);
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x2flux(IDN,k,j,i)
                        / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
          Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
          w_x2f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
        }
      }}
    }
  }

//----------------------------------------------------------------------------------------
// k-direction

  if (pmb->block_size.nx3 > 1) {

    // set the loop limits
    il=is, iu=ie, jl=js, ju=je;
    if (MAGNETIC_FIELDS_ENABLED)
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;

    // reconstruct L/R states at k
    if (order == 1) {
      pmb->precon->DonorCellX3(pmb,ks,ke+1,jl,ju,il,iu,w,bcc,wl,wr);
    } else if (order == 2) {
      pmb->precon->PiecewiseLinearX3(pmb,ks,ke+1,jl,ju,il,iu,w,bcc,wl,wr);
    } else {
      pmb->precon->PiecewiseParabolicX3(pmb,ks,ke+1,jl,ju,il,iu,w,bcc,wl,wr);
    }

    // compute fluxes, store directly into 3D arrays
    // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
    // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
    RiemannSolver(ks,ke+1,jl,ju,il,iu,IVZ,b3,wl,wr,x3flux,e2x3,e1x3);

    // compute weights for GS07 CT algorithm
    if (MAGNETIC_FIELDS_ENABLED) {
      for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->CenterWidth3(k,j,il,iu,dxw);
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x3flux(IDN,k,j,i)
                        / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
          Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
          w_x3f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
        }
      }}
    }
  }


  if (SELF_GRAVITY_ENABLED) AddGravityFlux(); // add gravity flux directly

  return;
}
