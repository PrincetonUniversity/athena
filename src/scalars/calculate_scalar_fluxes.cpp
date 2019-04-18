//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//  \brief Calculate hydro/MHD fluxes

// C headers

// C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"   // reapply floors to face-centered reconstructed states
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "scalars.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void PassiveScalars::CalculateFluxes(AthenaArray<Real> &s, const int order) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = s_flux[X1DIR];
  AthenaArray<Real> &x2flux = s_flux[X2DIR];
  AthenaArray<Real> &x3flux = s_flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  AthenaArray<Real> &flux_fc = scr1_nkji_;
  AthenaArray<Real> &laplacian_all_fc = scr2_nkji_;

  //--------------------------------------------------------------------------------------
  // i-direction

  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  // TODO(felker): fix loop limits for fourth-order hydro
  //  if (MAGNETIC_FIELDS_ENABLED) {
  if (pmb->block_size.nx2 > 1) {
    if (pmb->block_size.nx3 == 1) // 2D
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    else // 3D
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  }
  //  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, w, bcc, sl_, sr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, w, bcc, sl_, sr_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, w, bcc, sl_, sr_);
      }

      pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw_);
      // TODO(felker): add x1flux idn scala

      if (order == 4) {
        for (int n=0; n<NWAVE; n++) {
          for (int i=is; i<=ie+1; i++) {
            sl3d_(n,k,j,i)=sl_(n,i);
            sr3d_(n,k,j,i)=sr_(n,i);
          }
        }
      }
    }
  }

  if (order == 4) {
    // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so this should factor out
    // TODO(felker): also, this may need to be dx1v, since Laplacian is cell-centered
    Real h = pmb->pcoord->dx1f(is);  // pco->dx1f(i); inside loop
    Real C = (h*h)/24.0;

    // construct Laplacian from x1flux
    pmb->pcoord->LaplacianX1All(x1flux, laplacian_all_fc, 0, NHYDRO-1,
                                kl, ku, jl, ju, is, ie+1);

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // Compute Laplacian of primitive Riemann states on x1 faces
        for (int n=0; n<NWAVE; ++n) {
          pmb->pcoord->LaplacianX1(sl3d_, laplacian_l_fc_, n, k, j, is, ie+1);
          pmb->pcoord->LaplacianX1(sr3d_, laplacian_r_fc_, n, k, j, is, ie+1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            sl_(n,i) = sl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
            sr_(n,i) = sr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
          }
        }
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          // pmb->peos->ApplyPrimitiveFloors(sl_, k, j, i);
          // pmb->peos->ApplyPrimitiveFloors(sr_, k, j, i);
        }

        // Compute x1 interface fluxes from face-centered primitive variables
        // TODO(felker): check that e3x1,e2x1 arguments added in late 2017 work here
        pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw_);
        // RiemannSolver(k, j, is, ie+1, IVX, b1, sl_, sr_, flux_fc, e3x1, e2x1,
        //               w_x1f, dxw_);

        // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie+1; i++)
            x1flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
        }
      }
    }
  } // end if (order == 4)
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmb->block_size.nx2 > 1) {
    // set the loop limits
    il=is-1, iu=ie+1, kl=ks, ku=ke;
    // TODO(felker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->block_size.nx3 == 1) // 2D
      kl=ks, ku=ke;
    else // 3D
      kl=ks-1, ku=ke+1;
    //    }

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, w, bcc, sl_, sr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, w, bcc, sl_, sr_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, w, bcc, sl_, sr_);
      }
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, w, bcc, slb_, sr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, w, bcc, slb_, sr_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, w, bcc, slb_, sr_);
        }

        // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
        // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
        pmb->pcoord->CenterWidth2(k,j,il,iu,dxw_);
        //RiemannSolver(k, j, il, iu, IVY, b2, sl_, sr_, x2flux, e1x2, e3x2, w_x2f, dxw_);

        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              sl3d_(n,k,j,i)=sl_(n,i);
              sr3d_(n,k,j,i)=sr_(n,i);
            }
          }
        }

        // swap the arrays for the next step
        sl_.SwapAthenaArray(slb_);
      }
    }
    if (order == 4) {
      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(felker): also, this may need to be dx2v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx2f(js);  // pco->dx2f(j); inside loop
      Real C = (h*h)/24.0;

      // construct Laplacian from x2flux
      pmb->pcoord->LaplacianX2All(x2flux, laplacian_all_fc, 0, NHYDRO-1,
                                  kl, ku, js, je+1, il, iu);

      // Approximate x2 face-centered primitive Riemann states
      for (int k=kl; k<=ku; ++k) {
        for (int j=js; j<=je+1; ++j) {
          // Compute Laplacian of primitive Riemann states on x2 faces
          for (int n=0; n<NWAVE; ++n) {
            pmb->pcoord->LaplacianX2(sl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX2(sr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              sl_(n,i) = sl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              sr_(n,i) = sr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            // pmb->peos->ApplyPrimitiveFloors(sl_, k, j, i);
            // pmb->peos->ApplyPrimitiveFloors(sr_, k, j, i);
          }

          // Compute x2 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e1x2,e3x2 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth2(k,j,il,iu,dxw_);
          //RiemannSolver(k, j, il, iu, IVY, b2, sl_, sr_, flux_fc, e1x2, e3x2,
          //            w_x2f, dxw_);

          // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++)
              x2flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
          }
        }
      }
    } // end if (order == 4)
  }

  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->block_size.nx3 > 1) {
    // set the loop limits
    // TODO(felker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED)
    il=is-1, iu=ie+1, jl=js-1, ju=je+1;

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, w, bcc, sl_, sr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, w, bcc, sl_, sr_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, w, bcc, sl_, sr_);
      }
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, w, bcc, slb_, sr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, w, bcc, slb_, sr_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, w, bcc, slb_, sr_);
        }

        // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
        // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
        pmb->pcoord->CenterWidth3(k,j,il,iu,dxw_);
        //RiemannSolver(k, j, il, iu, IVZ, b3, sl_, sr_, x3flux, e2x3, e1x3, w_x3f, dxw_);

        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              sl3d_(n,k,j,i)=sl_(n,i);
              sr3d_(n,k,j,i)=sr_(n,i);
            }
          }
        }

        // swap the arrays for the next step
        sl_.SwapAthenaArray(slb_);
      }
    }
    if (order == 4) {
      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(felker): also, this may need to be dx3v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx3f(ks);  // pco->dx3f(j); inside loop
      Real C = (h*h)/24.0;

      // construct Laplacian from x3flux
      pmb->pcoord->LaplacianX3All(x3flux, laplacian_all_fc, 0, NHYDRO-1,
                                  ks, ke+1, jl, ju, il, iu);

      // Approximate x3 face-centered primitive Riemann states
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=jl; j<=ju; ++j) {
          // Compute Laplacian of primitive Riemann states on x3 faces
          for (int n=0; n<NWAVE; ++n) {
            pmb->pcoord->LaplacianX3(sl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX3(sr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              sl_(n,i) = sl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              sr_(n,i) = sr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            // pmb->peos->ApplyPrimitiveFloors(sl_, k, j, i);
            // pmb->peos->ApplyPrimitiveFloors(sr_, k, j, i);
          }

          // Compute x3 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e2x3,e1x3 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth3(k,j,il,iu,dxw_);
          // RiemannSolver(k, j, il, iu, IVZ, b3, sl_, sr_, flux_fc, e2x3, e1x3,
          //           w_x3f, dxw_);

          // Apply Laplacian of second-order accurate face-averaged flux on x3 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++)
              x3flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
          }
        }
      }
    } // end if (order == 4)
  }
  return;
}
