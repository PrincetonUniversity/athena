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
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"   // reapply floors to face-centered reconstructed states
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b, FaceField &b_fc,
                            AthenaArray<Real> &bcc, AthenaArray<Real> &bcc_center,
                            const int order) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  // b,bcc are passed as fn parameters becausse clients may want to pass different bcc1,
  // b1, b2, etc., but the remaining members of the Field class are accessed directly via
  // pointers because they are unique. NOTE: b, bcc are nullptrs if no MHD.
#if MAGNETIC_FIELDS_ENABLED
  // used only to pass to (up-to) 2x RiemannSolver() calls per dimension:
  // x1:
  AthenaArray<Real> &b1 = b.x1f, &w_x1f = pmb->pfield->wght.x1f,
                  &e3x1 = pmb->pfield->e3_x1f, &e2x1 = pmb->pfield->e2_x1f;
  // x2:
  AthenaArray<Real> &b2 = b.x2f, &w_x2f = pmb->pfield->wght.x2f,
                  &e1x2 = pmb->pfield->e1_x2f, &e3x2 = pmb->pfield->e3_x2f;
  // x3:
  AthenaArray<Real> &b3 = b.x3f, &w_x3f = pmb->pfield->wght.x3f,
                  &e1x3 = pmb->pfield->e1_x3f, &e2x3 = pmb->pfield->e2_x3f;
#endif
  AthenaArray<Real> &flux_fc = scr1_nkji_;
  AthenaArray<Real> &laplacian_all_fc = scr2_nkji_;


  // -------------- fourth-order MHD indices and variables ------
  int il_buf, iu_buf, jl_buf, ju_buf, kl_buf, ku_buf;
#if MAGNETIC_FIELDS_ENABLED
  AthenaArray<Real> &b1_fc = b_fc.x1f, &b2_fc = b_fc.x2f, &b3_fc = b_fc.x3f;

  AthenaArray<Real> flux_fc_IBY, flux_fc_IBZ;
  flux_fc_IBY.InitWithShallowSlice(flux_fc, 4, IBY, 1);
  flux_fc_IBZ.InitWithShallowSlice(flux_fc, 4, IBZ, 1);
  // Reconstruct (in x2) velocity L/R states along x1 faces to corners
  AthenaArray<Real> &v_SE = pmb->pfield->v_SE, &v_NE = pmb->pfield->v_NE,
                    &v_NW = pmb->pfield->v_NW, &v_SW = pmb->pfield->v_SW;
  // Reconstruct (in x2) the single state face-averaged b_x.x1f
  AthenaArray<Real> &bx_N = pmb->pfield->bx_N, &bx_S = pmb->pfield->bx_S;
  // Reconstruct (in x1) velocity L/R states along x2 faces to corners
  AthenaArray<Real> &vl_temp = pmb->pfield->vl_temp_, &vr_temp = pmb->pfield->vr_temp_;
  // Reconstruct (in x1) the single state face-averaged b_y.x2f
  AthenaArray<Real> &by_E = pmb->pfield->by_E, &by_W = pmb->pfield->by_W;

  // 3D UCT states:
  //  AthenaArray<Real> alpha_plus_x3_, alpha_minus_x3_;
  AthenaArray<Real> &bz_R1 = pmb->pfield->bz_R1, &bz_L1 = pmb->pfield->bz_L1,
                    &bz_R2 = pmb->pfield->bz_R2, &bz_L2 = pmb->pfield->bz_L2,
                    &by_R3 = pmb->pfield->by_R3, &by_L3 = pmb->pfield->by_L3,
                    &bx_R3 = pmb->pfield->bx_R3, &bx_L3 = pmb->pfield->bx_L3;
  AthenaArray<Real> &v_R3R2 = pmb->pfield->v_R3R2, &v_R3L2 = pmb->pfield->v_R3L2,
                    &v_L3R2 = pmb->pfield->v_L3R2, &v_L3L2 = pmb->pfield->v_L3L2,
                    &v_R3R1 = pmb->pfield->v_R3R1, &v_R3L1 = pmb->pfield->v_R3L1,
                    &v_L3R1 = pmb->pfield->v_L3R1, &v_L3L1 = pmb->pfield->v_L3L1;

  // -------------- end fourth-order MHD indices and variables ------
#endif
  //--------------------------------------------------------------------------------------
  // i-direction

  AthenaArray<Real> &x1flux = flux[X1DIR];
  // set the loop limits

  jl = js, ju = je, kl = ks, ku = ke;
  // TODO(felker): fix loop limits for fourth-order hydro
  //  if (MAGNETIC_FIELDS_ENABLED) {
  if (pmb->block_size.nx2 > 1) {
    if (pmb->block_size.nx3 == 1) // 2D
      jl = js-1, ju = je+1, kl = ks, ku = ke;
    else // 3D
      jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
  }

  // TODO(felker): reconcile 4th order MHD loop limit settings:
  //   jl=js, ju=je, kl=ks, ku=ke;

  // if (order != 4) {
  //   if (MAGNETIC_FIELDS_ENABLED) {
  //     if (pmb->block_size.nx2 > 1) {
  //       if (pmb->block_size.nx3 == 1) // 2D
  //         jl=js-1, ju=je+1, kl=ks, ku=ke;
  //       else // 3D
  //         jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  //     }
  //   }
  // } else { // fourth-order corrections
  //   if (MAGNETIC_FIELDS_ENABLED) {
  //     if (pmb->block_size.nx2 > 1) {
  //       if(pmb->block_size.nx3 == 1) { // 2D
  //         jl=js-3, ju=je+3, kl=ks, ku=ke;
  //       } else { // 3D
  //         jl=js-3, ju=je+3, kl=ks-3, ku=ke+3;
  //       }
  //     }
  //   } else { // fourth-order hydro has same reqs as second-order MHD
  //     if (pmb->block_size.nx2 > 1) {
  //       if (pmb->block_size.nx3 == 1) // 2D
  //         jl=js-1, ju=je+1, kl=ks, ku=ke;
  //       else // 3D
  //         jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  //     }
  //   }
  // } // end order==4

  // // Set transverse loop limits for quantities calculated to full accuracy
  // // TODO(kfelker): check fourth-order MHD dependence
  // jl_buf=jl, ju_buf=ju, kl_buf=kl, ku_buf=ku; // 1D defaults
  // if (order == 4) {
  //   if(pmb->block_size.nx2 > 1) {
  //     if(pmb->block_size.nx3 == 1) // 2D
  //       jl_buf+=1, ju_buf-=1;
  //     else // 3D
  //       jl_buf+=1, ju_buf-=1, kl_buf+=1, ku_buf-=1;
  //   }


  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      }

      pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
      RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, x1flux, dxw_);
#else  // MHD:
      // x1flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
      // x1flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
      RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_);
#endif

      if (order == 4) {
        for (int n=0; n<NWAVE; n++) {
          for (int i=is; i<=ie+1; i++) {
            wl3d_(n,k,j,i) = wl_(n,i);
            wr3d_(n,k,j,i) = wr_(n,i);
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
          pmb->pcoord->LaplacianX1(wl3d_, laplacian_l_fc_, n, k, j, is, ie+1);
          pmb->pcoord->LaplacianX1(wr3d_, laplacian_r_fc_, n, k, j, is, ie+1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
            wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            // KGF: temporary load of x1-sliced arrays into 3D arrays for MHD4
            if (MAGNETIC_FIELDS_ENABLED) {
              wl_fc_(n,k,j,i) = wl_(n,i);
              wr_fc_(n,k,j,i) = wr_(n,i);
            }
          }
        }
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          pmb->peos->ApplyPrimitiveFloors(wl_, k, j, i);
          pmb->peos->ApplyPrimitiveFloors(wr_, k, j, i);
        }

        // Compute x1 interface fluxes from face-centered primitive variables
        // TODO(felker): check that e3x1,e2x1 arguments added in late 2017 work here
        pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, flux_fc, dxw_);
#else  // MHD:
        RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, flux_fc, e3x1, e2x1,
                      w_x1f, dxw_);
#endif
        // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie+1; i++) {
            x1flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
            // TODO(felker): replace this loop-based deep copy with memcpy, or alternative
            if (n == IDN && NSCALARS > 0) {
              pmb->pscalars->mass_flux_fc[X1DIR](k,j,i) = flux_fc(n,k,j,i);
            }
          }
        }
      }
    }
  } // end if (order == 4)
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro and MHD

#if MAGNETIC_FIELDS_ENABLED
    //-------- begin fourth-order upwind constrained transport (UCT4x1)
    if (order == 4) {
      // Currently, 2D domain is assumed for UCT4. 1D domains work via trivial copying
      // of fluid fluxes, but cannot call PPMx2() on 1D domain
      if (pmb->block_size.nx2 > 1) {
        // Unlike standard Athena++ E_z^c upwinding, which requires loading
        // [is-1:ie+1] x [js-1:je+1]
        // face-states, the UCT e_NE, e_SE, ... quantities are centered on the corner,
        // so only the real range, including the uppermost corner are required
        // [is:ie+1] x [js,je+1]

        // Limited transverse reconstructions: call PPMx2() for vx, vy both L/R states
        // wl_{i-1/2}  is E side of interface --> discontinuous states in x2, L=N,  R=S
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je+1; ++j) {
            pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, wl3d_, v_NE, v_SE,
                                              IVX, IVX, 0);
            pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, wl3d_, v_NE, v_SE,
                                              IVY, IVY, 1);
            // wr_{i-1/2}  is W side of interface --> discontinuous states in x2, L=N, R=S
            pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, wr3d_, v_NW, v_SW,
                                              IVX, IVX, 0);
            pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, wr3d_, v_NW, v_SW,
                                              IVY, IVY, 1);
            // Limited transverse reconstructions: call PPMx2() for single-state b_x
            pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, b1, bx_N, bx_S,
                                              0, 0, 0);
          }
        }
        // Repeat calculation of x1 edge-centered wavespeeds as in HLL solver
        Real wli[NWAVE], wri[NWAVE];
        int ivx = IVX;
        int ivy = IVX + ((ivx-IVX)+1) % 3;
        int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=jl_buf; j<=ju_buf; ++j) {
            for (int i=is; i<=ie+1; ++i) {
              //--- Load L/R states into local variables
              // UCT with face-centered quantities
              wli[IDN] = wl_fc_(IDN,k,j,i);
              wli[IVX] = wl_fc_(ivx,k,j,i);
              wli[IVY] = wl_fc_(ivy,k,j,i);
              wli[IVZ] = wl_fc_(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wli[IPR] = wl_fc_(IPR,k,j,i);
              wli[IBY] = wl_fc_(IBY,k,j,i);
              wli[IBZ] = wl_fc_(IBZ,k,j,i);

              wri[IDN] = wr_fc_(IDN,k,j,i);
              wri[IVX] = wr_fc_(ivx,k,j,i);
              wri[IVY] = wr_fc_(ivy,k,j,i);
              wri[IVZ] = wr_fc_(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wri[IPR] = wr_fc_(IPR,k,j,i);
              wri[IBY] = wr_fc_(IBY,k,j,i);
              wri[IBZ] = wr_fc_(IBZ,k,j,i);
              Real bxi =  b1_fc(k,j,i);

              Real cl = pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr = pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al = std::min((wri[IVX]-cr),(wli[IVX] - cl));
              Real ar = std::max((wli[IVX] + cl),(wri[IVX] + cr));
              Real bp = ar > 0.0 ? ar : 0.0;
              Real bm = al < 0.0 ? al : 0.0;
              pmb->pfield->alpha_plus_x1_(k,j,i) = bp;
              pmb->pfield->alpha_minus_x1_(k,j,i) = bm;
            }
          }
        }
        // Compute states for 3D UCT
        if (pmb->block_size.nx3 > 1) {
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je; ++j) {
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wl3d_, v_L3L1, v_R3L1,
                                                IVX, IVX, 0);
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wl3d_, v_L3L1, v_R3L1,
                                                IVZ, IVZ, 2);
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wr3d_, v_L3R1, v_R3R1,
                                                IVX, IVX, 0);
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wr3d_, v_L3R1, v_R3R1,
                                                IVZ, IVZ, 2);
              // Limited transverse reconstructions: call PPMx3() for single-state b_x
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, b1, bx_L3, bx_R3,
                                                0, 0, 0);
            }
          }
        } // end UCT if 3D
      } // end if 2D or 3D
    }

//     else { // end if (order == 4) UCT4x1
//       // compute weights for GS07 CT algorithm
//       for (int k=kl_buf; k<=ku_buf; ++k) {
//         for (int j=jl_buf; j<=ju_buf; ++j) {
//           pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw);
// #pragma omp simd
//           for (int i=is; i<=ie+1; ++i) {
//             Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x1flux(IDN,k,j,i)
//                 / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
//             Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
//             w_x1f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
//           }
//         }
//       }
//     } // end GS07 CT algorithm
#endif  // MAGNETIC_FIELDS_ENABLED

//----------------------------------------------------------------------------------------
// j-direction
  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = flux[X2DIR];
    // set the loop limits
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    // TODO(felker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->block_size.nx3 == 1) // 2D
      kl = ks, ku = ke;
    else // 3D
      kl = ks-1, ku = ke+1;
    //    }

    // TODO(felker): reconcile 4th order MHD loop limit settings:
    // set the loop limits
    // il=is, iu=ie, kl=ks, ku=ke;

    // if (order != 4) {
    //   if (MAGNETIC_FIELDS_ENABLED) {
    //     if (pmb->block_size.nx3 == 1) // 2D
    //       il=is-1, iu=ie+1, kl=ks, ku=ke;
    //     else // 3D
    //       il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
    //   }
    // } else { // fourth-order corrections
    //   if (MAGNETIC_FIELDS_ENABLED) {
    //     if(pmb->block_size.nx3 == 1) { // 2D
    //       il=is-3, iu=ie+3, kl=ks, ku=ke;
    //     } else { // 3D
    //       il=is-3, iu=ie+3, kl=ks-3, ku=ke+3;
    //     }
    //   } else { // fourth-order hydro has same reqs as second-order MHD
    //     if (pmb->block_size.nx3 == 1) // 2D
    //       il=is-1, iu=ie+1, kl=ks, ku=ke;
    //     else // 3D
    //       il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
    //   }
    // } // end order==4

    // // Set transverse loop limits for quantities calculated to full accuracy
    // il_buf=il, iu_buf=iu, kl_buf=kl, ku_buf=ku;
    // if (order == 4) {
    //   if(pmb->block_size.nx3 == 1) // 2D
    //     il_buf+=1, iu_buf-=1;
    //   else // 3D
    //     il_buf+=1, iu_buf-=1, kl_buf+=1, ku_buf-=1;
    // }


    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      }
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, w, bcc, wlb_, wr_);
        }

        pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, il, iu, IVY, wl_, wr_, x2flux, dxw_);
#else  // MHD:
        // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
        // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
        RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, x2flux, e1x2, e3x2, w_x2f, dxw_);
#endif

       // begin x2 fourth-order hydro and MHD
    //------------------------------------------------------------------------------
        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              wl3d_(n,k,j,i) = wl_(n,i);
              wr3d_(n,k,j,i) = wr_(n,i);
            }
          }
        }
        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
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
            pmb->pcoord->LaplacianX2(wl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX2(wr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
              // KGF: temporary load of x1-sliced arrays into 3D arrays for MHD4
              if (MAGNETIC_FIELDS_ENABLED) {
                wl_fc_(n,k,j,i) = wl_(n,i);
                wr_fc_(n,k,j,i) = wr_(n,i);
              }
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            pmb->peos->ApplyPrimitiveFloors(wl_, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_, k, j, i);
          }

          // Compute x2 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e1x2,e3x2 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVY, wl_, wr_, flux_fc, dxw_);
#else  // MHD:
          RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, flux_fc, e1x2, e3x2,
                        w_x2f, dxw_);
#endif

          // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++) {
              x2flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
              if (n == IDN && NSCALARS > 0) {
                pmb->pscalars->mass_flux_fc[X2DIR](k,j,i) = flux_fc(n,k,j,i);
              }
            }
          }
        }
      }
    } // end if (order == 4)
    //------------------------------------------------------------------------------
    // end x2 fourth-order hydro and MHD

#if MAGNETIC_FIELDS_ENABLED
      //-------- begin fourth-order upwind constrained transport (UCT4x2)
      if (order == 4) {

        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je+1; ++j) {
            // Limited transverse reconstructions: call PPMx1() for vx, vy both L/R states
            // wl_{i,j-1/2} is N side of interface --> discontinuous states in x1, L=E,  R=W
            pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wl3d_,
                                              vl_temp, vr_temp, IVX, IVX, 0);
            pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wl3d_,
                                              vl_temp, vr_temp, IVY, IVY, 1);
          }
        }
        // Store temporary arrays as average of R_x[R_y[]] and R_y[R_x[]] reconstructions
        for (int n=0; n<2; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je+1; ++j) {
              for (int i=is; i<=ie+1; ++i) {
                v_NE(n,k,j,i) = 0.5*(v_NE(n,k,j,i) + vl_temp(n,k,j,i));
                v_NW(n,k,j,i) = 0.5*(v_NW(n,k,j,i) + vr_temp(n,k,j,i));
              }
            }
          }
        }

        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je+1; ++j) {
        // wr_{i,j-1/2}  is S side of interface --> discontinuous states in x1, L=E, R=W
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wr3d_,
                                             vl_temp, vr_temp, IVX, IVX, 0);
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wr3d_,
                                             vl_temp, vr_temp, IVY, IVY, 1);
          }
        }
        // Store temporary arrays as average of R_x[R_y[]] and R_y[R_x[]] reconstructions
        for (int n=0; n<2; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je+1; ++j) {
              for (int i=is; i<=ie+1; ++i) {
                v_SE(n,k,j,i) = 0.5*(v_SE(n,k,j,i) + vl_temp(n,k,j,i));
                v_SW(n,k,j,i) = 0.5*(v_SW(n,k,j,i) + vr_temp(n,k,j,i));
              }
            }
          }
        }
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je+1; ++j) {
            // Limited transverse reconstructions: call PPMx1() for single-state b_y
            pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, b2,
                                              by_E, by_W, 0, 0, 0);
          }
        }
        // Repeat calculation of x2 edge-centered wavespeeds as in HLL solver
        Real wli[NWAVE], wri[NWAVE];
        int ivx = IVY;
        int ivy = IVX + ((ivx-IVX)+1) % 3;
        int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=js; j<=je+1; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              //--- Load L/R states into local variables
              // UCT with face-centered quantities
              wli[IDN] = wl_fc_(IDN,k,j,i);
              wli[IVX] = wl_fc_(ivx,k,j,i);
              wli[IVY] = wl_fc_(ivy,k,j,i);
              wli[IVZ] = wl_fc_(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wli[IPR] = wl_fc_(IPR,k,j,i);
              wli[IBY] = wl_fc_(IBY,k,j,i);
              wli[IBZ] = wl_fc_(IBZ,k,j,i);

              wri[IDN] = wr_fc_(IDN,k,j,i);
              wri[IVX] = wr_fc_(ivx,k,j,i);
              wri[IVY] = wr_fc_(ivy,k,j,i);
              wri[IVZ] = wr_fc_(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wri[IPR] = wr_fc_(IPR,k,j,i);
              wri[IBY] = wr_fc_(IBY,k,j,i);
              wri[IBZ] = wr_fc_(IBZ,k,j,i);
              Real bxi =  b2_fc(k,j,i);

              Real cl = pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr = pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al = std::min((wri[IVX]-cr),(wli[IVX] - cl));
              Real ar = std::max((wli[IVX] + cl),(wri[IVX] + cr));
              Real bp = ar > 0.0 ? ar : 0.0;
              Real bm = al < 0.0 ? al : 0.0;
              pmb->pfield->alpha_plus_x2_(k,j,i) = bp;
              pmb->pfield->alpha_minus_x2_(k,j,i) = bm;
            }
          }
        }
        // Compute states for 3D UCT
        if (pmb->block_size.nx3 > 1) {
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je+1; ++j) {
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wl3d_,
                                                v_L3L2, v_R3L2, IVY, IVY, 1);
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wl3d_,
                                                v_L3L2, v_R3L2, IVZ, IVZ, 2);
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wr3d_,
                                                v_L3R2, v_R3R2, IVY, IVY, 1);
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wr3d_,
                                                v_L3R2, v_R3R2, IVZ, IVZ, 2);
              // Limited transverse reconstructions: call PPMx3() for single-state b_x
              pmb->precon->PiecewiseParabolicX3(k, j, is, ie, b2,
                                                by_L3, by_R3, 0, 0, 0);
            }
          }
        } // end UCT if 3D
      } // end if (order == 4) UCT4x1
//       else {
//         // compute weights for GS07 CT algorithm
//         for (int k=kl_buf; k<=ku_buf; ++k) {
//           for (int j=js; j<=je+1; ++j) {
//             pmb->pcoord->CenterWidth2(k,j,il_buf,iu_buf,dxw);
// #pragma omp simd
//             for (int i=il_buf; i<=iu_buf; ++i) {
//               Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x2flux(IDN,k,j,i)
//                   / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
//               Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
//               w_x2f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
//             }
//           }
//         }
//       } // end GS07 CT algorithm
#endif // MAGNETIC_FIELDS_ENABLED
  } // if 2D or 3D

  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = flux[X3DIR];
    // set the loop limits
    il=is, iu=ie, jl=js, ju=je;

    if (order != 4) {
      if (MAGNETIC_FIELDS_ENABLED)
        il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    } else { // fourth-order corrections
      if (MAGNETIC_FIELDS_ENABLED) {
        il=is-3, iu=ie+3, jl=js-3, ju=je+3;
      } else { // fourth-order hydro has same reqs as second-order MHD
        il=is-1, iu=ie+1, jl=js-1, ju=je+1;
      }
    } // end order==4

    // Set transverse loop limits for quantities calculated to full accuracy
    il_buf=il, iu_buf=iu, kl_buf=kl, ku_buf=ku;
    if (order == 4) {
      il_buf=il+1, iu_buf=iu-1, jl_buf=jl+1, ju_buf=ju-1;
    }

    // TODO(felker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED)
    // il = is-1, iu = ie+1, jl = js-1, ju = je+1;

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      }
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, w, bcc, wlb_, wr_);
        }

        pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, x3flux, dxw_);
#else  // MHD:
        // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
        // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
        RiemannSolver(k, j, il, iu, IVZ, b3, wl_, wr_, x3flux, e2x3, e1x3, w_x3f, dxw_);
#endif
        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              wl3d_(n,k,j,i) = wl_(n,i);
              wr3d_(n,k,j,i) = wr_(n,i);
            }
          }
        }
        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
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
            pmb->pcoord->LaplacianX3(wl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX3(wr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
              // KGF: temporary load of x1-sliced arrays into 3D arrays for MHD4
              if (MAGNETIC_FIELDS_ENABLED) {
                wl_fc_(n,k,j,i) = wl_(n,i);
                wr_fc_(n,k,j,i) = wr_(n,i);
              }
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            pmb->peos->ApplyPrimitiveFloors(wl_, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_, k, j, i);
          }

          // Compute x3 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e2x3,e1x3 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, flux_fc, dxw_);
#else  // MHD:
          RiemannSolver(k, j, il, iu, IVZ, b3_fc, wl_, wr_, flux_fc,
                        e2x3, e1x3, // flux_fc_IBY, flux_fc_IBZ
                        w_x3f, dxw_);
#endif
          // Apply Laplacian of second-order accurate face-averaged flux on x3 faces
          for (int n=0; n<NWAVE; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++) {
              x3flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
              if (n == IDN && NSCALARS > 0) {
                pmb->pscalars->mass_flux_fc[X3DIR](k,j,i) = flux_fc(n,k,j,i);
              }
            }
          }
        }
      }
    } // end if (order == 4)
    //------------------------------------------------------------------------------
    // end x3 fourth-order hydro and MHD
#if MAGNETIC_FIELDS_ENABLED
    if (order == 4) {
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
      //-------- begin fourth-order upwind constrained transport (UCT4x3)
        // Limited transverse reconstructions: call PPMx1() for x1x3 interfaces
        // L3 states:
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wl3d_,
                                             vl_temp, vr_temp, IVX, IVX, 0);
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wl3d_,
                                             vl_temp, vr_temp, IVZ, IVZ, 2);
        }
      }
        // Store temporary arrays as average of R_x[R_z[]] and R_z[R_x[]] reconstructions
        for (int n=0; n<3; n+=2) { // n=0, 2
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie+1; ++i) {
                v_L3L1(n,k,j,i) = 0.5*(v_L3L1(n,k,j,i) + vl_temp(n,k,j,i));
                v_L3R1(n,k,j,i) = 0.5*(v_L3R1(n,k,j,i) + vr_temp(n,k,j,i));
              }
            }
          }
        }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {

        // R3 states:
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wr3d_,
                                             vl_temp, vr_temp, IVX, IVX, 0);
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, wr3d_,
                                             vl_temp, vr_temp, IVZ, IVZ, 2);
        // Limited transverse reconstructions: call PPM for single-state b_z
        // PPMx1()
        pmb->precon->PiecewiseParabolicX1(k, j, is, ie+1, b3,
                                             bz_L1, bz_R1, 0, 0, 0);

        }
      }
        // Store temporary arrays as average of R_x[R_z[]] and R_z[R_x[]] reconstructions
        for (int n=0; n<3; n+=2) { // n=0, 2
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie+1; ++i) {
                v_R3L1(n,k,j,i) = 0.5*(v_R3L1(n,k,j,i) + vl_temp(n,k,j,i));
                v_R3R1(n,k,j,i) = 0.5*(v_R3R1(n,k,j,i) + vr_temp(n,k,j,i));
              }
            }
          }
        }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je+1; ++j) {
        // Limited transverse reconstructions: call PPMx2() for x2x3 interfaces
        // L3 states:
        pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wl3d_,
                                             vl_temp, vr_temp, IVY, IVY, 1 );
        pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wl3d_,
                                             vl_temp, vr_temp, IVZ, IVZ, 2);
        }
      }
        // Store temporary arrays as average of R_y[R_z[]] and R_z[R_y[]] reconstructions
        for (int n=1; n<3; ++n) {
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je+1; ++j) {
              for (int i=is; i<=ie; ++i) {
                v_L3L2(n,k,j,i) = 0.5*(v_L3L2(n,k,j,i) + vl_temp(n,k,j,i));
                v_L3R2(n,k,j,i) = 0.5*(v_L3R2(n,k,j,i) + vr_temp(n,k,j,i));
              }
            }
          }
        }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je+1; ++j) {

        // R3 states:
        pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wr3d_,
                                             vl_temp, vr_temp, IVY, IVY, 1);
        pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wr3d_,
                                             vl_temp, vr_temp, IVZ, IVZ, 2);
        // PPMx2()
        pmb->precon->PiecewiseParabolicX2(k, j, is, ie, b3, bz_L2, bz_R2, 0, 0, 0);

        }
      }
        // Store temporary arrays as average of R_y[R_z[]] and R_z[R_y[]] reconstructions
        for (int n=1; n<3; ++n) {
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je+1; ++j) {
              for (int i=is; i<=ie; ++i) {
                v_R3L2(n,k,j,i) = 0.5*(v_R3L2(n,k,j,i) + vl_temp(n,k,j,i));
                v_R3R2(n,k,j,i) = 0.5*(v_R3R2(n,k,j,i) + vr_temp(n,k,j,i));
              }
            }
          }
        }

        // Repeat calculation of x3 edge-centered wavespeeds as in HLL solver
        Real wli[NWAVE], wri[NWAVE];
        int ivx = IVZ;
        int ivy = IVX + ((ivx-IVX)+1) % 3;
        int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=il_buf; j<=ju_buf; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              //--- Load L/R states into local variables
              // UCT with face-centered quantities
              wli[IDN] = wl_fc_(IDN,k,j,i);
              wli[IVX] = wl_fc_(ivx,k,j,i);
              wli[IVY] = wl_fc_(ivy,k,j,i);
              wli[IVZ] = wl_fc_(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wli[IPR] = wl_fc_(IPR,k,j,i);
              wli[IBY] = wl_fc_(IBY,k,j,i);
              wli[IBZ] = wl_fc_(IBZ,k,j,i);

              wri[IDN] = wr_fc_(IDN,k,j,i);
              wri[IVX] = wr_fc_(ivx,k,j,i);
              wri[IVY] = wr_fc_(ivy,k,j,i);
              wri[IVZ] = wr_fc_(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wri[IPR] = wr_fc_(IPR,k,j,i);
              wri[IBY] = wr_fc_(IBY,k,j,i);
              wri[IBZ] = wr_fc_(IBZ,k,j,i);

              Real bxi = b3_fc(k,j,i);

              Real cl = pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr = pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al = std::min((wri[IVX]-cr),(wli[IVX] - cl));
              Real ar = std::max((wli[IVX] + cl),(wri[IVX] + cr));
              Real bp = ar > 0.0 ? ar : 0.0;
              Real bm = al < 0.0 ? al : 0.0;
              pmb->pfield->alpha_plus_x3_(k,j,i) = bp;
              pmb->pfield->alpha_minus_x3_(k,j,i) = bm;
            }
          }
        }
      }
//     else { // end if (order == 4) UCT4x1
//         // compute weights for GS07 CT algorithm
//         for (int k=ks; k<=ke+1; ++k) {
//           for (int j=jl_buf; j<=ju_buf; ++j) {
//             pmb->pcoord->CenterWidth3(k,j,il_buf,iu_buf,dxw);
// #pragma omp simd
//             for (int i=il_buf; i<=iu_buf; ++i) {
//               Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x3flux(IDN,k,j,i)
//                   / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
//               Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
//               w_x3f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
//             }
//           }
//         }
//       } // end if GS07 second-order MHD in x3
#endif  // MAGNETIC_FIELDS_ENABLED
  } // end if 3D

  if (SELF_GRAVITY_ENABLED) AddGravityFlux(); // add gravity flux directly

  if (!STS_ENABLED) { // add diffusion fluxes
    AddDiffusionFluxes();
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes_STS
//  \brief Calculate Hydrodynamic Diffusion Fluxes for STS

void Hydro::CalculateFluxes_STS() {
  AddDiffusionFluxes();
}

void Hydro::AddDiffusionFluxes() {
  Field *pf = pmy_block->pfield;
  // add diffusion fluxes
  if (hdif.hydro_diffusion_defined) {
    if (hdif.nu_iso > 0.0 || hdif.nu_aniso > 0.0)
      hdif.AddDiffusionFlux(hdif.visflx,flux);
    if (NON_BAROTROPIC_EOS) {
      if (hdif.kappa_iso > 0.0 || hdif.kappa_aniso > 0.0)
        hdif.AddDiffusionEnergyFlux(hdif.cndflx,flux);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED && NON_BAROTROPIC_EOS) {
    if (pf->fdif.field_diffusion_defined)
      pf->fdif.AddPoyntingFlux(pf->fdif.pflux);
  }
  return;
}
