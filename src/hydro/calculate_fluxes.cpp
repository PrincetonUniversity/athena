//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//! \brief Calculate hydro/MHD fluxes

// C headers
#include <iostream>   // endl
#include <sstream>    // stringstream

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
//! \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, const int order) {
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

  // -------------- fourth-order MHD (UCT4) indices and variables ------
  // transverse loop limits for quantities calculated to full accuracy
  int il_buf, iu_buf, jl_buf, ju_buf, kl_buf, ku_buf;
  // fourth-order approximations to face-centered (point-valued) magnetic fields,
  // computed in Field::FaceAveragedToCellAveragedField before this function
  AthenaArray<Real> &b1_fc = pmb->pfield->b_fc.x1f, &b2_fc = pmb->pfield->b_fc.x2f,
                    &b3_fc = pmb->pfield->b_fc.x3f;
  // corner-reconstructed states for the UCT4 corner EMF (see ComputeCornerE_UCT4):
  // (in x2) of the x1-face L/R velocity states and single-state b1
  AthenaArray<Real> &v_SE = pmb->pfield->v_SE, &v_NE = pmb->pfield->v_NE,
                    &v_NW = pmb->pfield->v_NW, &v_SW = pmb->pfield->v_SW;
  AthenaArray<Real> &bx_N = pmb->pfield->bx_N, &bx_S = pmb->pfield->bx_S;
  // (in x1) of the x2-face single-state b2
  AthenaArray<Real> &by_E = pmb->pfield->by_E, &by_W = pmb->pfield->by_W;
  // 3D corner states
  AthenaArray<Real> &bz_R1 = pmb->pfield->bz_R1, &bz_L1 = pmb->pfield->bz_L1,
                    &bz_R2 = pmb->pfield->bz_R2, &bz_L2 = pmb->pfield->bz_L2,
                    &by_R3 = pmb->pfield->by_R3, &by_L3 = pmb->pfield->by_L3,
                    &bx_R3 = pmb->pfield->bx_R3, &bx_L3 = pmb->pfield->bx_L3;
  AthenaArray<Real> &v_R3R2 = pmb->pfield->v_R3R2, &v_R3L2 = pmb->pfield->v_R3L2,
                    &v_L3R2 = pmb->pfield->v_L3R2, &v_L3L2 = pmb->pfield->v_L3L2,
                    &v_R3R1 = pmb->pfield->v_R3R1, &v_R3L1 = pmb->pfield->v_R3L1,
                    &v_L3R1 = pmb->pfield->v_L3R1, &v_L3L1 = pmb->pfield->v_L3L1;
  // 1D pencil scratch for the corner reconstruction sweeps:
  AthenaArray<Real> &v_SE_ = pmb->pfield->v_SE_, &v_NE_ = pmb->pfield->v_NE_,
                    &v_NW_ = pmb->pfield->v_NW_, &v_SW_ = pmb->pfield->v_SW_;
  AthenaArray<Real> &bx_N_ = pmb->pfield->bx_N_, &bx_S_ = pmb->pfield->bx_S_;
  AthenaArray<Real> &by_E_ = pmb->pfield->by_E_, &by_W_ = pmb->pfield->by_W_;
  AthenaArray<Real> &bz_R1_ = pmb->pfield->bz_R1_, &bz_L1_ = pmb->pfield->bz_L1_,
                    &bz_R2_ = pmb->pfield->bz_R2_, &bz_L2_ = pmb->pfield->bz_L2_,
                    &by_R3_ = pmb->pfield->by_R3_, &by_L3_ = pmb->pfield->by_L3_,
                    &bx_R3_ = pmb->pfield->bx_R3_, &bx_L3_ = pmb->pfield->bx_L3_;
  AthenaArray<Real> &v_R3R2_ = pmb->pfield->v_R3R2_, &v_R3L2_ = pmb->pfield->v_R3L2_,
                    &v_L3R2_ = pmb->pfield->v_L3R2_, &v_L3L2_ = pmb->pfield->v_L3L2_,
                    &v_R3R1_ = pmb->pfield->v_R3R1_, &v_R3L1_ = pmb->pfield->v_R3L1_,
                    &v_L3R1_ = pmb->pfield->v_L3R1_, &v_L3L1_ = pmb->pfield->v_L3L1_;
  AthenaArray<Real> &vl_temp2 = pmb->pfield->vl_temp2_,
                    &vr_temp2 = pmb->pfield->vr_temp2_;
  // "current-row" swap partners for the j-pencil decompositions
  AthenaArray<Real> &v_NEb_ = pmb->pfield->v_NEb_, &v_NWb_ = pmb->pfield->v_NWb_,
                    &bx_Nb_ = pmb->pfield->bx_Nb_;
  // scratch for the 2D 4th-order face-averaged EMF corrections (reused from Field)
  AthenaArray<Real> &laplacian_e2x1 = pmb->pfield->scr1_kji_x1fc_,
                    &laplacian_e1x2 = pmb->pfield->scr2_kji_x2fc_;
  // -------------- end fourth-order MHD indices and variables ------
#endif
  AthenaArray<Real> &flux_fc = scr1_nkji_;
  AthenaArray<Real> &laplacian_all_fc = scr2_nkji_;

  //--------------------------------------------------------------------------------------
  // i-direction

  AthenaArray<Real> &x1flux = flux[X1DIR];
  // set the loop limits
  jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED || order == 4) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl = js-1, ju = je+1, kl = ks, ku = ke;
      else // 3D
        jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
    }
  }
  // fourth-order MHD (UCT4): the transverse corner reconstructions of the x1-face
  // Riemann states require 2 additional rows of face states in each transverse direction
  if (MAGNETIC_FIELDS_ENABLED && order == 4) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl = js-3, ju = je+3, kl = ks, ku = ke;
      else // 3D
        jl = js-3, ju = je+3, kl = ks-3, ku = ke+3;
    }
  }
#if MAGNETIC_FIELDS_ENABLED
  // transverse limits of the rows on which the face-centered states/wavespeeds are
  // valid after the Laplacian correction (shrinks by 1 at each transverse edge)
  jl_buf = jl, ju_buf = ju, kl_buf = kl, ku_buf = ku;
  if (order == 4) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl_buf += 1, ju_buf -= 1;
      else // 3D
        jl_buf += 1, ju_buf -= 1, kl_buf += 1, ku_buf -= 1;
    }
  }
#endif

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else {
        if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
          pmb->precon->PiecewiseParabolicFastX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
          pmb->precon->WENOZX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
          pmb->precon->WENOMZX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        else
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

#if MAGNETIC_FIELDS_ENABLED
    // 2D UCT4: e2_x1f and e1_x2f are consumed directly as the edge-averaged EMFs in the
    // b.x3f CT update; correct them to 4th-order face averages like the hydro fluxes.
    // Compute the Laplacian of the face-averaged EMFY before it is overwritten below.
    if (pmb->block_size.nx2 > 1 && pmb->block_size.nx3 == 1) {
      pmb->pcoord->LaplacianX1All(e2x1, laplacian_e2x1, 0, 0,
                                  kl_buf, ku_buf, jl_buf, ju_buf, is, ie+1);
    }
#endif

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
            // cache face-centered L/R states for the UCT4 corner wavespeed estimates
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
        pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, flux_fc, dxw_);
#else  // MHD: (pass face-centered B1 for the point-valued Riemann problem)
        RiemannSolver(k, j, is, ie+1, IVX, b1_fc, wl_, wr_, flux_fc, e3x1, e2x1,
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

#if MAGNETIC_FIELDS_ENABLED
    // 2D UCT4: apply the 4th-order face-averaged correction to EMFY on x1 faces
    if (pmb->block_size.nx2 > 1 && pmb->block_size.nx3 == 1) {
      for (int j=jl_buf; j<=ju_buf; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          e2x1(ks,j,i) += C*laplacian_e2x1(ks,j,i);
        }
      }
    }
#endif
  } // end if (order == 4)
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro

#if MAGNETIC_FIELDS_ENABLED
  //-------- begin fourth-order upwind constrained transport (UCT4x1)
  if (order == 4) {
    // 1D domains work via trivial copying of the fluid fluxes in ComputeCornerE_UCT4;
    // the corner reconstructions require at least a 2D domain
    if (pmb->block_size.nx2 > 1) {
      // Unlike standard Athena++ E_z^c upwinding, which requires loading
      // [is-1:ie+1] x [js-1:je+1] face states, the UCT corner states are centered on
      // the corner, so only the real range including the uppermost corner is required:
      // [is:ie+1] x [js:je+1]

      // Limited transverse reconstructions: call PPMx2() on the x1-face L/R Riemann
      // states. Pencil decomposition in j: for the corner at x2_{j-1/2} (index j), the
      // below/N/L2 state is the PPM ql of row j-1 and the above/S/R2 state is the PPM
      // qr of row j; ql of the current row is buffered and swapped for the next row.
      // wl_{i-1/2} is the E (L1) side of the interface
      for (int k=ks; k<=ke; ++k) {
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie+1, wl3d_, v_NE_, v_SE_,
                                          IVX, IVY, 0);
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie+1, wr3d_, v_NW_, v_SW_,
                                          IVX, IVY, 0);
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie+1, b1, bx_N_, bx_S_,
                                          0, 0, 0);
        for (int j=js; j<=je+1; ++j) {
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, wl3d_, v_NEb_, v_SE_,
                                            IVX, IVY, 0);
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, wr3d_, v_NWb_, v_SW_,
                                            IVX, IVY, 0);
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie+1, b1, bx_Nb_, bx_S_,
                                            0, 0, 0);
          for (int n=0; n<2; n++) {
            for (int i=is; i<=ie+1; i++) {
              v_NE(n,k,j,i) = v_NE_(n,i);
              v_SE(n,k,j,i) = v_SE_(n,i);
              v_NW(n,k,j,i) = v_NW_(n,i);
              v_SW(n,k,j,i) = v_SW_(n,i);
            }
          }
          for (int i=is; i<=ie+1; i++) {
            bx_N(k,j,i) = bx_N_(i);
            bx_S(k,j,i) = bx_S_(i);
          }
          v_NE_.SwapAthenaArray(v_NEb_);
          v_NW_.SwapAthenaArray(v_NWb_);
          bx_N_.SwapAthenaArray(bx_Nb_);
        } // end of loop over j
      } // end of loop over k

      // Repeat calculation of x1 face-centered wavespeeds as in the HLL solver
      {
        Real wli[NWAVE], wri[NWAVE];
        const int ivx = IVX;
        const int ivy = IVX + ((ivx-IVX)+1) % 3;
        const int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=jl_buf; j<=ju_buf; ++j) {
            for (int i=is; i<=ie+1; ++i) {
              //--- Load face-centered L/R states into local variables
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
              Real bxi = b1_fc(k,j,i);

              Real cl = pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr = pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al = std::min((wri[IVX]-cr), (wli[IVX]-cl));
              Real ar = std::max((wli[IVX]+cl), (wri[IVX]+cr));
              pmb->pfield->alpha_plus_x1_(k,j,i) = ar > 0.0 ? ar : 0.0;
              pmb->pfield->alpha_minus_x1_(k,j,i) = al < 0.0 ? al : 0.0;
            }
          }
        }
      }

      // Compute 3D corner states: PPMx3() of the x1-face L/R states and b1.
      // Pencil decomposition in k: for the corner at x3_{k-1/2} (index k), the L3 state
      // is the PPM ql of row k-1 (stored with a +1 row offset below) and the R3 state
      // is the PPM qr of row k.
      if (pmb->block_size.nx3 > 1) {
        for (int k=ks-1; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wl3d_, v_L3L1_, v_R3L1_,
                                              IVX, IVX, 0);
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wl3d_, v_L3L1_, v_R3L1_,
                                              IVZ, IVZ, 2);
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wr3d_, v_L3R1_, v_R3R1_,
                                              IVX, IVX, 0);
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, wr3d_, v_L3R1_, v_R3R1_,
                                              IVZ, IVZ, 2);
            // Limited transverse reconstructions: call PPMx3() for single-state b_x
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie+1, b1, bx_L3_, bx_R3_,
                                              0, 0, 0);
            for (int n=0; n<3; n+=2) {
              for (int i=is; i<=ie+1; i++) {
                v_L3L1(n,k+1,j,i) = v_L3L1_(n,i);
                v_L3R1(n,k+1,j,i) = v_L3R1_(n,i);
                v_R3L1(n,k,j,i) = v_R3L1_(n,i);
                v_R3R1(n,k,j,i) = v_R3R1_(n,i);
              }
            }
            for (int i=is; i<=ie+1; i++) {
              bx_L3(k+1,j,i) = bx_L3_(i);
              bx_R3(k,j,i) = bx_R3_(i);
            }
          }
        }
      } // end UCT if 3D
    } // end if 2D or 3D
  }  // end if (order == 4) UCT4x1
#endif  // MAGNETIC_FIELDS_ENABLED

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = flux[X2DIR];
    // set the loop limits
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED || order == 4) {
      if (pmb->block_size.nx3 == 1) // 2D
        kl = ks, ku = ke;
      else // 3D
        kl = ks-1, ku = ke+1;
    }
    // fourth-order MHD (UCT4): extended transverse rows for corner reconstructions
    if (MAGNETIC_FIELDS_ENABLED && order == 4) {
      if (pmb->block_size.nx3 == 1) // 2D
        il = is-3, iu = ie+3, kl = ks, ku = ke;
      else // 3D
        il = is-3, iu = ie+3, kl = ks-3, ku = ke+3;
    }
#if MAGNETIC_FIELDS_ENABLED
    il_buf = il, iu_buf = iu, kl_buf = kl, ku_buf = ku;
    if (order == 4) {
      if (pmb->block_size.nx3 == 1) // 2D
        il_buf += 1, iu_buf -= 1;
      else // 3D
        il_buf += 1, iu_buf -= 1, kl_buf += 1, ku_buf -= 1;
    }
#endif

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else {
        if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
          pmb->precon->PiecewiseParabolicFastX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
          pmb->precon->WENOZX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
          pmb->precon->WENOMZX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        else
          pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      }
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
            pmb->precon->PiecewiseParabolicFastX2(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
            pmb->precon->WENOZX2(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
            pmb->precon->WENOMZX2(k, j, il, iu, w, bcc, wlb_, wr_);
          else
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

#if MAGNETIC_FIELDS_ENABLED
      // 2D UCT4: Laplacian of the face-averaged EMFX before it is overwritten below
      if (pmb->block_size.nx3 == 1) {
        pmb->pcoord->LaplacianX2All(e1x2, laplacian_e1x2, 0, 0,
                                    kl_buf, ku_buf, js, je+1, il_buf, iu_buf);
      }
#endif

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
              // cache face-centered L/R states for the UCT4 corner wavespeed estimates
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
          pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVY, wl_, wr_, flux_fc, dxw_);
#else  // MHD: (pass face-centered B2 for the point-valued Riemann problem)
          RiemannSolver(k, j, il, iu, IVY, b2_fc, wl_, wr_, flux_fc, e1x2, e3x2,
                        w_x2f, dxw_);
#endif

          // Apply Laplacian of second-order accurate face-averaged flux on x2 faces
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

#if MAGNETIC_FIELDS_ENABLED
      // 2D UCT4: apply the 4th-order face-averaged correction to EMFX on x2 faces
      if (pmb->block_size.nx3 == 1) {
        for (int j=js; j<=je+1; ++j) {
#pragma omp simd
          for (int i=il_buf; i<=iu_buf; ++i) {
            e1x2(ks,j,i) += C*laplacian_e1x2(ks,j,i);
          }
        }
      }
#endif
    } // end if (order == 4)
    //------------------------------------------------------------------------------
    // end x2 fourth-order hydro

#if MAGNETIC_FIELDS_ENABLED
    //-------- begin fourth-order upwind constrained transport (UCT4x2)
    if (order == 4) {
      // Limited transverse reconstructions: call PPMx1() on the x2-face L/R Riemann
      // states, and average with the transposed reconstruction ordering from UCT4x1:
      // corner state = 0.5*(R_x1[R_x2[w]] + R_x2[R_x1[w]]).
      // PPMx1 pencil output ql(i) already holds the L1 state of the corner at
      // x1_{i-1/2} (index i), so no swap buffering is needed in x1.
      // wl_{j-1/2} is the N (L2) side of the interface
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, wl3d_,
                                            vl_temp2, vr_temp2, IVX, IVY, 0);
          for (int n=0; n<2; n++) {
            for (int i=is; i<=ie+1; i++) {
              v_NE(n,k,j,i) = 0.5*(v_NE(n,k,j,i) + vl_temp2(n,i));
              v_NW(n,k,j,i) = 0.5*(v_NW(n,k,j,i) + vr_temp2(n,i));
            }
          }
          // wr_{j-1/2} is the S (R2) side of the interface
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, wr3d_,
                                            vl_temp2, vr_temp2, IVX, IVY, 0);
          for (int n=0; n<2; n++) {
            for (int i=is; i<=ie+1; i++) {
              v_SE(n,k,j,i) = 0.5*(v_SE(n,k,j,i) + vl_temp2(n,i));
              v_SW(n,k,j,i) = 0.5*(v_SW(n,k,j,i) + vr_temp2(n,i));
            }
          }
          // Limited transverse reconstructions: call PPMx1() for single-state b_y
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, b2,
                                            by_E_, by_W_, 0, 0, 0);
          for (int i=is; i<=ie+1; i++) {
            by_E(k,j,i) = by_E_(i);
            by_W(k,j,i) = by_W_(i);
          }
        }
      }

      // Repeat calculation of x2 face-centered wavespeeds as in the HLL solver
      {
        Real wli[NWAVE], wri[NWAVE];
        const int ivx = IVY;
        const int ivy = IVX + ((ivx-IVX)+1) % 3;
        const int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=js; j<=je+1; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              //--- Load face-centered L/R states into local variables
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
              Real bxi = b2_fc(k,j,i);

              Real cl = pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr = pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al = std::min((wri[IVX]-cr), (wli[IVX]-cl));
              Real ar = std::max((wli[IVX]+cl), (wri[IVX]+cr));
              pmb->pfield->alpha_plus_x2_(k,j,i) = ar > 0.0 ? ar : 0.0;
              pmb->pfield->alpha_minus_x2_(k,j,i) = al < 0.0 ? al : 0.0;
            }
          }
        }
      }

      // Compute 3D corner states: PPMx3() of the x2-face L/R states and b2, with the
      // same +1 row offset for the ql (L3) states as in UCT4x1
      if (pmb->block_size.nx3 > 1) {
        for (int k=ks-1; k<=ke+1; ++k) {
          for (int j=js; j<=je+1; ++j) {
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wl3d_, v_L3L2_, v_R3L2_,
                                              IVY, IVY, 1);
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wl3d_, v_L3L2_, v_R3L2_,
                                              IVZ, IVZ, 2);
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wr3d_, v_L3R2_, v_R3R2_,
                                              IVY, IVY, 1);
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie, wr3d_, v_L3R2_, v_R3R2_,
                                              IVZ, IVZ, 2);
            // Limited transverse reconstructions: call PPMx3() for single-state b_y
            pmb->precon->PiecewiseParabolicX3(k, j, is, ie, b2, by_L3_, by_R3_,
                                              0, 0, 0);
            for (int n=1; n<3; n++) {
              for (int i=is; i<=ie; i++) {
                v_L3L2(n,k+1,j,i) = v_L3L2_(n,i);
                v_L3R2(n,k+1,j,i) = v_L3R2_(n,i);
                v_R3L2(n,k,j,i) = v_R3L2_(n,i);
                v_R3R2(n,k,j,i) = v_R3R2_(n,i);
              }
            }
            for (int i=is; i<=ie; i++) {
              by_L3(k+1,j,i) = by_L3_(i);
              by_R3(k,j,i) = by_R3_(i);
            }
          }
        }
      } // end UCT if 3D
    } // end if (order == 4) UCT4x2
#endif // MAGNETIC_FIELDS_ENABLED
  }

  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = flux[X3DIR];
    // set the loop limits
    il = is, iu = ie, jl = js, ju = je;
    if (MAGNETIC_FIELDS_ENABLED || order == 4) {
      il = is-1, iu = ie+1, jl = js-1, ju = je+1;
    }
    // fourth-order MHD (UCT4): extended transverse rows for corner reconstructions
    if (MAGNETIC_FIELDS_ENABLED && order == 4) {
      il = is-3, iu = ie+3, jl = js-3, ju = je+3;
    }
#if MAGNETIC_FIELDS_ENABLED
    il_buf = il, iu_buf = iu, jl_buf = jl, ju_buf = ju;
    if (order == 4) {
      il_buf += 1, iu_buf -= 1, jl_buf += 1, ju_buf -= 1;
    }
#endif

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else {
        if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
          pmb->precon->PiecewiseParabolicFastX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
          pmb->precon->WENOZX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
          pmb->precon->WENOMZX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        else
          pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      }
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
            pmb->precon->PiecewiseParabolicFastX3(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
            pmb->precon->WENOZX3(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
            pmb->precon->WENOMZX3(k, j, il, iu, w, bcc, wlb_, wr_);
          else
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
              // cache face-centered L/R states for the UCT4 corner wavespeed estimates
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
          pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, flux_fc, dxw_);
#else  // MHD: (pass face-centered B3 for the point-valued Riemann problem)
          RiemannSolver(k, j, il, iu, IVZ, b3_fc, wl_, wr_, flux_fc, e2x3, e1x3,
                        w_x3f, dxw_);
#endif
          // Apply Laplacian of second-order accurate face-averaged flux on x3 faces
          for (int n=0; n<NHYDRO; ++n) {
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
    // end x3 fourth-order hydro

#if MAGNETIC_FIELDS_ENABLED
    //-------- begin fourth-order upwind constrained transport (UCT4x3)
    if (order == 4) {
      // (a) Limited transverse reconstructions: PPMx1() of the x3-face L/R Riemann
      // states and single-state b3, averaged with the x3-of-x1 reconstructions from
      // UCT4x1. PPMx1 pencil output ql(i) holds the L1 state at x1_{i-1/2} directly.
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          // wl_{k-1/2} is the L3 side of the interface
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, wl3d_,
                                            vl_temp2, vr_temp2, IVX, IVX, 0);
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, wl3d_,
                                            vl_temp2, vr_temp2, IVZ, IVZ, 2);
          for (int n=0; n<3; n+=2) {
            for (int i=is; i<=ie+1; i++) {
              v_L3L1(n,k,j,i) = 0.5*(v_L3L1(n,k,j,i) + vl_temp2(n,i));
              v_L3R1(n,k,j,i) = 0.5*(v_L3R1(n,k,j,i) + vr_temp2(n,i));
            }
          }
          // wr_{k-1/2} is the R3 side of the interface
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, wr3d_,
                                            vl_temp2, vr_temp2, IVX, IVX, 0);
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, wr3d_,
                                            vl_temp2, vr_temp2, IVZ, IVZ, 2);
          for (int n=0; n<3; n+=2) {
            for (int i=is; i<=ie+1; i++) {
              v_R3L1(n,k,j,i) = 0.5*(v_R3L1(n,k,j,i) + vl_temp2(n,i));
              v_R3R1(n,k,j,i) = 0.5*(v_R3R1(n,k,j,i) + vr_temp2(n,i));
            }
          }
          // Limited transverse reconstructions: call PPMx1() for single-state b_z
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, b3,
                                            bz_L1_, bz_R1_, 0, 0, 0);
          for (int i=is; i<=ie+1; i++) {
            bz_L1(k,j,i) = bz_L1_(i);
            bz_R1(k,j,i) = bz_R1_(i);
          }
        }
      }

      // (b) Limited transverse reconstructions: PPMx2() of the x3-face L/R Riemann
      // states and single-state b3, averaged with the x3-of-x2 reconstructions from
      // UCT4x2. Pencil decomposition in j with swap buffering, as in UCT4x1.
      for (int k=ks; k<=ke+1; ++k) {
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie, wl3d_, v_L3L2_, v_L3R2_,
                                          IVY, IVY, 1);
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie, wl3d_, v_L3L2_, v_L3R2_,
                                          IVZ, IVZ, 2);
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie, wr3d_, v_R3L2_, v_R3R2_,
                                          IVY, IVY, 1);
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie, wr3d_, v_R3L2_, v_R3R2_,
                                          IVZ, IVZ, 2);
        pmb->precon->PiecewiseParabolicX2(k, js-1, is, ie, b3, bz_L2_, bz_R2_,
                                          0, 0, 0);
        for (int j=js; j<=je+1; ++j) {
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wl3d_, v_NEb_, v_L3R2_,
                                            IVY, IVY, 1);
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wl3d_, v_NEb_, v_L3R2_,
                                            IVZ, IVZ, 2);
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wr3d_, v_NWb_, v_R3R2_,
                                            IVY, IVY, 1);
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie, wr3d_, v_NWb_, v_R3R2_,
                                            IVZ, IVZ, 2);
          pmb->precon->PiecewiseParabolicX2(k, j, is, ie, b3, bx_Nb_, bz_R2_,
                                            0, 0, 0);
          for (int n=1; n<3; n++) {
            for (int i=is; i<=ie; i++) {
              v_L3L2(n,k,j,i) = 0.5*(v_L3L2(n,k,j,i) + v_L3L2_(n,i));
              v_L3R2(n,k,j,i) = 0.5*(v_L3R2(n,k,j,i) + v_L3R2_(n,i));
              v_R3L2(n,k,j,i) = 0.5*(v_R3L2(n,k,j,i) + v_R3L2_(n,i));
              v_R3R2(n,k,j,i) = 0.5*(v_R3R2(n,k,j,i) + v_R3R2_(n,i));
            }
          }
          for (int i=is; i<=ie; i++) {
            bz_L2(k,j,i) = bz_L2_(i);
            bz_R2(k,j,i) = bz_R2_(i);
          }
          v_L3L2_.SwapAthenaArray(v_NEb_);
          v_R3L2_.SwapAthenaArray(v_NWb_);
          bz_L2_.SwapAthenaArray(bx_Nb_);
        } // end of loop over j
      } // end of loop over k

      // Repeat calculation of x3 face-centered wavespeeds as in the HLL solver
      {
        Real wli[NWAVE], wri[NWAVE];
        const int ivx = IVZ;
        const int ivy = IVX + ((ivx-IVX)+1) % 3;
        const int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=jl_buf; j<=ju_buf; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              //--- Load face-centered L/R states into local variables
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
              Real al = std::min((wri[IVX]-cr), (wli[IVX]-cl));
              Real ar = std::max((wli[IVX]+cl), (wri[IVX]+cr));
              pmb->pfield->alpha_plus_x3_(k,j,i) = ar > 0.0 ? ar : 0.0;
              pmb->pfield->alpha_minus_x3_(k,j,i) = al < 0.0 ? al : 0.0;
            }
          }
        }
      }
    } // end if (order == 4) UCT4x3
#endif  // MAGNETIC_FIELDS_ENABLED
  }

  if (!STS_ENABLED)
    AddDiffusionFluxes();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes_STS
//! \brief Calculate Hydrodynamic Diffusion Fluxes for STS

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
