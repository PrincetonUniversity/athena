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

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b, FaceField &b_fc,
                            AthenaArray<Real> &bcc, AthenaArray<Real> &bcc_center,
                            const int order) {
  MeshBlock *pmb=pmy_block;
  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  // fourth-order indices and variables:
  int il_buf, iu_buf, jl_buf, ju_buf, kl_buf, ku_buf;
  AthenaArray<Real> b1_fc, b2_fc, b3_fc;

  AthenaArray<Real> b1, b2, b3, w_x1f, w_x2f, w_x3f, e2x1, e3x1, e1x2, e3x2, e1x3, e2x3;
  if (MAGNETIC_FIELDS_ENABLED) {
    b1.InitWithShallowCopy(b.x1f);
    b2.InitWithShallowCopy(b.x2f);
    b3.InitWithShallowCopy(b.x3f);
    b1_fc.InitWithShallowCopy(b_fc.x1f);
    b2_fc.InitWithShallowCopy(b_fc.x2f);
    b3_fc.InitWithShallowCopy(b_fc.x3f);
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
  // fourth-order MHD quantities:
  AthenaArray<Real> flux_fc_IBY, flux_fc_IBZ;
  flux_fc_IBY.InitWithShallowSlice(flux_fc, 4, IBY, 1);
  flux_fc_IBZ.InitWithShallowSlice(flux_fc, 4, IBZ, 1);
  // Reconstruct (in x2) velocity L/R states along x1 faces to corners
  AthenaArray<Real> v_SE, v_NE, v_NW, v_SW;
  v_SE.InitWithShallowCopy(pmb->pfield->v_SE);
  v_NE.InitWithShallowCopy(pmb->pfield->v_NE);
  v_NW.InitWithShallowCopy(pmb->pfield->v_NW);
  v_SW.InitWithShallowCopy(pmb->pfield->v_SW);
  // Reconstruct (in x2) the single state face-averaged b_x.x1f
  AthenaArray<Real> bx_N,  bx_S;
  bx_N.InitWithShallowCopy(pmb->pfield->bx_N);
  bx_S.InitWithShallowCopy(pmb->pfield->bx_S);
  // Reconstruct (in x1) velocity L/R states along x2 faces to corners
  AthenaArray<Real> vl_temp;
  AthenaArray<Real> vr_temp;
  vl_temp.InitWithShallowCopy(pmb->pfield->vl_temp_);
  vr_temp.InitWithShallowCopy(pmb->pfield->vr_temp_);
  // Reconstruct (in x1) the single state face-averaged b_y.x2f
  AthenaArray<Real> by_E, by_W;
  by_E.InitWithShallowCopy(pmb->pfield->by_E);
  by_W.InitWithShallowCopy(pmb->pfield->by_W);

  // 3D UCT states:
  //  AthenaArray<Real> alpha_plus_x3_, alpha_minus_x3_;
  AthenaArray<Real> v_R3R2, v_R3L2, v_L3R2, v_L3L2;
  AthenaArray<Real> v_R3R1, v_R3L1, v_L3R1, v_L3L1;
  AthenaArray<Real> bz_R1, bz_L1;
  AthenaArray<Real> bz_R2, bz_L2;
  AthenaArray<Real> by_R3, by_L3;
  AthenaArray<Real> bx_R3, bx_L3;
  bz_R1.InitWithShallowCopy(pmb->pfield->bz_R1);
  bz_L1.InitWithShallowCopy(pmb->pfield->bz_L1);
  bz_R2.InitWithShallowCopy(pmb->pfield->bz_R2);
  bz_L2.InitWithShallowCopy(pmb->pfield->bz_L2);
  by_R3.InitWithShallowCopy(pmb->pfield->by_R3);
  by_L3.InitWithShallowCopy(pmb->pfield->by_L3);
  bx_R3.InitWithShallowCopy(pmb->pfield->bx_R3);
  bx_L3.InitWithShallowCopy(pmb->pfield->bx_L3);

  v_R3R2.InitWithShallowCopy(pmb->pfield->v_R3R2);
  v_R3L2.InitWithShallowCopy(pmb->pfield->v_R3L2);
  v_L3R2.InitWithShallowCopy(pmb->pfield->v_L3R2);
  v_L3L2.InitWithShallowCopy(pmb->pfield->v_L3L2);
  v_R3R1.InitWithShallowCopy(pmb->pfield->v_R3R1);
  v_R3L1.InitWithShallowCopy(pmb->pfield->v_R3L1);
  v_L3R1.InitWithShallowCopy(pmb->pfield->v_L3R1);
  v_L3L1.InitWithShallowCopy(pmb->pfield->v_L3L1);
//----------------------------------------------------------------------------------------
// i-direction

  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;

  if (order != 4) {
    if (MAGNETIC_FIELDS_ENABLED) {
      if (pmb->block_size.nx2 > 1) {
        if (pmb->block_size.nx3 == 1) // 2D
          jl=js-1, ju=je+1, kl=ks, ku=ke;
        else // 3D
          jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
      }
    }
  } else { // fourth-order corrections
    if (MAGNETIC_FIELDS_ENABLED) {
      if (pmb->block_size.nx2 > 1) {
        if(pmb->block_size.nx3 == 1) { // 2D
          jl=js-3, ju=je+3, kl=ks, ku=ke;
        } else { // 3D
          jl=js-3, ju=je+3, kl=ks-3, ku=ke+3;
        }
      }
    } else { // fourth-order hydro has same reqs as second-order MHD
      if (pmb->block_size.nx2 > 1) {
        if (pmb->block_size.nx3 == 1) // 2D
          jl=js-1, ju=je+1, kl=ks, ku=ke;
        else // 3D
          jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
      }
    }
  } // end order==4

  // Set transverse loop limits for quantities calculated to full accuracy
  // TODO(kfelker): check fourth-order MHD dependence
  jl_buf=jl, ju_buf=ju, kl_buf=kl, ku_buf=ku; // 1D defaults
  if (order == 4) {
    if(pmb->block_size.nx2 > 1) {
      if(pmb->block_size.nx3 == 1) // 2D
        jl_buf+=1, ju_buf-=1;
      else // 3D
        jl_buf+=1, ju_buf-=1, kl_buf+=1, ku_buf-=1;
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

  // begin x1 fourth-order hydro and MHD:
  //------------------------------------------------------------------------------
  if (order == 4) {
    if (MAGNETIC_FIELDS_ENABLED) {
      // Copy electric fields back into x1flux array
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            x1flux(IBY,k,j,i) = e3x1(k,j,i); // technically -e3x1
            x1flux(IBZ,k,j,i) = e2x1(k,j,i);
          }
        }
      }
    }
    // Compute Laplacian of primitive Riemann states on x1 faces
    pmb->pcoord->LaplacianX1(wl, laplacian_l_fc, is, ie+1, jl, ju, kl, ku, 0, NWAVE-1);
    pmb->pcoord->LaplacianX1(wr, laplacian_r_fc, is, ie+1, jl, ju, kl, ku, 0, NWAVE-1);

    // TODO(kfelker): assuming uniform mesh with dx1f=dx2f=dx3f, so this should factor out
    // TODO(kfelker): also, this may need to be dx1v, since Laplacian is cell-centered
    Real h = pmb->pcoord->dx1f(is);  // pco->dx1f(i); inside loop
    Real C = (h*h)/24.0;

    // Approximate x1 face-centered primitive Riemann states
    for (int n=0; n<NWAVE; ++n) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
          // pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw);
          for (int i=is; i<=ie+1; ++i) {
            wl_fc(n,k,j,i) = wl(n,k,j,i) - C*laplacian_l_fc(n,k,j,i);
            wr_fc(n,k,j,i) = wr(n,k,j,i) - C*laplacian_r_fc(n,k,j,i);
            // reapply primitive variable floors to face-centered L/R Riemann states
            // TODO(kfelker): only needs to be called 1x for all NWAVE
            pmb->peos->ApplyPrimitiveFloors(wl_fc, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_fc, k, j, i);
          }
        }
      }
    }

    // Compute x1 interface fluxes from face-centered primitive variables
    RiemannSolver(kl, ku, jl, ju, is, ie+1, IVX, b1_fc, wl_fc, wr_fc, flux_fc,
                  flux_fc_IBY, flux_fc_IBZ);

    // Compute Laplacian of second-order accurate face-averaged flux on x1 faces
    pmb->pcoord->LaplacianX1(x1flux, laplacian_l_fc, is, ie+1, jl, ju, kl, ku,
                             0, NWAVE-1);

    // Correct face-averaged fluxes (Guzik eq. 10)
    for(int n=0; n<NWAVE; n++) {
      for (int k=kl_buf; k<=ku_buf; ++k) {
        for (int j=jl_buf; j<=ju_buf; ++j) {
          // pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw);
          for(int i=is; i<=ie+1; i++) {
            x1flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_l_fc(n,k,j,i);
          }
        }
      }
    }
    if (MAGNETIC_FIELDS_ENABLED) {
      // Copy high-order magnetic x1flux back into electric field arrays
      for (int k=kl_buf; k<=ku_buf; ++k) {
        for (int j=jl_buf; j<=ju_buf; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            e3x1(k,j,i) = x1flux(IBY,k,j,i);
            e2x1(k,j,i) = x1flux(IBZ,k,j,i);
          }
        }
      }
    }

  } // end if (order == 4)
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro and MHD

  if (MAGNETIC_FIELDS_ENABLED) {
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
        pmb->precon->PiecewiseParabolicUCTx2(pmb, ks, ke, js, je+1, is, ie+1, wl, IVX, 0,
                                             v_NE, v_SE);
        pmb->precon->PiecewiseParabolicUCTx2(pmb, ks, ke, js, je+1, is, ie+1, wl, IVY, 1,
                                             v_NE, v_SE);
        // wr_{i-1/2}  is W side of interface --> discontinuous states in x2, L=N,  R=S
        pmb->precon->PiecewiseParabolicUCTx2(pmb, ks, ke, js, je+1, is, ie+1, wr, IVX, 0,
                                             v_NW, v_SW);
        pmb->precon->PiecewiseParabolicUCTx2(pmb, ks, ke, js, je+1, is, ie+1, wr, IVY, 1,
                                             v_NW, v_SW);
        // Limited transverse reconstructions: call PPMx2() for single-state b_x
        pmb->precon->PiecewiseParabolicUCTx2(pmb, ks, ke, js, je+1, is, ie+1, b1, 0, 0,
                                             bx_N, bx_S);
        // Repeat calculation of x1 edge-centered wavespeeds as in HLL solver
        Real wli[(NWAVE)], wri[(NWAVE)];
        int ivx = IVX;
        int ivy = IVX + ((ivx-IVX)+1) % 3;
        int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=jl_buf; j<=ju_buf; ++j) {
            for (int i=is; i<=ie+1; ++i) {
              //--- Load L/R states into local variables
              // UCT with face-centered quantities
              wli[IDN]=wl_fc(IDN,k,j,i);
              wli[IVX]=wl_fc(ivx,k,j,i);
              wli[IVY]=wl_fc(ivy,k,j,i);
              wli[IVZ]=wl_fc(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wli[IPR]=wl_fc(IPR,k,j,i);
              wli[IBY]=wl_fc(IBY,k,j,i);
              wli[IBZ]=wl_fc(IBZ,k,j,i);

              wri[IDN]=wr_fc(IDN,k,j,i);
              wri[IVX]=wr_fc(ivx,k,j,i);
              wri[IVY]=wr_fc(ivy,k,j,i);
              wri[IVZ]=wr_fc(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wri[IPR]=wr_fc(IPR,k,j,i);
              wri[IBY]=wr_fc(IBY,k,j,i);
              wri[IBZ]=wr_fc(IBZ,k,j,i);
              Real bxi= b1_fc(k,j,i);

              Real cl=pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr=pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al=std::min((wri[IVX]-cr),(wli[IVX] - cl));
              Real ar=std::max((wli[IVX] + cl),(wri[IVX] + cr));
              Real bp=ar > 0.0 ? ar : 0.0;
              Real bm=al < 0.0 ? al : 0.0;
              pmb->pfield->alpha_plus_x1_(k,j,i) = bp;
              pmb->pfield->alpha_minus_x1_(k,j,i) = bm;
            }
          }
        }
        // Compute states for 3D UCT
        if (pmb->block_size.nx3 > 1) {
          pmb->precon->PiecewiseParabolicUCTx3(pmb, ks, ke+1, js, je, is, ie+1, wl, IVX,
                                               0, v_L3L1, v_R3L1);
          pmb->precon->PiecewiseParabolicUCTx3(pmb, ks, ke+1, js, je, is, ie+1, wl, IVZ,
                                               3, v_L3L1, v_R3L1);
          pmb->precon->PiecewiseParabolicUCTx3(pmb, ks, ke+1, js, je, is, ie+1, wr, IVX,
                                               0, v_L3R1, v_R3R1);
          pmb->precon->PiecewiseParabolicUCTx3(pmb, ks, ke+1, js, je, is, ie+1, wr, IVZ,
                                               3, v_L3R1, v_R3R1);

          // Limited transverse reconstructions: call PPMx3() for single-state b_x
          pmb->precon->PiecewiseParabolicUCTx3(pmb, ks, ke+1, js, je, is, ie+1, b1, 0, 0,
                                               bx_L3, bx_R3);
        } // end UCT if 3D
      } // end if 2D or 3D
    } else { // end if (order == 4) UCT4x1
      // compute weights for GS07 CT algorithm
      for (int k=kl_buf; k<=ku_buf; ++k) {
        for (int j=jl_buf; j<=ju_buf; ++j) {
          pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x1flux(IDN,k,j,i)
                / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
            Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
            w_x1f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
          }
        }
      }
    } // end GS07 CT algorithm
  } // end if (MAGNETIC_FIELDS_ENABLED)

//----------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {
    // set the loop limits
    il=is, iu=ie, kl=ks, ku=ke;

    if (order != 4) {
      if (MAGNETIC_FIELDS_ENABLED) {
        if (pmb->block_size.nx3 == 1) // 2D
          il=is-1, iu=ie+1, kl=ks, ku=ke;
        else // 3D
          il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
      }
    } else { // fourth-order corrections
      if (MAGNETIC_FIELDS_ENABLED) {
        if(pmb->block_size.nx3 == 1) { // 2D
          il=is-3, iu=ie+3, kl=ks, ku=ke;
        } else { // 3D
          il=is-3, iu=ie+3, kl=ks-3, ku=ke+3;
        }
      } else { // fourth-order hydro has same reqs as second-order MHD
        if (pmb->block_size.nx3 == 1) // 2D
          il=is-1, iu=ie+1, kl=ks, ku=ke;
        else // 3D
          il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
      }
    } // end order==4

    // Set transverse loop limits for quantities calculated to full accuracy
    il_buf=il, iu_buf=iu, kl_buf=kl, ku_buf=ku;
    if (order == 4) {
      if(pmb->block_size.nx3 == 1) // 2D
        il_buf+=1, iu_buf-=1;
      else // 3D
        il_buf+=1, iu_buf-=1, kl_buf+=1, ku_buf-=1;
    }

    // reconstruct L/R states at j
    if (order == 1) {
      pmb->precon->DonorCellX2(pmb, kl, ku, js, je+1, il, iu, w, bcc, wl, wr);
    } else if (order == 2) {
      pmb->precon->PiecewiseLinearX2(pmb, kl, ku, js, je+1, il, iu, w, bcc, wl, wr);
    } else {
      pmb->precon->PiecewiseParabolicX2(pmb, kl, ku, js, je+1, il, iu, w, bcc, wl, wr);
    }

    // compute fluxes, store directly into 3D arrays
    // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
    // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
    RiemannSolver(kl, ku, js, je+1, il, iu, IVY, b2, wl, wr, x2flux, e1x2, e3x2);

    // begin x2 fourth-order hydro and MHD
    //------------------------------------------------------------------------------
    if (order == 4) {
      if (MAGNETIC_FIELDS_ENABLED) {
        // Copy electric fields back into x2flux array
        for (int k=kl; k<=ku; ++k) {
          for (int j=js; j<=je+1; ++j) {
            for (int i=il; i<=iu; ++i) {
              x2flux(IBY,k,j,i) = e1x2(k,j,i); // technically -e1x2
              x2flux(IBZ,k,j,i) = e3x2(k,j,i);
            }
          }
        }
      }
      // Compute Laplacian of primitive Riemann states on x2 faces
      pmb->pcoord->LaplacianX2(wl, laplacian_l_fc, il, iu, js, je+1, kl, ku, 0, NWAVE-1);
      pmb->pcoord->LaplacianX2(wr, laplacian_r_fc, il, iu, js, je+1, kl, ku, 0, NWAVE-1);

      // TODO(kfelker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(kfelker): also, this may need to be dx1v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx2f(js);  // pco->dx2f(j); inside loop
      Real C = (h*h)/24.0;

      // Approximate x2 face-centered primitive Riemann states
      for (int n=0; n<NWAVE; ++n) {
        for (int k=kl; k<=ku; ++k) {
          for (int j=js; j<=je+1; ++j) {
            //pmb->pcoord->CenterWidth2(k, j, il, iu, dxw);
            for (int i=il; i<=iu; ++i) {
              wl_fc(n,k,j,i) = wl(n,k,j,i) - C*laplacian_l_fc(n,k,j,i);
              wr_fc(n,k,j,i) = wr(n,k,j,i) - C*laplacian_r_fc(n,k,j,i);
              // reapply primitive variable floors to face-centered L/R Riemann states
              // TODO(kfelker): only needs to be called 1x for all NWAVE
              pmb->peos->ApplyPrimitiveFloors(wl_fc, k, j, i);
              pmb->peos->ApplyPrimitiveFloors(wr_fc, k, j, i);
            }
          }
        }
      }

      // Compute x2 interface fluxes from face-centered primitive variables
      RiemannSolver(kl, ku, js, je+1, il, iu, IVY, b2_fc, wl_fc, wr_fc, flux_fc,
                    flux_fc_IBY, flux_fc_IBZ);

      // Compute Laplacian of second-order accurate face-averaged flux on x1 faces
      pmb->pcoord->LaplacianX2(x2flux, laplacian_l_fc, il, iu, js, je+1, kl, ku,
                               0, NWAVE-1);

      // Correct face-averaged fluxes (Guzik eq. 10)
      for(int n=0; n<NWAVE; n++) {
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=js; j<=je+1; ++j) {
            //pmb->pcoord->CenterWidth2(k, j, il_buf, iu_buf, dxw);
            for(int i=il_buf; i<=iu_buf; i++) {
              x2flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_l_fc(n,k,j,i);
            }
          }
        }
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        // Copy electric fields back into x2flux array
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=js; j<=je+1; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              e1x2(k,j,i) = x2flux(IBY,k,j,i); // technically -e1x2
              e3x2(k,j,i) = x2flux(IBZ,k,j,i);
            }
          }
        }
      }
    } // end if (order == 4)
    //------------------------------------------------------------------------------
    // end x2 fourth-order hydro and MHD

    if (MAGNETIC_FIELDS_ENABLED) {
      //-------- begin fourth-order upwind constrained transport (UCT4x2)
      if (order == 4) {
        // Limited transverse reconstructions: call PPMx1() for vx, vy both L/R states
        // wl_{i,j-1/2} is N side of interface --> discontinuous states in x1, L=E,  R=W
        pmb->precon->PiecewiseParabolicUCTx1(pmb, ks, ke, js, je+1, is, ie+1, wl, IVX, 0,
                                             vl_temp, vr_temp);
        pmb->precon->PiecewiseParabolicUCTx1(pmb, ks, ke, js, je+1, is, ie+1, wl, IVY, 1,
                                             vl_temp, vr_temp);
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
        // wr_{i,j-1/2}  is S side of interface --> discontinuous states in x1, L=E, R=W
        pmb->precon->PiecewiseParabolicUCTx1(pmb, ks, ke, js, je+1, is, ie+1, wr, IVX, 0,
                                             vl_temp, vr_temp);
        pmb->precon->PiecewiseParabolicUCTx1(pmb, ks, ke, js, je+1, is, ie+1, wr, IVY, 1,
                                             vl_temp, vr_temp);

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
        // Limited transverse reconstructions: call PPMx1() for single-state b_y
        pmb->precon->PiecewiseParabolicUCTx1(pmb, ks, ke, js, je+1, is, ie+1, b2, 0, 0,
                                             by_E, by_W);
        // Repeat calculation of x2 edge-centered wavespeeds as in HLL solver
        Real wli[(NWAVE)], wri[(NWAVE)];
        int ivx = IVY;
        int ivy = IVX + ((ivx-IVX)+1) % 3;
        int ivz = IVX + ((ivx-IVX)+2) % 3;
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=js; j<=je+1; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              //--- Load L/R states into local variables
              // UCT with face-centered quantities
              wli[IDN]=wl_fc(IDN,k,j,i);
              wli[IVX]=wl_fc(ivx,k,j,i);
              wli[IVY]=wl_fc(ivy,k,j,i);
              wli[IVZ]=wl_fc(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wli[IPR]=wl_fc(IPR,k,j,i);
              wli[IBY]=wl_fc(IBY,k,j,i);
              wli[IBZ]=wl_fc(IBZ,k,j,i);

              wri[IDN]=wr_fc(IDN,k,j,i);
              wri[IVX]=wr_fc(ivx,k,j,i);
              wri[IVY]=wr_fc(ivy,k,j,i);
              wri[IVZ]=wr_fc(ivz,k,j,i);
              if (NON_BAROTROPIC_EOS) wri[IPR]=wr_fc(IPR,k,j,i);
              wri[IBY]=wr_fc(IBY,k,j,i);
              wri[IBZ]=wr_fc(IBZ,k,j,i);
              Real bxi= b2_fc(k,j,i);

              Real cl=pmb->peos->FastMagnetosonicSpeed(wli,bxi);
              Real cr=pmb->peos->FastMagnetosonicSpeed(wri,bxi);

              // eq 55 in Londrillo and Del Zanna 2004
              Real al=std::min((wri[IVX]-cr),(wli[IVX] - cl));
              Real ar=std::max((wli[IVX] + cl),(wri[IVX] + cr));
              Real bp=ar > 0.0 ? ar : 0.0;
              Real bm=al < 0.0 ? al : 0.0;
              pmb->pfield->alpha_plus_x2_(k,j,i) = bp;
              pmb->pfield->alpha_minus_x2_(k,j,i) = bm;
            }
          }
        }
      } else { // end if (order == 4) UCT4x1
        // compute weights for GS07 CT algorithm
        for (int k=kl_buf; k<=ku_buf; ++k) {
          for (int j=js; j<=je+1; ++j) {
            pmb->pcoord->CenterWidth2(k,j,il_buf,iu_buf,dxw);
#pragma omp simd
            for (int i=il_buf; i<=iu_buf; ++i) {
              Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*x2flux(IDN,k,j,i)
                  / (dxw(i)*(wl(IDN,k,j,i) + wr(IDN,k,j,i)));
              Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
              w_x2f(k,j,i) = 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
            }
          }
        }
      } // end GS07 CT algorithm
    } // end if (MAGNETIC_FIELDS_ENABLED)
  }

//----------------------------------------------------------------------------------------
// k-direction

  if (pmb->block_size.nx3 > 1) {

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

    // reconstruct L/R states at k
    if (order == 1) {
      pmb->precon->DonorCellX3(pmb, ks, ke+1, jl, ju, il, iu, w, bcc, wl, wr);
    } else if (order == 2) {
      pmb->precon->PiecewiseLinearX3(pmb, ks, ke+1, jl, ju, il, iu, w, bcc, wl, wr);
    } else {
      pmb->precon->PiecewiseParabolicX3(pmb, ks, ke+1, jl, ju, il, iu, w, bcc, wl, wr);
    }

    // compute fluxes, store directly into 3D arrays
    // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
    // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
    RiemannSolver(ks, ke+1, jl, ju, il, iu, IVZ, b3, wl, wr, x3flux, e2x3, e1x3);

    // begin x3 fourth-order hydro and MHD
    //------------------------------------------------------------------------------
    if (order == 4) {
      if (MAGNETIC_FIELDS_ENABLED) {
        // Copy electric fields back into x3flux array
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=jl; j<=ju; ++j) {
            for (int i=il; i<=iu; ++i) {
              x3flux(IBY,k,j,i) = e2x3(k,j,i); // technically -e2x3
              x3flux(IBZ,k,j,i) = e1x3(k,j,i);
            }
          }
        }
      }
      // Compute Laplacian of primitive Riemann states on x3 faces
      pmb->pcoord->LaplacianX3(wl, laplacian_l_fc, il, iu, jl, ju, ks, ke+1, 0, NWAVE-1);
      pmb->pcoord->LaplacianX3(wr, laplacian_r_fc, il, iu, jl, ju, ks, ke+1, 0, NWAVE-1);

      // TODO(kfelker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(kfelker): also, this may need to be dx1v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx3f(ks);  // pco->dx3f(k); inside loop
      Real C = (h*h)/24.0;

      // Approximate x3 face-centered primitive Riemann states
      for (int n=0; n<NWAVE; ++n) {
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=jl; j<=ju; ++j) {
            //pmb->pcoord->CenterWidth3(k, j, il, iu, dxw);
            for (int i=il; i<=iu; ++i) {
              wl_fc(n,k,j,i) = wl(n,k,j,i) - C*laplacian_l_fc(n,k,j,i);
              wr_fc(n,k,j,i) = wr(n,k,j,i) - C*laplacian_r_fc(n,k,j,i);
              // reapply primitive variable floors to face-centered L/R Riemann states
              // TODO(kfelker): only needs to be called 1x for all NWAVE
              pmb->peos->ApplyPrimitiveFloors(wl_fc, k, j, i);
              pmb->peos->ApplyPrimitiveFloors(wr_fc, k, j, i);
            }
          }
        }
      }

      // Compute x3 interface fluxes from face-centered primitive variables
      RiemannSolver(ks, ke+1, jl, ju, il, iu, IVZ, b3_fc, wl_fc, wr_fc, flux_fc,
                    flux_fc_IBY, flux_fc_IBZ);

      // Compute Laplacian of second-order accurate face-averaged flux on x1 faces
      pmb->pcoord->LaplacianX3(x3flux, laplacian_l_fc, il, iu, jl, ju, ks, ke+1,
                               0, NWAVE-1);

      // Correct face-averaged fluxes (Guzik eq. 10)
      for(int n=0; n<NWAVE; n++) {
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=jl_buf; j<=ju_buf; ++j) {
            //pmb->pcoord->CenterWidth3(k, j, il_buf, iu_buf, dxw);
            for(int i=il_buf; i<=iu_buf; i++) {
                x3flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_l_fc(n,k,j,i);
            }
          }
        }
      }

      if (MAGNETIC_FIELDS_ENABLED) {
        // Copy high-order magnetic x1flux back into electric field arrays
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=jl_buf; j<=ju_buf; ++j) {
            for (int i=il_buf; i<=iu_buf; ++i) {
              e2x3(k,j,i) = x3flux(IBY,k,j,i); // technically -e2x3
              e1x3(k,j,i) = x3flux(IBZ,k,j,i);
            }
          }
        }
      }
    } // end if (order == 4)
    //------------------------------------------------------------------------------
    // end x3 fourth-order hydro and MHD

    // compute weights for GS07 CT algorithm
    if (MAGNETIC_FIELDS_ENABLED) {
      for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl_buf; j<=ju_buf; ++j) {
        pmb->pcoord->CenterWidth3(k,j,il_buf,iu_buf,dxw);
#pragma omp simd
        for (int i=il_buf; i<=iu_buf; ++i) {
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
