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
#include "hydro_diffusion/hydro_diffusion.hpp"
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


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, const int order) {
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
  // TODO(kfelker): fix loop limits for fourth-order hydro
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
        pmb->precon->DonorCellX1(pmb, k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(pmb, k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX1(pmb, k, j, is-1, ie+1, w, bcc, wl_, wr_);
      }

      // compute fluxes, store directly into 3D arrays
      // x1flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
      // x1flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
      pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw_);
      RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_);

    }
  }
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro


//----------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {

    // set the loop limits
    il=is-1, iu=ie+1, kl=ks, ku=ke;
    // TODO(kfelker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->block_size.nx3 == 1) // 2D
      kl=ks, ku=ke;
    else // 3D
      kl=ks-1, ku=ke+1;
    //    }

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(pmb, k, js-1, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(pmb, k, js-1, il, iu, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX2(pmb, k, js-1, il, iu, w, bcc, wl_, wr_);
      }
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(pmb, k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(pmb, k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          pmb->precon->PiecewiseParabolicX2(pmb, k, j, il, iu, w, bcc, wlb_, wr_);
        }

        // compute fluxes, store directly into 3D arrays
        // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
        // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
        pmb->pcoord->CenterWidth2(k,j,il,iu,dxw_);
        RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, x2flux, e1x2, e3x2, w_x2f, dxw_);

        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
      }
    }
  }

//----------------------------------------------------------------------------------------
// k-direction

  if (pmb->block_size.nx3 > 1) {

    // set the loop limits
    // TODO(kfelker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED)
    il=is-1, iu=ie+1, jl=js-1, ju=je+1;

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(pmb, ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(pmb, ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else {
        pmb->precon->PiecewiseParabolicX3(pmb, ks-1, j, il, iu, w, bcc, wl_, wr_);
      }
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(pmb, k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(pmb, k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          pmb->precon->PiecewiseParabolicX3(pmb, k, j, il, iu, w, bcc, wlb_, wr_);
        }

        // compute fluxes, store directly into 3D arrays
        // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
        // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
        pmb->pcoord->CenterWidth3(k,j,il,iu,dxw_);
        RiemannSolver(k, j, il, iu, IVZ, b3, wl_, wr_, x3flux, e2x3, e1x3, w_x3f, dxw_);

        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
      }
    }
  }


  if (SELF_GRAVITY_ENABLED) AddGravityFlux(); // add gravity flux directly

  if (!STS_ENABLED) { // add diffusion fluxes
    if (phdif->hydro_diffusion_defined) {
      if (phdif->nu_iso > 0.0 || phdif->nu_aniso > 0.0)
        phdif->AddHydroDiffusionFlux(phdif->visflx,flux);

      if (NON_BAROTROPIC_EOS) {
        if (phdif->kappa_iso > 0.0 || phdif->kappa_aniso > 0.0)
          phdif->AddHydroDiffusionEnergyFlux(phdif->cndflx,flux);
      }
    }

    if (MAGNETIC_FIELDS_ENABLED && NON_BAROTROPIC_EOS) {
      if (pmb->pfield->pfdif->field_diffusion_defined)
        pmb->pfield->pfdif->AddPoyntingFlux(pmb->pfield->pfdif->pflux);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes_STS
//  \brief Calculate Hydrodynamic Diffusion Fluxes for STS

void Hydro::CalculateFluxes_STS() {
  MeshBlock *pmb=pmy_block;
  // add diffusion fluxes
  if (phdif->hydro_diffusion_defined) {
    if (phdif->nu_iso > 0.0 || phdif->nu_aniso > 0.0)
      phdif->AddHydroDiffusionFlux(phdif->visflx,flux);

    if (NON_BAROTROPIC_EOS) {
      if (phdif->kappa_iso > 0.0 || phdif->kappa_aniso > 0.0)
        phdif->AddHydroDiffusionEnergyFlux(phdif->cndflx,flux);
    }
  }

  if (MAGNETIC_FIELDS_ENABLED && NON_BAROTROPIC_EOS) {
    if (pmb->pfield->pfdif->field_diffusion_defined)
      pmb->pfield->pfdif->AddPoyntingFlux(pmb->pfield->pfdif->pflux);
  }

  return;
}
