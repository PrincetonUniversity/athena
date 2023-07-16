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
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file cr_transport.cpp
//  \brief implementation of cosmic ray transport
//======================================================================================

// C headers

// C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../reconstruct/reconstruction.hpp"
#include "../../utils/utils.hpp"
#include "../cr.hpp"

// class header
#include "./cr_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

static constexpr Real tau_asymptotic_lim=1.0e-3;

void CRIntegrator::CalculateFluxes(
    AthenaArray<Real> &w, AthenaArray<Real> &bcc,
    AthenaArray<Real> &cr, const int order) {
  CosmicRay *pcr=pmy_cr;
  MeshBlock *pmb=pcr->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Real invlim = 1.0/pcr->vmax;

  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2,
      ncells3 = pmb->ncells3;

  AthenaArray<Real> &x1flux=pcr->flux[X1DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if (ncells2 > 1) {
    if (ncells3 == 1) {
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    } else {
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }

  //---------------------------------------------------------------------
  for (int k=0; k<ncells3; ++k) {
    for (int j=0; j<ncells2; ++j) {
      // diffusion velocity along the direction of sigma vector
      // We first assume B is along x coordinate
      // Then rotate according to B direction to the actual acooridnate
      for (int i=0; i<ncells1; ++i) {
        Real eddxx=1.0/3.0;
        Real totsigma = pcr->sigma_diff(0,k,j,i);
        if (pcr->stream_flag)
          totsigma = 1.0/(1.0/pcr->sigma_diff(0,k,j,i)
                          + 1.0/pcr->sigma_adv(0,k,j,i));
        Real taux = taufact_ * totsigma * pco->dx1f(i);
        taux = taux * taux/(2.0 * eddxx);
        Real diffv = 1.0;

        if (taux < tau_asymptotic_lim) {
          diffv = std::sqrt((1.0 - 0.5* taux));
        } else {
          diffv = std::sqrt((1.0 - std::exp(-taux)) / taux);
        }
        pcr->v_diff(0,k,j,i) = pcr->vmax * std::sqrt(eddxx) * diffv;
      }

      // y direction
      if (ncells2 >1) {
        pco->CenterWidth2(k,j,0,ncells1-1,cwidth2_);
        // get the optical depth across the cell
        for (int i=0; i<ncells1; ++i) {
          Real eddyy=1.0/3.0;
          Real totsigma = pcr->sigma_diff(1,k,j,i);
          if (pcr->stream_flag)
            totsigma = 1.0/(1.0/pcr->sigma_diff(1,k,j,i)
                            + 1.0/pcr->sigma_adv(1,k,j,i));
          Real tauy = taufact_ * totsigma * cwidth2_(i);
          tauy = tauy * tauy/(2.0 * eddyy);

          Real diffv = 1.0;

          if (tauy < tau_asymptotic_lim) {
            diffv = std::sqrt((1.0 - 0.5* tauy));
          } else {
            diffv = std::sqrt((1.0 - std::exp(-tauy)) / tauy);
          }
          pcr->v_diff(1,k,j,i) = pcr->vmax * std::sqrt(eddyy) * diffv;
        }
      } else {
        for (int i=0; i<ncells1; ++i)
          pcr->v_diff(1,k,j,i) = 0.0;
      }

      // z direction
      if (ncells3 > 1) {
        pco->CenterWidth3(k,j,0,ncells1-1,cwidth3_);
        // get the optical depth across the cell
        for (int i=0; i<ncells1; ++i) {
          Real eddzz=1.0/3.0;
          Real totsigma = pcr->sigma_diff(2,k,j,i);
          if (pcr->stream_flag)
            totsigma = 1.0/(1.0/pcr->sigma_diff(2,k,j,i)
                            + 1.0/pcr->sigma_adv(2,k,j,i));
          Real tauz = taufact_ * totsigma * cwidth3_(i);
          tauz = tauz * tauz/(2.0 * eddzz);

          Real diffv = 1.0;

          if (tauz < tau_asymptotic_lim) {
            diffv = std::sqrt((1.0 - 0.5* tauz));
          } else {
            diffv = std::sqrt((1.0 - std::exp(-tauz)) / tauz);
          }
          pcr->v_diff(2,k,j,i) = pcr->vmax * std::sqrt(eddzz) * diffv;
        }
      } else {
        for (int i=0; i<ncells1; ++i)
          pcr->v_diff(2,k,j,i) = 0.0;
      }

      // rotate the v_diff vector to the local coordinate
      if (MAGNETIC_FIELDS_ENABLED) {
        for (int i=0; i<ncells1; ++i) {
          InvRotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                       pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),
                       pcr->v_diff(0,k,j,i),pcr->v_diff(1,k,j,i),
                       pcr->v_diff(2,k,j,i));
          // take the absolute value
          // Also add the Alfven velocity for the streaming flux
          pcr->v_diff(0,k,j,i) = std::abs(pcr->v_diff(0,k,j,i));

          pcr->v_diff(1,k,j,i) = std::abs(pcr->v_diff(1,k,j,i));

          pcr->v_diff(2,k,j,i) = std::abs(pcr->v_diff(2,k,j,i));
        }
      }
      // need to add additional sound speed for stability
      for (int i=0; i<ncells1; ++i) {
        Real cr_sound_x = vel_flx_flag_ * std::sqrt((4.0/9.0)
                                                    * cr(CRE,k,j,i)/w(IDN,k,j,i));
        pcr->v_diff(0,k,j,i) += cr_sound_x;
        if (ncells2 > 1) {
          Real cr_sound_y = vel_flx_flag_ * std::sqrt((4.0/9.0)
                                                      * cr(CRE,k,j,i)/w(IDN,k,j,i));
          pcr->v_diff(1,k,j,i) += cr_sound_y;
        }

        if (ncells3 > 1) {
          Real cr_sound_z = vel_flx_flag_ * std::sqrt((4.0/9.0)
                                                      * cr(CRE,k,j,i)/w(IDN,k,j,i));

          pcr->v_diff(2,k,j,i) += cr_sound_z;
        }
      }
    }
  }
  // prepare Array for reconstruction
  for (int n=0; n<NCR; ++n) {
    for (int k=0; k<ncells3; ++k) {
      for (int j=0; j<ncells2; ++j) {
        for (int i=0; i<ncells1; ++i) {
          ucr_vel_(n,k,j,i) = cr(n,k,j,i);
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------
  // i-direction

  // add vx velocity
  for (int k=0; k<ncells3; ++k) {
    for (int j=0; j<ncells2; ++j) {
      for (int i=0; i<ncells1; ++i) {
        ucr_vel_(NCR,k,j,i) = w(IVX,k,j,i);
      }
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // First, need to do reconstruction
      // to reconstruct Ec, Fc, vel, v_a and
      // return Ec,Fc and signal speed at left and right state
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, ucr_vel_, ucr_l_, ucr_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, ucr_vel_, ucr_l_, ucr_r_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, ucr_vel_, ucr_l_, ucr_r_);
      }

      // get the optical depth across the cell
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        vdiff_l_(i) = pcr->v_diff(0,k,j,i-1);
        vdiff_r_(i) = pcr->v_diff(0,k,j,i);
      }

      // calculate the flux
      CRFlux(CRF1, is, ie+1, ucr_l_, ucr_r_, vdiff_l_, vdiff_r_, dflx_);
      // store the flux
      for (int n=0; n<NCR; ++n) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          x1flux(n,k,j,i) = dflx_(n,i);
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------
  // j-direction
  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
    // add vy velocity
    for (int k=0; k<ncells3; ++k) {
      for (int j=0; j<ncells2; ++j) {
        for (int i=0; i<ncells1; ++i) {
          ucr_vel_(NCR,k,j,i) = w(IVY,k,j,i);
        }
      }
    }

    il=is-1; iu=ie+1; kl=ks; ku=ke;
    if (ncells3 ==  1) // 2D
      kl = ks, ku = ke;
    else // 3D
      kl = ks-1, ku = ke+1;
    for (int k=kl; k<=ku; ++k) {
      //reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      }


      for (int j=js; j<=je+1; ++j) {
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        }

        // get the optical depth across the cell
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          vdiff_l_(i) = pcr->v_diff(1,k,j-1,i);
          vdiff_r_(i) = pcr->v_diff(1,k,j,i);
        }
        // calculate the flux
        CRFlux(CRF2, il, iu, ucr_l_, ucr_r_, vdiff_l_, vdiff_r_, dflx_);
        // store the flux
        for (int n=0; n<NCR; ++n) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            x2flux(n,k,j,i) = dflx_(n,i);
          }
        }
        // swap the array for next cycle
        ucr_l_.SwapAthenaArray(ucr_lb_);
      }
    }
  }
  //  k-direction
  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
    il =is-1, iu=ie+1, jl=js-1, ju=je+1;
    // add vz velocity
    for (int k=0; k<ncells3; ++k) {
      for (int j=0; j<ncells2; ++j) {
        for (int i=0; i<ncells1; ++i) {
          ucr_vel_(NCR,k,j,i) = w(IVZ,k,j,i);
        }// end i
      }// end j
    }// end k

    for (int j=jl; j<=ju; ++j) {
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      }

      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        }

#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          vdiff_l_(i) = pcr->v_diff(2,k-1,j,i);
          vdiff_r_(i) = pcr->v_diff(2,k,j,i);
        }
        // calculate the flux
        CRFlux(CRF3, il, iu, ucr_l_, ucr_r_, vdiff_l_, vdiff_r_, dflx_);
        for (int n=0; n<NCR; ++n) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            x3flux(n,k,j,i) = dflx_(n,i);
          }
        }
        // swap the array for next cycle
        ucr_l_.SwapAthenaArray(ucr_lb_);
      }
    }
  }

  //---------------------------------------------------------------------------------------
  // Now calculate Grad Pc and the associated heating term
  // the flux divergence term is Grad P_c for \partial F_c/\partial t
  // only do this for the MHD case and along direction perpendicular
  // to the magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,x1face_area_);
        pmb->pcoord->CellVolume(k,j,is,ie,cell_volume_);
        // x1 direction
        for (int n=0; n<3; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            grad_pc_(n,k,j,i) = (x1face_area_(i+1)*x1flux(CRF1+n,k,j,i+1)
                                 - x1face_area_(i)*x1flux(CRF1+n,k,j,i))/cell_volume_(i);
          }
        }

        if (pmb->block_size.nx2 > 1) {
          AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
          pmb->pcoord->Face2Area(k,j  ,is,ie,x2face_area_   );
          pmb->pcoord->Face2Area(k,j+1,is,ie,x2face_area_p1_);
          for (int n=0; n<3; ++n) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              grad_pc_(n,k,j,i) += (x2face_area_p1_(i)*x2flux(CRF1+n,k,j+1,i)
                                    - x2face_area_(i)*x2flux(CRF1+n,k,j,i))
                                   /cell_volume_(i);
            }
          }
        }
        if (pmb->block_size.nx3 > 1) {
          AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
          pmb->pcoord->Face3Area(k  ,j,is,ie,x3face_area_);
          pmb->pcoord->Face3Area(k+1,j,is,ie,x3face_area_p1_);
          for (int n=0; n<3; ++n) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              grad_pc_(n,k,j,i) += (x3face_area_p1_(i) *x3flux(CRF1+n,k+1,j,i)
                                    - x3face_area_(i)*x3flux(CRF1+n,k,j,i))
                                   /cell_volume_(i);
            }
          }
        }
        for (int n=0; n<3; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            grad_pc_(n,k,j,i) *= invlim;
          }
        }

        // need to subtract the coordinate source term to get the actual grad Pc for c
        // curlinear coordinate system
        pmb->pcoord->CRGradPcCoordTermsDivergence(cr, grad_pc_);

        // update streaming with grad_pc

        pcr->UpdateStreaming(pmb, cr,w,bcc,grad_pc_,k,j,is,ie);

        // calculate streaming velocity with magnetic field
        for (int i=is; i<=ie; ++i) {
          Real v1 = w(IVX,k,j,i);
          Real v2 = w(IVY,k,j,i);
          Real v3 = w(IVZ,k,j,i);

          Real dpcdx = grad_pc_(0,k,j,i);
          Real dpcdy = grad_pc_(1,k,j,i);
          Real dpcdz = grad_pc_(2,k,j,i);

          RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                    pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),dpcdx,dpcdy,dpcdz);

          RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                    pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),v1,v2,v3);

          // only calculate v_dot_gradpc perpendicular to B
          // perpendicular direction only has flow velocity, no streaming velocity
          Real v_dot_gradpc = v2 * dpcdy + v3 * dpcdz;

          ec_source_(k,j,i) = v_dot_gradpc;
        }
      }
    }
  }
  //-----------------------------------------------------------------------
  // calculate coordinate source terms for Cosmic ray
  pco->AddCRCoordTermsDivergence(cr,coord_source_);
}


void CRIntegrator::FluxDivergence(const Real wght, AthenaArray<Real> &cr_out) {
  CosmicRay *pcr=pmy_cr;
  MeshBlock *pmb = pcr->pmy_block;
  //  Coordinates *pco = pmb->pcoord;

  AthenaArray<Real> &x1flux=pcr->flux[X1DIR];
  AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
  AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);

      for (int n=0; n<NCR; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }// end n
      }// End i

      // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for (int n=0; n<NCR; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }// end nx2

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);

        for (int n=0; n<NCR; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }// end nx3
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for (int n=0; n<NCR; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          cr_out(n,k,j,i) -= wght*dflx(n,i)/vol(i);
        }
      }
    }
  }

  // Add coordinate source term
  for (int n=0; n<NCR; ++n)
    for (int k=ks; k<=ke; ++k)
      for (int j=js; j<=je; ++j)
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          cr_out(n,k,j,i) += wght * coord_source_(n,k,j,i);

  // check Ec is positive
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        if (cr_out(CRE,k,j,i) < TINY_NUMBER)
          cr_out(CRE,k,j,i) = TINY_NUMBER;
      }
}
