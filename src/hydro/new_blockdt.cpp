//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_blockdt.cpp
//  \brief computes timestep using CFL condition on a MEshBlock

// C/C++ headers
#include <algorithm>  // max(), min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
// \!fn Real Hydro::NewBlockTimeStep(void)
// \brief calculate the minimum timestep within a MeshBlock

Real Hydro::NewBlockTimeStep(void) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> w,bcc,b_x1f,b_x2f,b_x3f;
  w.InitWithShallowCopy(pmb->phydro->w);
  if (MAGNETIC_FIELDS_ENABLED) {
    bcc.InitWithShallowCopy(pmb->pfield->bcc);
    b_x1f.InitWithShallowCopy(pmb->pfield->b.x1f);
    b_x2f.InitWithShallowCopy(pmb->pfield->b.x2f);
    b_x3f.InitWithShallowCopy(pmb->pfield->b.x3f);
  }

  AthenaArray<Real> dt1, dt2, dt3;
  dt1.InitWithShallowCopy(dt1_);
  dt2.InitWithShallowCopy(dt2_);
  dt3.InitWithShallowCopy(dt3_);
  Real wi[(NWAVE)];

  Real min_dt = (FLT_MAX);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->CenterWidth1(k,j,is,ie,dt1);
      pmb->pcoord->CenterWidth2(k,j,is,ie,dt2);
      pmb->pcoord->CenterWidth3(k,j,is,ie,dt3);

      // Newtonian dynamics
      if (!RELATIVISTIC_DYNAMICS) {
#pragma ivdep
        for (int i=is; i<=ie; ++i) {
          wi[IDN]=w(IDN,k,j,i);
          wi[IVX]=w(IVX,k,j,i);
          wi[IVY]=w(IVY,k,j,i);
          wi[IVZ]=w(IVZ,k,j,i);
          if (NON_BAROTROPIC_EOS) wi[IPR]=w(IPR,k,j,i);

          // Newtonian MHD
          if (MAGNETIC_FIELDS_ENABLED) {

            Real bx = bcc(IB1,k,j,i) + fabs(b_x1f(k,j,i)-bcc(IB1,k,j,i));
            wi[IBY] = bcc(IB2,k,j,i);
            wi[IBZ] = bcc(IB3,k,j,i);
            Real cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
            dt1(i) /= (fabs(wi[IVX]) + cf);

            wi[IBY] = bcc(IB3,k,j,i);
            wi[IBZ] = bcc(IB1,k,j,i);
            bx = bcc(IB2,k,j,i) + fabs(b_x2f(k,j,i)-bcc(IB2,k,j,i));
            cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
            dt2(i) /= (fabs(wi[IVY]) + cf);

            wi[IBY] = bcc(IB1,k,j,i);
            wi[IBZ] = bcc(IB2,k,j,i);
            bx = bcc(IB3,k,j,i) + fabs(b_x3f(k,j,i)-bcc(IB3,k,j,i));
            cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
            dt3(i) /= (fabs(wi[IVZ]) + cf);

          // Newtonian hydro
          } else {

            Real cs = pmb->peos->SoundSpeed(wi);
            dt1(i) /= (fabs(wi[IVX]) + cs);
            dt2(i) /= (fabs(wi[IVY]) + cs);
            dt3(i) /= (fabs(wi[IVZ]) + cs);

          }
        }

      // GR
      } else if (GENERAL_RELATIVITY) {

        // Extract adiabatic index and prepare metric
        Real gamma_adi = pmy_block->peos->GetGamma();
        Real gamma_prime = gamma_adi / (gamma_adi - 1.0);
        pmy_block->pcoord->CellMetric(k, j, is, ie, g_, gi_);

        // GR MHD
        if (MAGNETIC_FIELDS_ENABLED) {
          for (int i = is; i <= ie; ++i) {

            // Extract and calculate thermodynamic state
            Real rho = w(IDN,k,j,i);
            Real pgas = w(IPR,k,j,i);
            Real wgas = rho + gamma_prime * pgas;

            // Extract and calculate velocity
            Real alpha = std::sqrt(-1.0 / gi_(I00,i));
            Real uu1 = w(IVX,k,j,i);
            Real uu2 = w(IVY,k,j,i);
            Real uu3 = w(IVZ,k,j,i);
            Real uu_sq = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real gamma_rel = std::sqrt(1.0 + uu_sq);
            Real u0 = gamma_rel / alpha;
            Real u1 = uu1 - alpha * gamma_rel * gi_(I01,i);
            Real u2 = uu2 - alpha * gamma_rel * gi_(I02,i);
            Real u3 = uu3 - alpha * gamma_rel * gi_(I03,i);

            // Extract and calculate magnetic field
            Real bb1 = bcc(IB1,k,j,i);
            Real bb2 = bcc(IB2,k,j,i);
            Real bb3 = bcc(IB3,k,j,i);
            Real b0 = g_(I01,i) * u0 * bb1 + g_(I02,i) * u0 * bb2 + g_(I03,i) * u0 * bb3
                + g_(I11,i) * u1 * bb1 + g_(I12,i) * (u1 * bb2 + u2 * bb1)
                + g_(I13,i) * (u1 * bb3 + u3 * bb1) + g_(I22,i) * u2 * bb2
                + g_(I23,i) * (u2 * bb3 + u3 * bb2) + g_(I33,i) * u3 * bb3;
            Real b1 = (bb1 + b0 * u1) / u0;
            Real b2 = (bb2 + b0 * u2) / u0;
            Real b3 = (bb3 + b0 * u3) / u0;
            Real b_sq = g_(I00,i) * SQR(b0) + 2.0 * g_(I01,i) * b0 * b1
                + 2.0 * g_(I02,i) * b0 * b2 + 2.0 * g_(I03,i) * b0 * b3
                + g_(I11,i) * SQR(b1) + 2.0 * g_(I12,i) * b1 * b2
                + 2.0 * g_(I13,i) * b1 * b3 + g_(I22,i) * SQR(b2)
                + 2.0 * g_(I23,i) * b2 * b3 + g_(I33,i) * SQR(b3);

            // Calculate signal crossing times
            Real lambda_p, lambda_m;
            pmy_block->peos->FastMagnetosonicSpeedsGR(wgas, pgas, u0, u1, b_sq,
                gi_(I00,i), gi_(I01,i), gi_(I11,i), &lambda_p, &lambda_m);
            dt1(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
            pmy_block->peos->FastMagnetosonicSpeedsGR(wgas, pgas, u0, u2, b_sq,
                gi_(I00,i), gi_(I02,i), gi_(I22,i), &lambda_p, &lambda_m);
            dt2(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
            pmy_block->peos->FastMagnetosonicSpeedsGR(wgas, pgas, u0, u3, b_sq,
                gi_(I00,i), gi_(I03,i), gi_(I33,i), &lambda_p, &lambda_m);
            dt3(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
          }

        // GR hydro
        } else {
          for (int i = is; i <= ie; ++i) {

            // Extract and calculate thermodynamic state
            Real rho = w(IDN,k,j,i);
            Real pgas = w(IPR,k,j,i);
            Real wgas = rho + gamma_prime * pgas;

            // Extract and calculate velocity
            Real alpha = std::sqrt(-1.0 / gi_(I00,i));
            Real uu1 = w(IVX,k,j,i);
            Real uu2 = w(IVY,k,j,i);
            Real uu3 = w(IVZ,k,j,i);
            Real uu_sq = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
                + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22,i) * SQR(uu2)
                + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
            Real gamma_rel = std::sqrt(1.0 + uu_sq);
            Real u0 = gamma_rel / alpha;
            Real u1 = uu1 - alpha * gamma_rel * gi_(I01,i);
            Real u2 = uu2 - alpha * gamma_rel * gi_(I02,i);
            Real u3 = uu3 - alpha * gamma_rel * gi_(I03,i);

            // Calculate signal crossing times
            Real lambda_p, lambda_m;
            pmy_block->peos->SoundSpeedsGR(wgas, pgas, u0, u1, gi_(I00,i), gi_(I01,i),
                gi_(I11,i), &lambda_p, &lambda_m);
            dt1(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
            pmy_block->peos->SoundSpeedsGR(wgas, pgas, u0, u2, gi_(I00,i), gi_(I02,i),
                gi_(I22,i), &lambda_p, &lambda_m);
            dt2(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
            pmy_block->peos->SoundSpeedsGR(wgas, pgas, u0, u3, gi_(I00,i), gi_(I03,i),
                gi_(I33,i), &lambda_p, &lambda_m);
            dt3(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
          }
        }

      // SR
      } else {

        // SR MHD
        if (MAGNETIC_FIELDS_ENABLED) {

          // Place primitives in new array
          for (int n = 0; n < NHYDRO; ++n) {
            for (int i = is; i <= ie; ++i) {
              prim_field_(n,0,0,i) = w(n,k,j,i);
            }
          }

          // Calculate x1-crossing times
          for (int i = is; i <= ie; ++i) {
            bb_normal_(i) = bcc(IB1,k,j,i);
            prim_field_(IBY,0,0,i) = bcc(IB2,k,j,i);
            prim_field_(IBY,0,0,i) = bcc(IB3,k,j,i);
          }
          pmy_block->peos->FastMagnetosonicSpeedsSR(prim_field_, bb_normal_, k, j, is, ie,
              IVX, lambdas_p_, lambdas_m_);
          for (int i = is; i <= ie; ++i) {
            dt1(i) /= std::max(std::fabs(lambdas_p_(i)), std::fabs(lambdas_m_(i)));
          }

          // Calculate x2-crossing times
          for (int i = is; i <= ie; ++i) {
            bb_normal_(i) = bcc(IB2,k,j,i);
            prim_field_(IBY,i) = bcc(IB3,k,j,i);
            prim_field_(IBY,i) = bcc(IB1,k,j,i);
          }
          pmy_block->peos->FastMagnetosonicSpeedsSR(prim_field_, bb_normal_, k, j, is, ie,
              IVY, lambdas_p_, lambdas_m_);
          for (int i = is; i <= ie; ++i) {
            dt2(i) /= std::max(std::fabs(lambdas_p_(i)), std::fabs(lambdas_m_(i)));
          }

          // Calculate x3-crossing times
          for (int i = is; i <= ie; ++i) {
            bb_normal_(i) = bcc(IB3,k,j,i);
            prim_field_(IBY,i) = bcc(IB1,k,j,i);
            prim_field_(IBY,i) = bcc(IB2,k,j,i);
          }
          pmy_block->peos->FastMagnetosonicSpeedsSR(prim_field_, bb_normal_, k, j, is, ie,
              IVZ, lambdas_p_, lambdas_m_);
          for (int i = is; i <= ie; ++i) {
            dt3(i) /= std::max(std::fabs(lambdas_p_(i)), std::fabs(lambdas_m_(i)));
          }

        // SR hydro
        } else {

          // Extract adiabatic index
          Real gamma_adi = pmy_block->peos->GetGamma();
          Real gamma_prime = gamma_adi / (gamma_adi - 1.0);

          // Go through cells in x1-direction
          for (int i = is; i <= ie; ++i) {

            // Extract and calculate thermodynamic state
            Real rho = w(IDN,k,j,i);
            Real pgas = w(IPR,k,j,i);
            Real wgas = rho + gamma_prime * pgas;

            // Extract and calculate velocity
            Real v1 = w(IVX,k,j,i);
            Real v2 = w(IVY,k,j,i);
            Real v3 = w(IVZ,k,j,i);
            Real gamma_rel_sq = 1.0 / (1.0 - SQR(v1) - SQR(v2) - SQR(v3));

            // Calculate signal crossing times
            Real lambda_p, lambda_m;
            pmy_block->peos->SoundSpeedsSR(wgas, pgas, v1, gamma_rel_sq, &lambda_p,
                &lambda_m);
            dt1(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
            pmy_block->peos->SoundSpeedsSR(wgas, pgas, v2, gamma_rel_sq, &lambda_p,
                &lambda_m);
            dt2(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
            pmy_block->peos->SoundSpeedsSR(wgas, pgas, v3, gamma_rel_sq, &lambda_p,
                &lambda_m);
            dt3(i) /= std::max(std::fabs(lambda_p), std::fabs(lambda_m));
          }
        }
      }

      // compute minimum of (v1 +/- C)
      for (int i=is; i<=ie; ++i) {
        Real& dt_1 = dt1(i);
        min_dt = std::min(min_dt,dt_1);
      }

      // if grid is 2D/3D, compute minimum of (v2 +/- C)
      if (pmb->block_size.nx2 > 1) {
        for (int i=is; i<=ie; ++i) {
          Real& dt_2 = dt2(i);
          min_dt = std::min(min_dt,dt_2);
        }
      }

      // if grid is 3D, compute minimum of (v3 +/- C)
      if (pmb->block_size.nx3 > 1) {
        for (int i=is; i<=ie; ++i) {
          Real& dt_3 = dt3(i);
          min_dt = std::min(min_dt,dt_3);
        }
      }

    }
  }

// calculate the timestep limited by the diffusion process
  if (phdif->hydro_diffusion_defined) {
    Real mindt_vis, mindt_cnd;
    phdif->NewHydroDiffusionDt(mindt_vis, mindt_cnd);
    min_dt = std::min(min_dt,mindt_vis);
    min_dt = std::min(min_dt,mindt_cnd);
  } // hydro diffusion

  if(MAGNETIC_FIELDS_ENABLED &&
     pmb->pfield->pfdif->field_diffusion_defined) {
    Real mindt_oa, mindt_h;
    pmb->pfield->pfdif->NewFieldDiffusionDt(mindt_oa, mindt_h);
    min_dt = std::min(min_dt,mindt_oa);
    min_dt = std::min(min_dt,mindt_h);
  } // field diffusion

  min_dt *= pmb->pmy_mesh->cfl_number;

  if (UserTimeStep_!=NULL) {
    min_dt = std::min(min_dt, UserTimeStep_(pmb));
  }

  pmb->new_block_dt=min_dt;
  return min_dt;
}
