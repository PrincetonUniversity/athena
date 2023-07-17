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
//! \file cr_source.cpp
//  \brief Add cosmic ray source terms to  gas
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/utils.hpp"
#include "../cr.hpp"

// class header
#include "cr_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// add the source terms implicitly
void CRIntegrator::AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc,
        AthenaArray<Real> &u_cr) {
  CosmicRay *pcr=pmb->pcr;
  Real vlim = pcr->vmax;
  Real invlim = 1.0/vlim;
  Real rho_floor = pmb->peos->GetDensityFloor();

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
         // Real fxx = 1.0/3.0;
         // Real fyy = 1.0/3.0;
         // Real fzz = 1.0/3.0;
         // Real fxy = 0.0;
         // Real fxz = 0.0;
         // Real fyz = 0.0;

         Real *ec = &(u_cr(CRE,k,j,0));
         Real *fc1 = &(u_cr(CRF1,k,j,0));
         Real *fc2 = &(u_cr(CRF2,k,j,0));
         Real *fc3 = &(u_cr(CRF3,k,j,0));

         // The angle of B
         Real *sint_b = &(pcr->b_angle(0,k,j,0));
         Real *cost_b = &(pcr->b_angle(1,k,j,0));
         Real *sinp_b = &(pcr->b_angle(2,k,j,0));
         Real *cosp_b = &(pcr->b_angle(3,k,j,0));

         // adv1 is dPc/dx, adv2 is dPc/dy, adv3 is dPc/dz

      for (int i=is; i<=ie; ++i) {
         Real rho = u(IDN,k,j,i);
         rho = std::max(rho,rho_floor);
         Real v1 = u(IM1,k,j,i)/rho;
         Real v2 = u(IM2,k,j,i)/rho;
         Real v3 = u(IM3,k,j,i)/rho;
         Real vtot1 = v1;
         Real vtot2 = v2;
         Real vtot3 = v3;

         // add the streaming velocity
         if (pcr->stream_flag) {
           vtot1 += pcr->v_adv(0,k,j,i);
           vtot2 += pcr->v_adv(1,k,j,i);
           vtot3 += pcr->v_adv(2,k,j,i);
         }

         Real fr1 = fc1[i];
         Real fr2 = fc2[i];
         Real fr3 = fc3[i];

         // in the case with magnetic field
        // rotate the vectors to oriante to the B direction
         if (MAGNETIC_FIELDS_ENABLED) {
           // Apply rotation of the vectors
           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],v1,v2,v3);

           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],fr1,fr2,fr3);

           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],vtot1,vtot2,vtot3);

           // perpendicular energy source term is already added via ec_source_
           // we will still need to include this in general for diffusion case
           vtot2 = 0.0;
           vtot3 = 0.0;
         }

         Real sigma_x = pcr->sigma_diff(0,k,j,i);
         Real sigma_y = pcr->sigma_diff(1,k,j,i);
         Real sigma_z = pcr->sigma_diff(2,k,j,i);

         if (pcr->stream_flag) {
           sigma_x = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) +
                             1.0/pcr->sigma_adv(0,k,j,i));

           sigma_y = 1.0/(1.0/pcr->sigma_diff(1,k,j,i) +
                             1.0/pcr->sigma_adv(1,k,j,i));

           sigma_z = 1.0/(1.0/pcr->sigma_diff(2,k,j,i) +
                             1.0/pcr->sigma_adv(2,k,j,i));
         }

         // Now update the momentum equation
         //\partial F/\partial t=-V_m\sigma (F-v(E+Pc_)/v_m))
         // And the energy equation
         //\partial E_c/\partial t = -vtot sigma (F- v(E_c+P_c)/v_m)

         Real rhs1 = ec[i];
         Real rhs2 = fr1;
         Real rhs3 = fr2;
         Real rhs4 = fr3;

         Real coef_11 = 1.0 - dt * sigma_x * vtot1 * v1 * invlim * 4.0/3.0
                            - dt * sigma_y * vtot2 * v2 * invlim * 4.0/3.0
                            - dt * sigma_z * vtot3 * v3 * invlim * 4.0/3.0;
         Real coef_12 = dt * sigma_x * vtot1;
         Real coef_13 = dt * sigma_y * vtot2;
         Real coef_14 = dt * sigma_z * vtot3;

         Real coef_21 = -dt * v1 * sigma_x * 4.0/3.0;
         Real coef_22 = 1.0 + dt * vlim * sigma_x;

         Real coef_31 = -dt * v2 * sigma_y * 4.0/3.0;
         Real coef_33 = 1.0 + dt * vlim * sigma_y;

         Real coef_41 = -dt * v3 * sigma_z * 4.0/3.0;
         Real coef_44 = 1.0 + dt * vlim * sigma_z;

        //newfr1 = (rhs2 - coef21 * newEc)/coef22
        // newfr2= (rhs3 - coef31 * newEc)/coef33
        // newfr3 = (rhs4 - coef41 * newEc)/coef44
        // coef11 - coef21 * coef12 /coef22 - coef13 * coef31 /
         //              coef33 - coef41 * coef14 /coef44)* newec
        //    =rhs1 - coef12 *rhs2/coef22 - coef13 * rhs3/coef33 - coef14 * rhs4/coef44

        Real e_coef = coef_11 - coef_12 * coef_21/coef_22 - coef_13 * coef_31/coef_33
                       - coef_14 * coef_41/coef_44;
        Real new_ec = rhs1 - coef_12 * rhs2/coef_22 - coef_13 * rhs3/coef_33
                      - coef_14 * rhs4/coef_44;
        new_ec /= e_coef;

        Real newfr1 = (rhs2 - coef_21 * new_ec)/coef_22;
        Real newfr2 = (rhs3 - coef_31 * new_ec)/coef_33;
        Real newfr3 = (rhs4 - coef_41 * new_ec)/coef_44;


        // Now apply the invert rotation
        if (MAGNETIC_FIELDS_ENABLED) {
         // Apply rotation of the vectors
          InvRotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],
                                         newfr1,newfr2,newfr3);
          new_ec += dt * ec_source_(k,j,i);
        }

         // Add the energy source term
         if (NON_BAROTROPIC_EOS && (pcr->src_flag > 0)) {
           Real new_eg = u(IEN,k,j,i) - (new_ec - ec[i]);
           if (new_eg < 0.0) new_eg = u(IEN,k,j,i);
           u(IEN,k,j,i) = new_eg;
         }

         if (new_ec < 0.0) new_ec = ec[i];

         if (pcr->src_flag > 0) {
           u(IM1,k,j,i) += (-(newfr1 - fc1[i]) * invlim);
           u(IM2,k,j,i) += (-(newfr2 - fc2[i]) * invlim);
           u(IM3,k,j,i) += (-(newfr3 - fc3[i]) * invlim);
         }
         u_cr(CRE,k,j,i) = new_ec;
         u_cr(CRF1,k,j,i) = newfr1;
         u_cr(CRF2,k,j,i) = newfr2;
         u_cr(CRF3,k,j,i) = newfr3;
      }
    }
  }
  // Add user defined source term for cosmic rays
  if (pcr->cr_source_defined)
    pcr->UserSourceTerm_(pmb, pmb->pmy_mesh->time, dt, w,pmb->pfield->b, u_cr);
}
