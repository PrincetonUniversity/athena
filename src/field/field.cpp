//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field.cpp
//! \brief implementation of functions in class Field

// C headers

// C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "field.hpp"
#include "field_diffusion/field_diffusion.hpp"

//! constructor, initializes data structures and parameters

Field::Field(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb), b(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    b1(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    bcc(NFIELD, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    e(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    wght(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    e2_x1f( pmb->ncells3   , pmb->ncells2   ,(pmb->ncells1+1)),
    e3_x1f( pmb->ncells3   , pmb->ncells2   ,(pmb->ncells1+1)),
    e1_x2f( pmb->ncells3   ,(pmb->ncells2+1), pmb->ncells1   ),
    e3_x2f( pmb->ncells3   ,(pmb->ncells2+1), pmb->ncells1   ),
    e1_x3f((pmb->ncells3+1), pmb->ncells2   , pmb->ncells1   ),
    e2_x3f((pmb->ncells3+1), pmb->ncells2   , pmb->ncells1   ),
    coarse_bcc_(3, pmb->ncc3, pmb->ncc2, pmb->ncc1,
                (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
                 AthenaArray<Real>::DataStatus::empty)),
    coarse_b_(pmb->ncc3, pmb->ncc2, pmb->ncc1+1,
              (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
               AthenaArray<Real>::DataStatus::empty)),
    fbvar(pmb, &b, coarse_b_, e),
    fdif(pmb, pin) {
  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, ncells3 = pmb->ncells3;
  Mesh *pm = pmy_block->pmy_mesh;

  pmb->RegisterMeshBlockData(b);

  // If user-requested time integrator is type 3S*, allocate additional memory registers
  // Note the extra cell in each longitudinal direction for interface fields
  std::string integrator = pin->GetOrAddString("time","integrator","vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED) {
    // future extension may add "int nregister" to Hydro class
    b2.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b2.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b2.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  }

  if (STS_ENABLED) {
    std::string sts_integrator = pin->GetOrAddString("time", "sts_integrator", "rkl2");
    if (sts_integrator == "rkl2") {
      b0.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
      b0.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
      b0.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
      ct_update.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
      ct_update.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
      ct_update.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    }
  }

  // Allocate memory for scratch vectors
  if (!pm->f3)
    cc_e_.NewAthenaArray(ncells3, ncells2, ncells1);
  else
    cc_e_.NewAthenaArray(3, ncells3, ncells2, ncells1);

  face_area_.NewAthenaArray(ncells1);
  edge_length_.NewAthenaArray(ncells1);
  edge_length_p1_.NewAthenaArray(ncells1);
  if (GENERAL_RELATIVITY) {
    g_.NewAthenaArray(NMETRIC, ncells1);
    gi_.NewAthenaArray(NMETRIC, ncells1);
  }

  if (pm->multilevel) {
    // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
    refinement_idx = pmy_block->pmr->AddToRefinement(&b, &coarse_b_);
  }

  // enroll FaceCenteredBoundaryVariable object
  fbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&fbvar);
  pmb->pbval->bvars_main_int.push_back(&fbvar);
  if (STS_ENABLED) {
    if (fdif.field_diffusion_defined) {
      if (!pmb->phydro->hdif.hydro_diffusion_defined && NON_BAROTROPIC_EOS) {
        pmb->pbval->bvars_sts.push_back(pmb->pbval->bvars_main_int[0]);
      }
      pmb->pbval->bvars_sts.push_back(&fbvar);
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void Field::CalculateCellCenteredField
//! \brief cell center B-fields are defined as spatial interpolation at the volume center

void Field::CalculateCellCenteredField(
    const FaceField &bf, AthenaArray<Real> &bc, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  // Defer to Reconstruction class to check if uniform Cartesian formula can be used
  // (unweighted average)
  const bool uniform_ave_x1 = pmy_block->precon->uniform[X1DIR];
  const bool uniform_ave_x2 = pmy_block->precon->uniform[X2DIR];
  const bool uniform_ave_x3 = pmy_block->precon->uniform[X3DIR];

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // calc cell centered fields first
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        const Real& b1_i   = bf.x1f(k,j,i  );
        const Real& b1_ip1 = bf.x1f(k,j,i+1);
        const Real& b2_j   = bf.x2f(k,j  ,i);
        const Real& b2_jp1 = bf.x2f(k,j+1,i);
        const Real& b3_k   = bf.x3f(k  ,j,i);
        const Real& b3_kp1 = bf.x3f(k+1,j,i);

        Real& bcc1 = bc(IB1,k,j,i);
        Real& bcc2 = bc(IB2,k,j,i);
        Real& bcc3 = bc(IB3,k,j,i);
        Real lw, rw; // linear interpolation coefficients from lower and upper cell faces

        // cell center B-fields are defined as spatial interpolation at the volume center
        if (uniform_ave_x1) {
          lw = 0.5;
          rw = 0.5;
        } else {
          const Real& x1f_i  = pco->x1f(i);
          const Real& x1f_ip = pco->x1f(i+1);
          const Real& x1v_i  = pco->x1v(i);
          const Real& dx1_i  = pco->dx1f(i);
          lw = (x1f_ip - x1v_i)/dx1_i;
          rw = (x1v_i  - x1f_i)/dx1_i;
        }
        bcc1 = lw*b1_i + rw*b1_ip1;

        if (uniform_ave_x2) {
          lw = 0.5;
          rw = 0.5;
        } else {
          const Real& x2f_j  = pco->x2f(j);
          const Real& x2f_jp = pco->x2f(j+1);
          const Real& x2v_j  = pco->x2v(j);
          const Real& dx2_j  = pco->dx2f(j);
          lw = (x2f_jp - x2v_j)/dx2_j;
          rw = (x2v_j  - x2f_j)/dx2_j;
        }
        bcc2 = lw*b2_j + rw*b2_jp1;
        if (uniform_ave_x3) {
          lw = 0.5;
          rw = 0.5;
        } else {
          const Real& x3f_k  = pco->x3f(k);
          const Real& x3f_kp = pco->x3f(k+1);
          const Real& x3v_k  = pco->x3v(k);
          const Real& dx3_k  = pco->dx3f(k);
          lw = (x3f_kp - x3v_k)/dx3_k;
          rw = (x3v_k  - x3f_k)/dx3_k;
        }
        bcc3 = lw*b3_k + rw*b3_kp1;
      }
    }
  }
  return;
}
