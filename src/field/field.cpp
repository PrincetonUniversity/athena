//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field.cpp
//  \brief implementation of functions in class Field

// C headers

// C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "field.hpp"
#include "field_diffusion/field_diffusion.hpp"

// constructor, initializes data structures and parameters

Field::Field(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for interface fields
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  //  Note the extra cell in each longitudinal dirn for interface fields
  b.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
  b.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
  b.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  pmb->RegisterMeshBlockData(b);

  b1.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
  b1.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
  b1.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time","integrator","vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED) {
    // future extension may add "int nregister" to Hydro class
    b2.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b2.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b2.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  }

  bcc.NewAthenaArray (NFIELD,ncells3,ncells2,ncells1);

  e.x1e.NewAthenaArray((ncells3+1),(ncells2+1), ncells1   );
  e.x2e.NewAthenaArray((ncells3+1), ncells2   ,(ncells1+1));
  e.x3e.NewAthenaArray( ncells3   ,(ncells2+1),(ncells1+1));

  wght.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
  wght.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
  wght.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

  e2_x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
  e3_x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
  e1_x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
  e3_x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
  e1_x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  e2_x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

  // Allocate memory for scratch vectors
  if (pmb->block_size.nx3 == 1)
    cc_e_.NewAthenaArray(ncells3,ncells2,ncells1);
  else
    cc_e_.NewAthenaArray(3,ncells3,ncells2,ncells1);

  face_area_.NewAthenaArray(ncells1);
  edge_length_.NewAthenaArray(ncells1);
  edge_length_p1_.NewAthenaArray(ncells1);
  if (GENERAL_RELATIVITY) {
    g_.NewAthenaArray(NMETRIC, ncells1);
    gi_.NewAthenaArray(NMETRIC, ncells1);
  }
  // allocate prolongation buffers
  if (pmy_block->pmy_mesh->multilevel == true) {
    int ncc1 = pmb->block_size.nx1/2 + 2*NGHOST;
    int ncc2 = 1;
    if (pmb->block_size.nx2 > 1) ncc2 = pmb->block_size.nx2/2 + 2*NGHOST;
    int ncc3 = 1;
    if (pmb->block_size.nx3 > 1) ncc3 = pmb->block_size.nx3/2 + 2*NGHOST;
    coarse_b_.x1f.NewAthenaArray(ncc3, ncc2, ncc1+1);
    coarse_b_.x2f.NewAthenaArray(ncc3, ncc2+1, ncc1);
    coarse_b_.x3f.NewAthenaArray(ncc3+1, ncc2, ncc1);
    coarse_bcc_.NewAthenaArray(3, ncc3, ncc2, ncc1);
    // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
    pmy_block->pmr->AddToRefinement(&b, &coarse_b_);
  }
  // ptr to diffusion object
  pfdif = new FieldDiffusion(pmb, pin);
  // create object to interface with BoundaryValues
  pfbval  = new FaceCenteredBoundaryVariable(pmy_block, &b, coarse_b_, e);
  pfbval->bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(pfbval);
  pmb->pbval->bvars_main_int.push_back(pfbval);
}

// destructor

Field::~Field() {
  delete pfdif;
  delete pfbval;
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

void Field::CalculateCellCenteredField(
    const FaceField &bf, AthenaArray<Real> &bc, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  // Defer to Reconstruction class to check if uniform Cartesian formula can be used
  // (unweighted average)
  const bool uniform_ave_x1 = pmy_block->precon->uniform_limiter[X1DIR];
  const bool uniform_ave_x2 = pmy_block->precon->uniform_limiter[X2DIR];
  const bool uniform_ave_x3 = pmy_block->precon->uniform_limiter[X3DIR];

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
        if (uniform_ave_x1 == true) {
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

        if (uniform_ave_x2 == true) {
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
        if (uniform_ave_x3 == true) {
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
