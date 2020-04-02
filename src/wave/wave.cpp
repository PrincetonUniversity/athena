//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave.cpp
//  \brief implementation of functions in the Wave class

// C++ headers
// #include <algorithm>  // min()
// #include <cmath>      // fabs(), sqrt()
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

#include "wave.hpp"

// constructor, initializes data structures and parameters

Wave::Wave(MeshBlock *pmb, ParameterInput *pin) :
  pmy_block(pmb),
#if PREFER_VC
  u(NWAVE_CPT,
    pmb->nverts3,
    pmb->nverts2,
    pmb->nverts1),
  coarse_u_(NWAVE_CPT,
            pmb->ncv3,
            pmb->ncv2,
            pmb->ncv1,
            (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
             AthenaArray<Real>::DataStatus::empty)),
#else
  u(NWAVE_CPT,
    pmb->ncells3,
    pmb->ncells2,
    pmb->ncells1),
  coarse_u_(NWAVE_CPT,
            pmb->ncc3,
            pmb->ncc2,
            pmb->ncc1,
            (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
             AthenaArray<Real>::DataStatus::empty)),
#endif
  empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
  ubvar(pmb, &u, &coarse_u_, empty_flux)
{
  printf("Wave spawned\n");

  Mesh *pm = pmb->pmy_mesh;
  Coordinates * pco = pmb->pcoord;
  // each wave obj. should have an associated block [now in initializer]
  // pmy_block = pmb;


  // dimensions required for data allocation
  if (PREFER_VC) {
    mbi.nn1 = pmb->nverts1;
    mbi.nn2 = pmb->nverts2;
    mbi.nn3 = pmb->nverts3;
  } else {
    mbi.nn1 = pmb->ncells1;
    mbi.nn2 = pmb->ncells2;
    mbi.nn3 = pmb->ncells3;
  }
  int nn1 = mbi.nn1, nn2 = mbi.nn2, nn3 = mbi.nn3;

  // convenience for per-block iteration (private Wave scope)
  mbi.il = pmb->is; mbi.jl = pmb->js; mbi.kl = pmb->ks;
  if (PREFER_VC) {
    mbi.iu = pmb->ive; mbi.ju = pmb->jve; mbi.ku = pmb->kve;
  } else {
    mbi.iu = pmb->ie; mbi.ju = pmb->je; mbi.ku = pmb->ke;
  }


  // point to appropriate grid
  if (PREFER_VC) {
    printf("slice coords for vc\n");
    mbi.x1.InitWithShallowSlice(pco->x1f, 1, 0, nn1);
    mbi.x2.InitWithShallowSlice(pco->x2f, 1, 0, nn2);
    mbi.x3.InitWithShallowSlice(pco->x3f, 1, 0, nn3);
  } else {
    printf("slice coords for cc\n");
    mbi.x1.InitWithShallowSlice(pco->x1v, 1, 0, nn1);
    mbi.x2.InitWithShallowSlice(pco->x2v, 1, 0, nn2);
    mbi.x3.InitWithShallowSlice(pco->x3v, 1, 0, nn3);
  }

  if (pm->multilevel) {
    // analogously construct coarse grid information as required
    mbi.cil = pmb->civs; mbi.cjl = pmb->cjvs; mbi.ckl = pmb->ckvs;

    if (PREFER_VC) {
      mbi.ciu = pmb->cive; mbi.cju = pmb->cjve; mbi.cku = pmb->ckve;

      mbi.cnn1 = pmb->ncv1; mbi.cnn2 = pmb->ncv2; mbi.cnn3 = pmb->ncv3;
    } else {
      mbi.ciu = pmb->cie; mbi.cju = pmb->cje; mbi.cku = pmb->cke;

      mbi.cnn1 = pmb->ncc1; mbi.cnn2 = pmb->ncc2; mbi.cnn3 = pmb->ncc3;
    }
  }

  // point to appropriate grid
  if (FILL_WAVE_COARSE_P)
  if (PREFER_VC) {
    printf("slice coords for vc\n");
    mbi.cx1.InitWithShallowSlice(pmb->pmr->pcoarsec->x1f, 0, 1);
    mbi.cx2.InitWithShallowSlice(pmb->pmr->pcoarsec->x2f, 0, 1);
    mbi.cx3.InitWithShallowSlice(pmb->pmr->pcoarsec->x3f, 0, 1);
  } else {
    printf("slice coords for cc\n");
    mbi.cx1.InitWithShallowSlice(pmb->pmr->pcoarsec->x1v, 0, 1);
    mbi.cx2.InitWithShallowSlice(pmb->pmr->pcoarsec->x2v, 0, 1);
    mbi.cx3.InitWithShallowSlice(pmb->pmr->pcoarsec->x3v, 0, 1);
  }

  // inform MeshBlock that this array is the "primary" representation
  // Used for:
  // (1) load-balancing
  // (2) (future) dumping to restart file
  pmb->RegisterMeshBlockData(u);


  // Allocate memory for the solution and its time derivative
  // u.NewAthenaArray(NWAVE_CPT, nn3, nn2, nn1);
  u1.NewAthenaArray(NWAVE_CPT, nn3, nn2, nn1);
  rhs.NewAthenaArray(NWAVE_CPT, nn3, nn2, nn1);

  exact.NewAthenaArray(nn3, nn2, nn1);
  error.NewAthenaArray(nn3, nn2, nn1);

  // If user-requested time integrator is type 3S* allocate additional memory
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4")
    u2.NewAthenaArray(NWAVE_CPT, nn3, nn2, nn1);

  c = pin->GetOrAddReal("wave", "c", 1.0);
  use_Sommerfeld = pin->GetOrAddInteger("wave", "use_Sommerfeld", 0);

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmb->pmr->AddToRefinement(&u, &coarse_u_);
  }

  // enroll CellCenteredBoundaryVariable / VertexCenteredBoundaryVariable object
  ubvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&ubvar);
  pmb->pbval->bvars_main_int.push_back(&ubvar);

  // Allocate memory for scratch arrays
  dt1_.NewAthenaArray(nn1);
  dt2_.NewAthenaArray(nn1);
  dt3_.NewAthenaArray(nn1);

  // Set up finite difference operators
  Real dx1, dx2, dx3;
  if (PREFER_VC) {
    dx1 = pco->dx1f(0); dx2 = pco->dx2f(0); dx3 = pco->dx3f(0);
  } else {
    dx1 = pco->dx1v(0); dx2 = pco->dx2v(0); dx3 = pco->dx3v(0);
  }

  FD.stride[0] = 1;
  FD.stride[1] = 0;
  FD.stride[2] = 0;
  FD.idx[0] = 1.0 / dx1;
  FD.idx[1] = 0.0;
  FD.idx[2] = 0.0;
  if(nn2 > 1) {
    FD.stride[1] = nn1;
    FD.idx[1] = 1.0 / dx2;
  }

  if(nn3 > 1) {
    FD.stride[2] = nn2 * nn1;
    FD.idx[2] = 1.0 / dx3;
  }

}


// destructor
Wave::~Wave()
{
  u.DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  u2.DeleteAthenaArray(); // only allocated in case of 3S*-type of integrator

  rhs.DeleteAthenaArray();

  exact.DeleteAthenaArray();
  error.DeleteAthenaArray();

  // note: do not include x1_, x2_, x3_
}


//----------------------------------------------------------------------------------------
//! \fn  void Wave::AddWaveRHS
//  \brief Adds RHS to weighted average of variables from
//  previous step(s) of time integrator algorithm

void Wave::AddWaveRHS(const Real wght, AthenaArray<Real> &u_out) {
  MeshBlock *pmb=pmy_block;

  for (int k=mbi.kl; k<=mbi.ku; ++k) {
    for (int j=mbi.jl; j<=mbi.ju; ++j) {
      // update variables
      for (int n=0; n<NWAVE_CPT; ++n) {
#pragma omp simd
        for (int i=mbi.il; i<=mbi.iu; ++i) {
          u_out(n, k, j, i) += wght*(pmb->pmy_mesh->dt)*rhs(n, k, j, i);
        }
        }
      }
    }

  return;
}
