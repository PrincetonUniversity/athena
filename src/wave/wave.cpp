//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave.cpp
//  \brief implementation of functions in the Wave class

// C++ headers
#include <iostream>
#include <string>
// #include <algorithm>  // min()
// #include <cmath>      // fabs(), sqrt()
#ifdef COMPACT_FD
#include <vector>
#endif // COMPACT_FD
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
  empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
  ubvar(pmb, &u, &coarse_u_, empty_flux)
{
  Mesh *pm = pmb->pmy_mesh;
  Coordinates * pco = pmb->pcoord;

  // dimensions required for data allocation
  mbi.nn1 = pmb->nverts1;
  mbi.nn2 = pmb->nverts2;
  mbi.nn3 = pmb->nverts3;
  int nn1 = mbi.nn1, nn2 = mbi.nn2, nn3 = mbi.nn3;

  // convenience for per-block iteration (private Wave scope)
  mbi.il = pmb->is; mbi.jl = pmb->js; mbi.kl = pmb->ks;
  mbi.iu = pmb->ive; mbi.ju = pmb->jve; mbi.ku = pmb->kve;

  // point to appropriate grid
  mbi.x1.InitWithShallowSlice(pco->x1f, 1, 0, nn1);
  mbi.x2.InitWithShallowSlice(pco->x2f, 1, 0, nn2);
  mbi.x3.InitWithShallowSlice(pco->x3f, 1, 0, nn3);

  // BD: debug remove after tests pass
  // if (pm->multilevel) {
  //   // analogously construct coarse grid information as required
  //   mbi.cil = pmb->civs; mbi.cjl = pmb->cjvs; mbi.ckl = pmb->ckvs;

  //   mbi.ciu = pmb->cive; mbi.cju = pmb->cjve; mbi.cku = pmb->ckve;

  //   mbi.cnn1 = pmb->ncv1; mbi.cnn2 = pmb->ncv2; mbi.cnn3 = pmb->ncv3;
  // }

  // point to appropriate grid
#ifdef FILL_WAVE_COARSE_P
  printf("slice coords for vc\n");
  mbi.cx1.InitWithShallowSlice(pmb->pmr->pcoarsec->x1f, 0, 1);
  mbi.cx2.InitWithShallowSlice(pmb->pmr->pcoarsec->x2f, 0, 1);
  mbi.cx3.InitWithShallowSlice(pmb->pmr->pcoarsec->x3f, 0, 1);
#endif // FILL_WAVE_COARSE_P
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

  // debug control
  debug_inspect_error = pin->GetOrAddBoolean("wave", "debug_inspect_err", false);
  debug_abort_threshold = pin->GetOrAddReal("wave", "debug_abort_threshold",
    std::numeric_limits<Real>::max());

  // Additional boundary condition control
  std::string boundary_type = pin->GetOrAddString("wave", "boundary_type", "none");

  if (boundary_type == "Dirichlet") {
    use_Dirichlet = true;

    if (pmb->pmy_mesh->ndim == 1) {

      Lx1_ = pin->GetOrAddReal("wave", "Lx1", 1.);
      M_ = pin->GetOrAddInteger("wave", "M_cutoff", 1.);

      A_.NewAthenaArray(M_);
      B_.NewAthenaArray(M_);

      for (int m=1; m<=M_; ++m) {
        std::string m_str = std::to_string(m);
        A_(m-1) = pin->GetOrAddReal("wave", "A_" + m_str, 0.);
        B_(m-1) = pin->GetOrAddReal("wave", "B_" + m_str, 0.);
      }

    } else if (pmb->pmy_mesh->ndim == 2) {

      Lx1_ = pin->GetOrAddReal("wave", "Lx1", 1.);
      Lx2_ = pin->GetOrAddReal("wave", "Lx2", 1.);

      M_ = pin->GetOrAddInteger("wave", "M_cutoff", 1.);
      N_ = pin->GetOrAddInteger("wave", "N_cutoff", 1.);

      A_.NewAthenaArray(M_, N_);
      B_.NewAthenaArray(M_, N_);

      for (int n=1; n<=N_; ++n) {
        std::string n_str = std::to_string(n);

        for (int m=1; m<=M_; ++m) {
          std::string m_str = std::to_string(m);
          A_(m-1, n-1) = pin->GetOrAddReal("wave", "A_" + m_str + n_str, 0.);
          B_(m-1, n-1) = pin->GetOrAddReal("wave", "B_" + m_str + n_str, 0.);
        }
      }

    } else {
      Lx1_ = pin->GetOrAddReal("wave", "Lx1", 1.);
      Lx2_ = pin->GetOrAddReal("wave", "Lx2", 1.);
      Lx3_ = pin->GetOrAddReal("wave", "Lx3", 1.);

      M_ = pin->GetOrAddInteger("wave", "M_cutoff", 1.);
      N_ = pin->GetOrAddInteger("wave", "N_cutoff", 1.);
      O_ = pin->GetOrAddInteger("wave", "O_cutoff", 1.);

      A_.NewAthenaArray(M_, N_, O_);
      B_.NewAthenaArray(M_, N_, O_);

      for (int o=1; o<=O_; ++o) {
        std::string o_str = std::to_string(o);
        for (int n=1; n<=N_; ++n) {
          std::string n_str = std::to_string(n);
          for (int m=1; m<=M_; ++m) {
            std::string m_str = std::to_string(m);
            A_(m-1, n-1, o-1) = pin->GetOrAddReal(
              "wave", "A_" + m_str + n_str + o_str, 0.);
            B_(m-1, n-1, o-1) = pin->GetOrAddReal(
              "wave", "B_" + m_str + n_str + o_str, 0.);
          }
        }
      }

    }

  } else if (boundary_type == "Sommerfeld")
    use_Sommerfeld = true;

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmb->pmr->AddToRefinement(&u, &coarse_u_);
  }

  // enroll CellCenteredBoundaryVariable / VertexCenteredBoundaryVariable object
  ubvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&ubvar);
  pmb->pbval->bvars_main_int_vc.push_back(&ubvar);

  // Allocate memory for scratch arrays
  dt1_.NewAthenaArray(nn1);
  dt2_.NewAthenaArray(nn1);
  dt3_.NewAthenaArray(nn1);

  // Set up finite difference operators
  Real dx1, dx2, dx3;
  dx1 = pco->dx1f(0); dx2 = pco->dx2f(0); dx3 = pco->dx3f(0);

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

#ifdef COMPACT_FD
  // Set up compact finite difference operators
  const unsigned int order = 2;
  if (pmb->pmy_mesh->ndim == 3) {
    const std::vector<unsigned int> dims_N {
      (unsigned int) nn1, (unsigned int) nn2, (unsigned int) nn3
    };
    const std::vector<Real> dx {
      dx1, dx2, dx3
    };
    FDC2 = new FDCompact(order, dims_N, dx);

  } else if (pmb->pmy_mesh->ndim == 2) {
    const std::vector<unsigned int> dims_N {
      (unsigned int) nn1, (unsigned int) nn2
    };
    const std::vector<Real> dx {
      dx1, dx2
    };

    FDC2 = new FDCompact(order, dims_N, dx);
  } else if (pmb->pmy_mesh->ndim == 1) {
    const std::vector<unsigned int> dims_N {(unsigned int) nn1};
    const std::vector<Real> dx {
      dx1
    };
    FDC2 = new FDCompact(order, dims_N, dx);
  }
  scratch_wu.NewAthenaArray(nn3, nn2, nn1);
#endif // COMPACT_FD

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

  if (use_Dirichlet) {
    A_.DeleteAthenaArray();
    B_.DeleteAthenaArray();
  }

#ifdef COMPACT_FD
  delete FDC2;
  scratch_wu.DeleteAthenaArray();
#endif // COMPACT_FD

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
