//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ct.cpp
//  \brief

// C++ headers
#include <algorithm>  // max(), min()

// Athena++ headers
#include "field.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../bvals/bvals.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Field::CT
//  \brief Constrained Transport implementation of dB/dt = -Curl(E), where E=-(v X B)

void Field::CT(const Real wght, FaceField &b_out) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1,e2,e3;
  e1.InitWithShallowCopy(pmb->pfield->e.x1e);
  e2.InitWithShallowCopy(pmb->pfield->e.x2e);
  e3.InitWithShallowCopy(pmb->pfield->e.x3e);

  AthenaArray<Real> area,len,len_p1;
  area.InitWithShallowCopy(face_area_);
  len.InitWithShallowCopy(edge_length_);
  len_p1.InitWithShallowCopy(edge_length_p1_);

//---- update B1

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {

    // add curl(E) in 2D and 3D problem
    if (pmb->block_size.nx2 > 1) {
      pmb->pcoord->Face1Area(k,j,is,ie+1,area);
      pmb->pcoord->Edge3Length(k,j  ,is,ie+1,len);
      pmb->pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b_out.x1f(k,j,i) -= wght*
           ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e3(k,j+1,i) - len(i)*e3(k,j,i));
      }

      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Edge2Length(k  ,j,is,ie+1,len);
        pmb->pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b_out.x1f(k,j,i) += wght*
             ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e2(k+1,j,i) -len(i)*e2(k,j,i));
        }
      }
    }
  }}

//---- update B2 (curl terms in 1D and 3D problems)

  for (int k=ks; k<=ke; ++k) {
    // reset loop limits for polar boundary
    int jl=js; int ju=je+1;
    if (pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY
     || pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) jl=js+1;
    if (pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY
     || pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) ju=je;
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->Face2Area(k,j,is,ie,area);
      pmb->pcoord->Edge3Length(k,j,is,ie+1,len);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_out.x2f(k,j,i) += (wght*(pmb->pmy_mesh->dt)/area(i))*(len(i+1)*e3(k,j,i+1)
                                                                  - len(i)*e3(k,j,i));
      }
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Edge1Length(k  ,j,is,ie,len);
        pmb->pcoord->Edge1Length(k+1,j,is,ie,len_p1);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x2f(k,j,i) -= wght*
             ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e1(k+1,j,i) - len(i)*e1(k,j,i));
        }
      }
    }
  }

//---- update B3 (curl terms in 1D and 2D problems)

  for (int k=ks; k<=ke+1; ++k) {
  for (int j=js; j<=je; ++j) {
    pmb->pcoord->Face3Area(k,j,is,ie,area);
    pmb->pcoord->Edge2Length(k,j,is,ie+1,len);
#pragma omp simd
    for (int i=is; i<=ie; ++i) {
      b_out.x3f(k,j,i) -= (wght*(pmb->pmy_mesh->dt)/area(i))*(len(i+1)*e2(k,j,i+1) -
                                                                len(i)*e2(k,j,i));
    }
    if (pmb->block_size.nx2 > 1) {
      pmb->pcoord->Edge1Length(k,j  ,is,ie,len);
      pmb->pcoord->Edge1Length(k,j+1,is,ie,len_p1);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_out.x3f(k,j,i) += wght*
           ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e1(k,j+1,i) - len(i)*e1(k,j,i));
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Field::WeightedAveB
//  \brief Compute weighted average of face-averaged B in time integrator step

void Field::WeightedAveB(FaceField &b_out, FaceField &b_in1, FaceField &b_in2,
                         const Real wght[3]) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // Note: these loops can be combined now that they avoid curl terms
  // Only need to separately account for the final longitudinal face in each loop limit

  // b_in2 may be an unallocated AthenaArray if using a 2S time integrator
  if (wght[2] != 0.0) {
//---- B1
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i)
              + wght[2]*b_in2.x1f(k,j,i);
        }
      }}

//---- B2

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      // move these limit modifications outside the loop
      if (pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY
          || pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) jl=js+1;
      if (pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY
          || pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) ju=je;
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i)
              + wght[2]*b_in2.x2f(k,j,i);
        }
      }
    }

//---- B3

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i)
              + wght[2]*b_in2.x3f(k,j,i);
        }
      }}
  } else { // do not derefernce b_in2
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i);
        }
      }}

//---- B2

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      // move these limit modifications outside the loop
      if (pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY
          || pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) jl=js+1;
      if (pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY
          || pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) ju=je;
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i);
        }
      }
    }

//---- B3

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i);
        }
      }}
  }
  return;
}
