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

void Field::CT(const IntegratorWeight wght, FaceField &b_out)
{
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e1,e2,e3;
  e1.InitWithShallowCopy(pmb->pfield->e.x1e);
  e2.InitWithShallowCopy(pmb->pfield->e.x2e);
  e3.InitWithShallowCopy(pmb->pfield->e.x3e);

  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> area,len,len_p1;
  area.InitWithShallowSlice(face_area_,2,tid,1);
  len.InitWithShallowSlice(edge_length_,2,tid,1);
  len_p1.InitWithShallowSlice(edge_length_p1_,2,tid,1);

//---- update B1

  for (int k=ks; k<=ke; ++k) {
#pragma omp for schedule(static)
  for (int j=js; j<=je; ++j) {

    // add curl(E) in 2D and 3D problem
    if (pmb->block_size.nx2 > 1) {
      pmb->pcoord->Face1Area(k,j,is,ie+1,area);
      pmb->pcoord->Edge3Length(k,j  ,is,ie+1,len);
      pmb->pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b_out.x1f(k,j,i) -= wght.c*
           ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e3(k,j+1,i) - len(i)*e3(k,j,i));
      }

      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Edge2Length(k  ,j,is,ie+1,len);
        pmb->pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
#pragma simd
        for (int i=is; i<=ie+1; ++i) {
          b_out.x1f(k,j,i) += wght.c*
             ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e2(k+1,j,i) -len(i)*e2(k,j,i));
        }
      }
    }
  }}

//---- update B2 (in 1D and 3D problems)

  for (int k=ks; k<=ke; ++k) {
    // reset loop limits for polar boundary
    int jl=js; int ju=je+1;
    if (pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY
     || pmb->pbval->block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) jl=js+1;
    if (pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY
     || pmb->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) ju=je;
#pragma omp for schedule(static)
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->Face2Area(k,j,is,ie,area);
      pmb->pcoord->Edge3Length(k,j,is,ie+1,len);
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b_out.x2f(k,j,i) += (wght.c*(pmb->pmy_mesh->dt)/area(i))*(len(i+1)*e3(k,j,i+1)
                                                                  - len(i)*e3(k,j,i));
      }
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Edge1Length(k  ,j,is,ie,len);
        pmb->pcoord->Edge1Length(k+1,j,is,ie,len_p1);
#pragma simd
        for (int i=is; i<=ie; ++i) {
          b_out.x2f(k,j,i) -= wght.c*
             ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e1(k+1,j,i) - len(i)*e1(k,j,i));
        }
      }
    }
  }

//---- update B3 (in 1D and 2D problems)

  for (int k=ks; k<=ke+1; ++k) {
#pragma omp for schedule(static)
  for (int j=js; j<=je; ++j) {
    pmb->pcoord->Face3Area(k,j,is,ie,area);
    pmb->pcoord->Edge2Length(k,j,is,ie+1,len);
#pragma simd
    for (int i=is; i<=ie; ++i) {
      b_out.x3f(k,j,i) -= (wght.c*(pmb->pmy_mesh->dt)/area(i))*(len(i+1)*e2(k,j,i+1) -
                                                                len(i)*e2(k,j,i));
    }
    if (pmb->block_size.nx2 > 1) {
      pmb->pcoord->Edge1Length(k,j  ,is,ie,len);
      pmb->pcoord->Edge1Length(k,j+1,is,ie,len_p1);
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b_out.x3f(k,j,i) += wght.c*
           ((pmb->pmy_mesh->dt)/area(i))*(len_p1(i)*e1(k,j+1,i) - len(i)*e1(k,j,i));
      }
    }
  }}

} // end of OMP parallel region

  return;
}
