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
//! \file new_blockdt.cpp
//  \brief computes timestep using CFL condition on a MEshBlock
//======================================================================================

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"

// this class header
#include "hydro.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
// \!fn Real Hydro::NewBlockTimeStep(void)
// \brief calculate the minimum timestep within a MeshBlock 

Real Hydro::NewBlockTimeStep(void)
{
  MeshBlock *pmb=pmy_block;
  int tid=0;
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

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
  Real *pthread_min_dt;
  pthread_min_dt = new Real [nthreads];

  for (int n=0; n<nthreads; ++n) pthread_min_dt[n] = (FLT_MAX);

#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> dt1, dt2, dt3;
  dt1.InitWithShallowSlice(dt1_,2,tid,1);
  dt2.InitWithShallowSlice(dt2_,2,tid,1);
  dt3.InitWithShallowSlice(dt3_,2,tid,1);
  Real wi[(NWAVE)];

  for (int k=ks; k<=ke; ++k){

#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j){
//#pragma simd
      for (int i=is; i<=ie; ++i){
        wi[IDN]=w(IDN,k,j,i);
        wi[IVX]=w(IVX,k,j,i);
        wi[IVY]=w(IVY,k,j,i);
        wi[IVZ]=w(IVZ,k,j,i);
        if (NON_BAROTROPIC_EOS) wi[IPR]=w(IPR,k,j,i);

        if (RELATIVISTIC_DYNAMICS) {

          dt1(i) = pmy_block->pcoord->CenterWidth1(k,j,i);
          dt2(i) = pmy_block->pcoord->CenterWidth2(k,j,i);
          dt3(i) = pmy_block->pcoord->CenterWidth3(k,j,i);

        } else if (MAGNETIC_FIELDS_ENABLED) {

          Real bx = bcc(IB1,k,j,i) + fabs(b_x1f(k,j,i)-bcc(IB1,k,j,i));
          wi[IBY] = bcc(IB2,k,j,i);
          wi[IBZ] = bcc(IB3,k,j,i);
          Real cf = pmy_block->peos->FastMagnetosonicSpeed(wi,bx);
          dt1(i)= pmy_block->pcoord->CenterWidth1(k,j,i)/(fabs(wi[IVX]) + cf);

          wi[IBY] = bcc(IB3,k,j,i);
          wi[IBZ] = bcc(IB1,k,j,i);
          bx = bcc(IB2,k,j,i) + fabs(b_x2f(k,j,i)-bcc(IB2,k,j,i));
          cf = pmy_block->peos->FastMagnetosonicSpeed(wi,bx);
          dt2(i)= pmy_block->pcoord->CenterWidth2(k,j,i)/(fabs(wi[IVY]) + cf);

          wi[IBY] = bcc(IB1,k,j,i);
          wi[IBZ] = bcc(IB2,k,j,i);
          bx = bcc(IB3,k,j,i) + fabs(b_x3f(k,j,i)-bcc(IB3,k,j,i));
          cf = pmy_block->peos->FastMagnetosonicSpeed(wi,bx);
          dt3(i)= pmy_block->pcoord->CenterWidth3(k,j,i)/(fabs(wi[IVZ]) + cf);

        } else {

          Real cs = pmy_block->peos->SoundSpeed(wi);
          dt1(i)= pmy_block->pcoord->CenterWidth1(k,j,i)/(fabs(wi[IVX]) + cs);
          dt2(i)= pmy_block->pcoord->CenterWidth2(k,j,i)/(fabs(wi[IVY]) + cs);
          dt3(i)= pmy_block->pcoord->CenterWidth3(k,j,i)/(fabs(wi[IVZ]) + cs);

        }
      }

// compute minimum of (v1 +/- C)

      for (int i=is; i<=ie; ++i){
        Real& dt_1 = dt1(i);
        pthread_min_dt[tid] = std::min(pthread_min_dt[tid],dt_1);
      }
    
// if grid is 2D/3D, compute minimum of (v2 +/- C)

      if (pmb->block_size.nx2 > 1) {
        for (int i=is; i<=ie; ++i){
          Real& dt_2 = dt2(i);
          pthread_min_dt[tid] = std::min(pthread_min_dt[tid],dt_2);
        }
      }

// if grid is 3D, compute minimum of (v3 +/- C)

      if (pmb->block_size.nx3 > 1) {
        for (int i=is; i<=ie; ++i){
          Real& dt_3 = dt3(i);
          pthread_min_dt[tid] = std::min(pthread_min_dt[tid],dt_3);
        }
      }

    }
  }

} // end of omp parallel region

// compute minimum across all threads
  Real min_dt = pthread_min_dt[0];
  for (int n=1; n<nthreads; ++n) min_dt = std::min(min_dt,pthread_min_dt[n]);

  delete[] pthread_min_dt;

  pmb->new_block_dt=min_dt;
  return min_dt;
}
