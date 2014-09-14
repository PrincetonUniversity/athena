//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "fluid.hpp"

// C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena headers
#include "../athena.hpp"                // array access, macros, Real
#include "../athena_arrays.hpp"         // AthenaArray
#include "bvals/bvals.hpp"              // FluidBoundaryConditions
#include "eos/eos.hpp"                  // FluidEqnOfState
#include "srcterms/srcterms.hpp"        // FluidSourceTerms
#include "integrators/integrators.hpp"  // FluidIntegrator
#include "../mesh.hpp"                  // MeshBlock, Mesh
#include "../coordinates/coordinates.hpp" // VolumeCenterWidth()

#include <omp.h>

//======================================================================================
//! \file fluid.cpp
//  \brief implementation of functions in class Fluid
//======================================================================================

// constructor, initializes data structures and parameters

Fluid::Fluid(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;

// Allocate memory for primitive/conserved variables

  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  u.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

// Allocate memory for primitive/conserved variables at intermediate-time step

  u1.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w1.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

  // Allocate memory for metric
  // TODO: this should only be done if we are in GR
  g.NewAthenaArray(NMETRIC, ncells1);
  g_inv.NewAthenaArray(NMETRIC, ncells1);

// Allocate memory for scratch arrays

  dt1_.NewAthenaArray(ATHENA_MAX_NUM_THREADS,ncells1);
  dt2_.NewAthenaArray(ATHENA_MAX_NUM_THREADS,ncells1);
  dt3_.NewAthenaArray(ATHENA_MAX_NUM_THREADS,ncells1);

// Construct ptrs to objects of various classes needed to integrate fluid eqns 

  pf_integrator = new FluidIntegrator(this);
  pf_bcs = new FluidBoundaryConditions(this,pin);
  pf_eos = new FluidEqnOfState(this,pin);
  pf_srcterms = new FluidSourceTerms(this,pin);
}

// destructor

Fluid::~Fluid()
{
  pmy_block = NULL; // MeshBlock destructor will free this memory
  u.DeleteAthenaArray();
  w.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  w1.DeleteAthenaArray();
  g.DeleteAthenaArray();
  g_inv.DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();

  delete pf_integrator;
  delete pf_bcs;
  delete pf_eos;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Fluid::NewTimeStep(MeshBlock *pmb)
{
  int tid=0;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real tmin_dt[ATHENA_MAX_NUM_THREADS];

  AthenaArray<Real> w = pmb->pfluid->w.ShallowCopy();

  for (int i=0; i<(ATHENA_MAX_NUM_THREADS); ++i) tmin_dt[i] = (FLT_MAX);

#pragma omp parallel default(shared) private(tid) num_threads(ATHENA_MAX_NUM_THREADS)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> *pdt1 = dt1_.ShallowSlice(tid,1);
  AthenaArray<Real> *pdt2 = dt2_.ShallowSlice(tid,1);
  AthenaArray<Real> *pdt3 = dt3_.ShallowSlice(tid,1);
  Real wi[NVAR];

  for (int k=ks; k<=ke; ++k){

#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j){
      Real& dx2 = pmb->dx2f(j);
      Real& dx3 = pmb->dx3f(k);
#pragma simd
      for (int i=is; i<=ie; ++i){
        wi[IDN]=w(IDN,k,j,i);
        wi[IVX]=w(IVX,k,j,i);
        wi[IVY]=w(IVY,k,j,i);
        wi[IVZ]=w(IVZ,k,j,i);
        if (NON_BAROTROPIC_EOS) wi[IEN]=w(IEN,k,j,i);
        Real& dx1  = pmb->dx1f(i);
        Real& dt_1 = (*pdt1)(i);
        Real& dt_2 = (*pdt2)(i);
        Real& dt_3 = (*pdt3)(i);

        if (RELATIVISTIC_DYNAMICS) {
          dt_1 = dx1;
          dt_2 = dx2;
          dt_3 = dx3;
        } else {
          Real cs = pf_eos->SoundSpeed(wi);
          dt_1 = pmy_block->pcoord->VolumeCenterWidth1(k,j,i)/(fabs(wi[IVX]) + cs);
          dt_2 = pmy_block->pcoord->VolumeCenterWidth2(k,j,i)/(fabs(wi[IVY]) + cs);
          dt_3 = pmy_block->pcoord->VolumeCenterWidth3(k,j,i)/(fabs(wi[IVZ]) + cs);
        }
      }

// compute minimum of (v1 +/- C)

      for (int i=is; i<=ie; ++i){
        Real& dt_1 = (*pdt1)(i);
        tmin_dt[tid] = std::min(tmin_dt[tid],dt_1);
      }
    
// if grid is 2D/3D, compute minimum of (v2 +/- C)

      if (pmb->block_size.nx2 > 1) {
        for (int i=is; i<=ie; ++i){
          Real& dt_2 = (*pdt2)(i);
          tmin_dt[tid] = std::min(tmin_dt[tid],dt_2);
        }
      }

// if grid is 3D, compute minimum of (v3 +/- C)

      if (pmb->block_size.nx3 > 1) {
        for (int i=is; i<=ie; ++i){
          Real& dt_3 = (*pdt3)(i);
          tmin_dt[tid] = std::min(tmin_dt[tid],dt_3);
        }
      }

    }
  }
} // end of omp parallel region

// compute minimum across all threads
  Real min_dt = tmin_dt[0];
  for (int i=1; i<(ATHENA_MAX_NUM_THREADS); ++i) min_dt = std::min(min_dt,tmin_dt[i]);

// compute new global timestep
  Mesh *pm = pmb->pmy_domain->pmy_mesh;
  Real old_dt = pm->dt;
  pm->dt = std::min( ((pm->cfl_number)*min_dt) , (2.0*old_dt) );

  return;
}
