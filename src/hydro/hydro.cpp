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
//! \file hydro.cpp
//  \brief implementation of functions in class Hydro
//======================================================================================


// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "../athena.hpp"                // array access, macros, Real
#include "../athena_arrays.hpp"         // AthenaArray
#include "eos/eos.hpp"
#include "srcterms/srcterms.hpp"
#include "viscosity/viscosity.hpp"      // Viscosity 
#include "integrators/hydro_integrator.hpp"
#include "../mesh.hpp"                  // MeshBlock, Mesh
#include "../coordinates/coordinates.hpp" // CenterWidth()
#include "../field/field.hpp"             // B-fields

// this class header
#include "hydro.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// constructor, initializes data structures and parameters

Hydro::Hydro(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;

// Allocate memory for primitive/conserved variables

  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  u.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);
  w.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);

// Allocate memory for primitive/conserved variables at intermediate-time step

  u1.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);
  w1.NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1);

  flux[x1face].NewAthenaArray(NHYDRO,ncells3,ncells2,ncells1+1);
  if (pmy_block->block_size.nx2 > 1) 
    flux[x2face].NewAthenaArray(NHYDRO,ncells3,ncells2+1,ncells1);
  if (pmy_block->block_size.nx3 > 1) 
    flux[x3face].NewAthenaArray(NHYDRO,ncells3+1,ncells2,ncells1);

  // Allocate memory for metric
  // TODO: this should only be done if we are in GR
  g.NewAthenaArray(NMETRIC, ncells1);
  g_inv.NewAthenaArray(NMETRIC, ncells1);

// Allocate memory for scratch arrays

  int nthreads = pmy_block->pmy_mesh->GetNumMeshThreads();
  dt1_.NewAthenaArray(nthreads,ncells1);
  dt2_.NewAthenaArray(nthreads,ncells1);
  dt3_.NewAthenaArray(nthreads,ncells1);

// Allocate memory for internal hydro output variables (if needed)

  ifov.NewAthenaArray(NIFOV,ncells3,ncells2,ncells1);

// Construct ptrs to objects of various classes needed to integrate hydro eqns 

  pintegrator = new HydroIntegrator(this,pin);
  peos = new HydroEqnOfState(this,pin);
  pf_srcterms = new HydroSourceTerms(this,pin);
  if(VISCOSITY) pf_viscosity = new Viscosity(this,pin);
}

// destructor

Hydro::~Hydro()
{
  u.DeleteAthenaArray();
  w.DeleteAthenaArray();
  u1.DeleteAthenaArray();
  w1.DeleteAthenaArray();
  g.DeleteAthenaArray();
  g_inv.DeleteAthenaArray();

  flux[x1face].DeleteAthenaArray();
  if (pmy_block->block_size.nx2 > 1) flux[x2face].DeleteAthenaArray();
  if (pmy_block->block_size.nx3 > 1) flux[x3face].DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();

  ifov.DeleteAthenaArray();

  delete pintegrator;
  delete peos;
  delete pf_srcterms;
  if(VISCOSITY) delete pf_viscosity;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

Real Hydro::NewBlockTimeStep(MeshBlock *pmb)
{
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
#pragma simd
      for (int i=is; i<=ie; ++i){
        wi[IDN]=w(IDN,k,j,i);
        wi[IVX]=w(IVX,k,j,i);
        wi[IVY]=w(IVY,k,j,i);
        wi[IVZ]=w(IVZ,k,j,i);
        if (NON_BAROTROPIC_EOS) wi[IEN]=w(IEN,k,j,i);

        if (RELATIVISTIC_DYNAMICS) {

          dt1(i) = pmy_block->pcoord->CenterWidth1(k,j,i);
          dt2(i) = pmy_block->pcoord->CenterWidth2(k,j,i);
          dt3(i) = pmy_block->pcoord->CenterWidth3(k,j,i);

        } else if (MAGNETIC_FIELDS_ENABLED) {

          Real bx = bcc(IB1,k,j,i) + fabs(b_x1f(k,j,i)-bcc(IB1,k,j,i));
          wi[IBY] = bcc(IB2,k,j,i);
          wi[IBZ] = bcc(IB3,k,j,i);
          Real cf = peos->FastMagnetosonicSpeed(wi,bx);
          dt1(i)= pmy_block->pcoord->CenterWidth1(k,j,i)/(fabs(wi[IVX]) + cf);

          wi[IBY] = bcc(IB3,k,j,i);
          wi[IBZ] = bcc(IB1,k,j,i);
          bx = bcc(IB2,k,j,i) + fabs(b_x2f(k,j,i)-bcc(IB2,k,j,i));
          cf = peos->FastMagnetosonicSpeed(wi,bx);
          dt2(i)= pmy_block->pcoord->CenterWidth2(k,j,i)/(fabs(wi[IVY]) + cf);

          wi[IBY] = bcc(IB1,k,j,i);
          wi[IBZ] = bcc(IB2,k,j,i);
          bx = bcc(IB3,k,j,i) + fabs(b_x3f(k,j,i)-bcc(IB3,k,j,i));
          cf = peos->FastMagnetosonicSpeed(wi,bx);
          dt3(i)= pmy_block->pcoord->CenterWidth3(k,j,i)/(fabs(wi[IVZ]) + cf);

        } else {

          Real cs = peos->SoundSpeed(wi);
          dt1(i)= pmy_block->pcoord->CenterWidth1(k,j,i)/(fabs(wi[IVX]) + cs);
          dt2(i)= pmy_block->pcoord->CenterWidth2(k,j,i)/(fabs(wi[IVY]) + cs);
          dt3(i)= pmy_block->pcoord->CenterWidth3(k,j,i)/(fabs(wi[IVZ]) + cs);

        }
        if(VISCOSITY){
          dt1(i)=std::min(pmb->phydro->pf_viscosity->VisDt(pmy_block->pcoord->CenterWidth1(k,j,i),k,j,i),dt1(i));
          dt2(i)=std::min(pmb->phydro->pf_viscosity->VisDt(pmy_block->pcoord->CenterWidth2(k,j,i),k,j,i),dt2(i));
          dt3(i)=std::min(pmb->phydro->pf_viscosity->VisDt(pmy_block->pcoord->CenterWidth3(k,j,i),k,j,i),dt3(i));
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
