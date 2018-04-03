//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_blockdt.cpp
//  \brief computes timestep using CFL condition on a MEshBlock

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
//[diffusion
#include "diffusion/diffusion.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
//diffusion]

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
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

  AthenaArray<Real> dt1, dt2, dt3;
  dt1.InitWithShallowCopy(dt1_);
  dt2.InitWithShallowCopy(dt2_);
  dt3.InitWithShallowCopy(dt3_);
  Real wi[(NWAVE)];

  Real min_dt = (FLT_MAX);

  for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je; ++j){
      pmb->pcoord->CenterWidth1(k,j,is,ie,dt1);
      pmb->pcoord->CenterWidth2(k,j,is,ie,dt2);
      pmb->pcoord->CenterWidth3(k,j,is,ie,dt3);
      if(!RELATIVISTIC_DYNAMICS) {
#pragma ivdep
        for (int i=is; i<=ie; ++i){
          wi[IDN]=w(IDN,k,j,i);
          wi[IVX]=w(IVX,k,j,i);
          wi[IVY]=w(IVY,k,j,i);
          wi[IVZ]=w(IVZ,k,j,i);
          if (NON_BAROTROPIC_EOS) wi[IPR]=w(IPR,k,j,i);

          if (MAGNETIC_FIELDS_ENABLED) {

            Real bx = bcc(IB1,k,j,i) + fabs(b_x1f(k,j,i)-bcc(IB1,k,j,i));
            wi[IBY] = bcc(IB2,k,j,i);
            wi[IBZ] = bcc(IB3,k,j,i);
            Real cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
            dt1(i) /= (fabs(wi[IVX]) + cf);

            wi[IBY] = bcc(IB3,k,j,i);
            wi[IBZ] = bcc(IB1,k,j,i);
            bx = bcc(IB2,k,j,i) + fabs(b_x2f(k,j,i)-bcc(IB2,k,j,i));
            cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
            dt2(i) /= (fabs(wi[IVY]) + cf);

            wi[IBY] = bcc(IB1,k,j,i);
            wi[IBZ] = bcc(IB2,k,j,i);
            bx = bcc(IB3,k,j,i) + fabs(b_x3f(k,j,i)-bcc(IB3,k,j,i));
            cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
            dt3(i) /= (fabs(wi[IVZ]) + cf);

          } else {

            Real cs = pmb->peos->SoundSpeed(wi);
            dt1(i) /= (fabs(wi[IVX]) + cs);
            dt2(i) /= (fabs(wi[IVY]) + cs);
            dt3(i) /= (fabs(wi[IVZ]) + cs);

          }
          //[diffusion
          if (pmb->phydro->pdif->hydro_diffusion_defined) {
            dt1(i)=std::min(pmb->phydro->pdif->NewDtDiff(w,pmb->pcoord->dx1f(i),k,j,i),dt1(i));
            dt2(i)=std::min(pmb->phydro->pdif->NewDtDiff(w,pmb->pcoord->dx2f(j)*pmb->pcoord->h2v(i),k,j,i),dt2(i));
            dt3(i)=std::min(pmb->phydro->pdif->NewDtDiff(w,pmb->pcoord->dx3f(k)*pmb->pcoord->h31v(i)*fabs(pmb->pcoord->h32v(j)),k,j,i),dt3(i));
          }
          if (MAGNETIC_FIELDS_ENABLED && pmb->pfield->pdif->field_diffusion_defined) {
            dt1(i)=std::min(pmb->pfield->pdif->NewDtFldDiff(pmb->pcoord->dx1f(i),k,j,i),dt1(i));
            dt2(i)=std::min(pmb->pfield->pdif->NewDtFldDiff(pmb->pcoord->dx2f(j)*pmb->pcoord->h2v(i),k,j,i),dt2(i));
            dt3(i)=std::min(pmb->pfield->pdif->NewDtFldDiff(pmb->pcoord->dx3f(k)*pmb->pcoord->h31v(i)*fabs(pmb->pcoord->h32v(j)),k,j,i),dt3(i));
          }
          //diffusion]
        }
      }

      // compute minimum of (v1 +/- C)
      for (int i=is; i<=ie; ++i){
        Real& dt_1 = dt1(i);
        min_dt = std::min(min_dt,dt_1);
      }

      // if grid is 2D/3D, compute minimum of (v2 +/- C)
      if (pmb->block_size.nx2 > 1) {
        for (int i=is; i<=ie; ++i){
          Real& dt_2 = dt2(i);
          min_dt = std::min(min_dt,dt_2);
        }
      }

      // if grid is 3D, compute minimum of (v3 +/- C)
      if (pmb->block_size.nx3 > 1) {
        for (int i=is; i<=ie; ++i){
          Real& dt_3 = dt3(i);
          min_dt = std::min(min_dt,dt_3);
        }
      }

    }
  }

  min_dt *= pmb->pmy_mesh->cfl_number;

  if(UserTimeStep_!=NULL) {
    min_dt = std::min(min_dt, UserTimeStep_(pmb));
  }

  pmb->new_block_dt=min_dt;
  return min_dt;
}
