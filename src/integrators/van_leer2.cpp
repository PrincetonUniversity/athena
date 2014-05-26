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
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdio.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"
#include "../geometry/geometry.hpp"
#include "integrators.hpp"

//using namespace COORDINATE_SYSTEM;
//using namespace RIEMANN_SOLVER;
//using namespace RECONSTRUCTION_ALGORITHM;

#define SUM_ON 0

//======================================================================================
/*! \file van_leer2.cpp
 *  \brief van-Leer (MUSCL-Hancock) second-order integrator
 *====================================================================================*/

//--------------------------------------------------------------------------------------
///*! \fn  void FluidIntegrator::PredictVanLeer2
// *  \brief predictor step for 2nd order VL integrator */

void FluidIntegrator::Predict(Block *pb)
{
  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;
  Real dt = pb->pparent_domain->pparent_mesh->dt;

  Real sum=0.0; Real sum_2=0.0;
  int ndata = (pb->block_size.nx1)*(pb->block_size.nx2)*(pb->block_size.nx3)*(NVAR);
 
  AthenaArray<Real> u = pb->pfluid->u.ShallowCopy();
  AthenaArray<Real> w = pb->pfluid->w.ShallowCopy();
  AthenaArray<Real> u1 = pb->pfluid->u1.ShallowCopy();
  AthenaArray<Real> w1 = pb->pfluid->w1.ShallowCopy();

  AthenaArray<Real> wl = wl_.ShallowCopy();
  AthenaArray<Real> wr = wr_.ShallowCopy();
  AthenaArray<Real> flx = flx_.ShallowCopy();
 
  AthenaArray<Real> area = pb->pgeometry->face_area.ShallowCopy();
  AthenaArray<Real> vol  = pb->pgeometry->cell_volume.ShallowCopy();

#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
    sum += pb->pfluid->u(i);
    sum_2 += pb->pfluid->w(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif

//--------------------------------------------------------------------------------------
// i-direction 

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    ReconstructionFunc(k,j,is,ie+1,1,w,wl,wr);

    RiemannSolver(is,ie+1,wl,wr,flx);

    pb->pgeometry->Area1Face(k,j,is,ie+1,area);
    pb->pgeometry->CellVolume(k,j,is,ie,vol);

    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& ui  = u (n,k,j,i);
        Real& u1i = u1(n,k,j,i);
        Real& flxi   = flx(n,  i);
        Real& flxip1 = flx(n,i+1);

        Real& area_i   = area(i);
        Real& area_ip1 = area(i+1);
        Real& dvol = vol(i);
 
        u1i = ui - 0.5*dt*(area_ip1*flxip1 - area_i*flxi)/dvol;
      }
    }

  }}

#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
      sum += pb->pfluid->u(i);
      sum_2 += u1(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif


//--------------------------------------------------------------------------------------
// j-direction

  if (pb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je+1; ++j){

      ReconstructionFunc(k,j,is,ie,2,w,wl,wr);

      RiemannSolver(is,ie,wl,wr,flx); 

      pb->pgeometry->Area2Face(k,j,is,ie,area);

      if (j>js) {
        pb->pgeometry->CellVolume(k,j-1,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1jm1 = u1(n,k,j-1,i);
            Real& flxj  = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);
  
            u1jm1 -= 0.5*dt*area_i*flxj/dvol;
          }
        }
      }

      if (j<(je+1)) {
        pb->pgeometry->CellVolume(k,j,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1j   = u1(n,k,j  ,i);
            Real& flxj  = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);

            u1j   += 0.5*dt*area_i*flxj/dvol;
          }
        }
      }

    }}
  }


#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
    sum += pb->pfluid->u(i);
    sum_2 += u1(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif



//--------------------------------------------------------------------------------------
// k-direction 

  if (pb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke+1; ++k){
    for (int j=js; j<=je; ++j){

      ReconstructionFunc(k,j,is,ie,3,w,wl,wr);

      RiemannSolver(is,ie,wl,wr,flx);

      pb->pgeometry->Area3Face(k,j,is,ie,area);

      if (k>ks) {
        pb->pgeometry->CellVolume(k-1,j,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1km1 = u1(n,k-1,j,i);
            Real& flxk = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);

            u1km1 -= 0.5*dt*area_i*flxk/dvol;
          }
        }
      }

      if (k<(ke+1)) {
        pb->pgeometry->CellVolume(k,j,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& u1k   = u1(n,k  ,j,i);
            Real& flxk = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);

            u1k += 0.5*dt*area_i*flxk/dvol;
          }
        }
      }

    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void FluidIntegrator::Correct(Block *pb)
{
  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;
  Real dt = pb->pparent_domain->pparent_mesh->dt;

  Real sum=0.0; Real sum_2=0.0;
  int ndata = (pb->block_size.nx1)*(pb->block_size.nx2)*(pb->block_size.nx3)*(NVAR);
 
  AthenaArray<Real> u = pb->pfluid->u.ShallowCopy();
  AthenaArray<Real> w = pb->pfluid->w.ShallowCopy();
  AthenaArray<Real> u1 = pb->pfluid->u1.ShallowCopy();
  AthenaArray<Real> w1 = pb->pfluid->w1.ShallowCopy();

  AthenaArray<Real> wl = wl_.ShallowCopy();
  AthenaArray<Real> wr = wr_.ShallowCopy();
  AthenaArray<Real> flx = flx_.ShallowCopy();

  AthenaArray<Real> area = pb->pgeometry->face_area.ShallowCopy();
  AthenaArray<Real> vol  = pb->pgeometry->cell_volume.ShallowCopy();
 
#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
    sum += pb->pfluid->u(i);
    sum_2 += pb->pfluid->w(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif

//--------------------------------------------------------------------------------------
// i-direction 

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    ReconstructionFunc(k,j,is,ie+1,1,w1,wl,wr);

    RiemannSolver(is,ie+1,wl,wr,flx); 

    pb->pgeometry->Area1Face(k,j,is,ie+1,area);
    pb->pgeometry->CellVolume(k,j,is,ie,vol);

    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& ui  = u (n,k,j,i);
        Real& flxi   = flx(n,  i);
        Real& flxip1 = flx(n,i+1);

        Real& area_i   = area(i);
        Real& area_ip1 = area(i+1);
        Real& dvol = vol(i);
 
        ui -= dt*(area_ip1*flxip1 - area_i*flxi)/dvol;
      }
    }

  }}

#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
      sum += pb->pfluid->u(i);
      sum_2 += u1(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif


//--------------------------------------------------------------------------------------
// j-direction

  if (pb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k){
    for (int j=js; j<=je+1; ++j){

      ReconstructionFunc(k,j,is,ie,2,w1,wl,wr);

      RiemannSolver(is,ie,wl,wr,flx); 

      pb->pgeometry->Area2Face(k,j,is,ie,area);

      if (j>js){
        pb->pgeometry->CellVolume(k,j-1,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& ujm1 = u(n,k,j-1,i);
            Real& flxj  = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);
  
            ujm1 -= dt*area_i*flxj/dvol;
          }
        }
      }

      if (j>(je+1)){
        pb->pgeometry->CellVolume(k,j,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& uj   = u(n,k,j  ,i);
            Real& flxj  = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);
  
            uj += dt*area_i*flxj/dvol;
          }
        }
      }

    }}
  }

#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
    sum += pb->pfluid->u(i);
    sum_2 += u1(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif



//--------------------------------------------------------------------------------------
// k-direction 

  if (pb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke+1; ++k){
    for (int j=js; j<=je; ++j){

      ReconstructionFunc(k,j,is,ie,3,w1,wl,wr);

      RiemannSolver(is,ie,wl,wr,flx);

      pb->pgeometry->Area3Face(is,ie,j,k,area);

      if (k<ks){
        pb->pgeometry->CellVolume(k-1,j,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& ukm1 = u(n,k-1,j,i);
            Real& flxk = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);

            ukm1 -= dt*area_i*flxk/dvol;
          }
        }
      }

      if (k>(ke+1)){
        pb->pgeometry->CellVolume(k,j,is,ie,vol);
        for (int n=0; n<NVAR; ++n){
#pragma simd
          for (int i=is; i<=ie; ++i){
            Real& uk   = u(n,k  ,j,i);
            Real& flxk = flx(n,i);
            Real& area_i   = area(i);
            Real& dvol = vol(i);

            uk   += dt*area_i*flxk/dvol;
          }
        }
      }

    }}
  }

  return;
}
