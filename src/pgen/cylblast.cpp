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
//! \file cylblast.cpp
//  \brief Problem generator for the blast wave problem in cylindrical coordinates.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string, c_str()

#if MAGNETIC_FIELDS_ENABLED
#error "This test cannot be run with magnetic fields."
#endif
#ifndef NON_BAROTROPIC_EOS
#error "This test requires non-isothermal EOS."
#endif

#define CUBE(x) ((x)*(x)*(x))

static Real xo, yo, zo, ro, po;

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  const Real da=1.0, pa=1.0, r0=0.25, pb=100.0;
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  if((block_size.nx2 ==1 || block_size.nx3==1) ||
    COORDINATE_SYSTEM!="cylindrical") {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "This test must be run in 3D cylindrical coordinates."
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(Globals::nranks!=1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "This test must be run with one process. OpenMP is OK."
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  for (int k=ks-NGHOST; k<=ke+NGHOST; k++) {
    for (int j=js-NGHOST; j<=je+NGHOST; j++) {
      for (int i=is-NGHOST; i<=ie+NGHOST; i++) {
        phydro->u(IDN,k,j,i) = da;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pa/gm1;
      }
    }
  }

  ro = 2.0, po=0.0, zo=0.0;
  xo=ro*std::cos(po);
  yo=ro*std::sin(po);
  int nthreads = pmy_mesh->GetNumMeshThreads();
#pragma omp parallel for num_threads(nthreads)
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real cvol=pcoord->GetCellVolume(k,j,i);
        for(int kk=0;kk<10;kk++) {
          Real zm=pcoord->x3f(k)+kk*0.1*pcoord->dx3f(k);
          Real zp=pcoord->x3f(k)+(kk+1)*0.1*pcoord->dx3f(k);
          Real zt=(zm+zp)*0.5;
          for(int jj=0;jj<10;jj++) {
            Real pm=pcoord->x2f(j)+jj*0.1*pcoord->dx2f(j);
            Real pp=pcoord->x2f(j)+(jj+1)*0.1*pcoord->dx2f(j);
            Real p=(pm+pp)*0.5;
            for(int ii=0;ii<10;ii++) {
              Real rm=pcoord->x1f(i)+ii*0.1*pcoord->dx1f(i);
              Real rp=pcoord->x1f(i)+(ii+1)*0.1*pcoord->dx1f(i);
              Real r=((CUBE(rp)-CUBE(rm))/3.0)/((SQR(rp)-SQR(rm))/2.0);
              Real xt=r*std::cos(p);
              Real yt=r*std::sin(p);
              Real rad = std::sqrt(SQR(xt-xo)+SQR(yt-yo)+SQR(zt-zo));
              Real fvol=(1.0/2.0)*(SQR(rp) - SQR(rm))*(zp-zm)*(pp-pm);
              Real frac=fvol/cvol;
              if(rad<r0)
                phydro->u(IEN,k,j,i) += frac*(pb-pa)/gm1;
            }
          }
        }
      }
    }
  }
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  // analysis - check deformation of the spherical blast wave by sampling
  int is=pblock->is, ie=pblock->ie,
      js=pblock->js, je=pblock->je,
      ks=pblock->ks, ke=pblock->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pblock->phydro->w,4,IPR,1);

  // find the center
  int ic, jc, kc;
  for(ic=is; ic<=ie; ic++)
    if(pblock->pcoord->x1f(ic)>ro) break;
  ic--;
  for(jc=pblock->js; jc<=pblock->je; jc++)
    if(pblock->pcoord->x2f(jc)>po) break;
  jc--;
  for(kc=pblock->ks; kc<=pblock->ke; kc++)
    if(pblock->pcoord->x3f(kc)>zo) break;
  kc--;

  // find the lowest resolution
  Real drmax=0.0, dpmax=0.0, dzmax=0.0;
  for(int i=is; i<=ie; i++)
    drmax=std::max(drmax,pblock->pcoord->dx1f(i));
  for(int j=js; j<=je; j++)
    dpmax=std::max(dpmax,pblock->pcoord->dx2f(j));
  for(int k=ks; k<=ke; k++)
    dzmax=std::max(dzmax,pblock->pcoord->dx3f(k));
  drmax=std::max(drmax,std::max(3.0*dpmax,dzmax));

  // search pressure maximum in each direction
  Real rmax=0.0, rmin=100.0, rave=0.0;
  int nr=0;
  for(int o=0; o<=6; o++) {
    int ios=0, jos=0, kos=0;
    if(o==1) ios=-10;
    else if(o==2) ios= 10;
    else if(o==3) jos=-10;
    else if(o==4) jos= 10;
    else if(o==5) kos=-10;
    else if(o==6) kos= 10;
    for(int d=0; d<6; d++) {
      Real rm, tm, pm, xm, ym, zm, pmax=0.0;
      int imax, jmax, kmax;
      if(d==0) {
        if(ios!=0) continue;
        jmax=jc+jos, kmax=kc+kos;
        for(int i=ic; i>=is; i--) {
          if(pr(kmax,jmax,i)>pmax) {
            pmax=pr(kmax,jmax,i);
            imax=i;
          }
        }
      }
      else if(d==1) {
        if(ios!=0) continue;
        jmax=jc+jos, kmax=kc+kos;
        for(int i=ic; i<=ie; i++) {
          if(pr(kmax,jmax,i)>pmax) {
            pmax=pr(kmax,jmax,i);
            imax=i;
          }
        }
      }
      else if(d==2) {
        if(jos!=0) continue;
        imax=ic+ios, kmax=kc+kos;
        for(int j=jc; j>=js; j--) {
          if(pr(kmax,j,imax)>pmax) {
            pmax=pr(kmax,j,imax);
            jmax=j;
          }
        }
      }
      else if(d==3) {
        if(jos!=0) continue;
        imax=ic+ios, kmax=kc+kos;
        for(int j=jc; j<=je; j++) {
          if(pr(kmax,j,imax)>pmax) {
            pmax=pr(kmax,j,imax);
            jmax=j;
          }
        }
      }
      else if(d==4) {
        if(kos!=0) continue;
        imax=ic+ios, jmax=jc+jos;
        for(int k=kc; k>=ks; k--) {
          if(pr(k,jmax,imax)>pmax) {
            pmax=pr(k,jmax,imax);
            kmax=k;
          }
        }
      }
      else if(d==5) {
        if(kos!=0) continue;
        imax=ic+ios, jmax=jc+jos;
        for(int k=kc; k<=ke; k++) {
          if(pr(k,jmax,imax)>pmax) {
            pmax=pr(k,jmax,imax);
            kmax=k;
          }
        }
      }
      rm=pblock->pcoord->x1v(imax);
      pm=pblock->pcoord->x2v(jmax);
      zm=pblock->pcoord->x3v(kmax);
      xm=rm*std::cos(pm);
      ym=rm*std::sin(pm);
      Real rad = std::sqrt(SQR(xm-xo)+SQR(ym-yo)+SQR(zm-zo));
      if(rad>rmax) rmax=rad;
      if(rad<rmin) rmin=rad;
      rave+=rad;
      nr++;
    }
  }
  rave/=(Real)nr;
  Real deform=(rmax-rmin)/drmax;
  std::cout << "Offset blast wave test in cylindrical coordinates:" << std::endl 
            << "rmax = " << rmax << ", rmin = " << rmin << ", rave = " << rave << std::endl
            << "deformation = " <<  deform << std::endl;
  if(deform < 1.0)
    std::cout << "\"It looks like spherical,\" said Athena." << std::endl; 
  else if(deform < 2.0)
    std::cout << "\"It looks like more or less spherical,\" said Athena, clearly not being satisfied." << std::endl; 
  else
    std::cout << "You may lose the grace of Athena." << std::endl; 
  pr.DeleteAthenaArray();
  return;
}

