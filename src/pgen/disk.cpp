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

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

// File scope variables
static Real gm0=0.0, r0 = 1.0;
static Real rho0, rho_floor;
static int dflag, vflag;
static Real dslope, pslope, p0_over_r0;
static Real ifield,b0;

// Function Declarations
static Real A3(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A1(const Real x1, const Real x2, const Real x3);

//======================================================================================
//! \file disk.cpp
//  \brief Initializes Keplerian accretion disk in spherical polar coords
//======================================================================================

void Mesh::ProblemGenerator(Hydro *phyd, Field *pfld, ParameterInput *pin)
{
  MeshBlock *pmb = phyd->pmy_block;
  Coordinates *pco = pmb->pcoord;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","gm",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density
  rho0 = pin->GetReal("problem","rho0");
  rho_floor = pin->GetReal("problem","rho_floor"); 
  dflag = pin->GetInteger("problem","dflag");
  vflag = pin->GetInteger("problem","vflag");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure
  if(NON_BAROTROPIC_EOS){
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
  }else{
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }

  // Initialize the magnetic fields

  if (MAGNETIC_FIELDS_ENABLED){

    // Get parameters of inital magnetic fields
    ifield = pin->GetInteger("problem","ifield");
    Real beta = pin->GetReal("problem","beta");
    b0=sqrt(2.*p0_over_r0*rho0/beta);

    // Compute vector potential
    AthenaArray<Real> a1,a2,a3;
    int nx1 = (ie-is)+1 + 2*(NGHOST);
    int nx2 = (je-js)+1 + 2*(NGHOST);
    int nx3 = (ke-ks)+1 + 2*(NGHOST);
    a1.NewAthenaArray(nx3,nx2,nx1);
    a2.NewAthenaArray(nx3,nx2,nx1);
    a3.NewAthenaArray(nx3,nx2,nx1);

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          a1(k,j,i) = A1(pco->x1v(i), pco->x2f(j), pco->x3f(k));
          a2(k,j,i) = A2(pco->x1f(i), pco->x2v(j), pco->x3f(k));
          a3(k,j,i) = A3(pco->x1f(i), pco->x2f(j), pco->x3v(k));
        }
      }
    }

    // Initialize interface fields
    AthenaArray<Real> area,len,len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      if (pmb->block_bcs[inner_x2] == 5) jl=js+1; 
      if (pmb->block_bcs[outer_x2] == 5) ju=je;
      for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->Face2Area(k,j,is,ie,area);
        pmb->pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfld->b.x2f(k,j,i) = -1.0*(len(i+1)*a3(k,j,i+1) - len(i)*a3(k,j,i))/area(i);
        }
        pmb->pcoord->Face2Area(k,j,is,ie,area);
        pmb->pcoord->Edge1Length(k  ,j,is,ie,len);
        pmb->pcoord->Edge1Length(k+1,j,is,ie,len_p1);
        for (int i=is; i<=ie; ++i) {
          pfld->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face3Area(k,j,is,ie,area);
      pmb->pcoord->Edge2Length(k,j,is,ie+1,len);
      for (int i=is; i<=ie; ++i) {
        pfld->b.x3f(k,j,i) = (len(i+1)*a2(k,j,i+1) - len(i)*a2(k,j,i))/area(i);
      }
    }}

    if (pmb->block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,area);
        pmb->pcoord->Edge3Length(k,j  ,is,ie+1,len);
        pmb->pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
        for (int i=is; i<=ie+1; ++i) {
          pfld->b.x1f(k,j,i) = (len_p1(i)*a3(k,j+1,i) - len(i)*a3(k,j,i))/area(i);
        }
      }}

      for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face3Area(k,j,is,ie,area);
        pmb->pcoord->Edge1Length(k,j  ,is,ie,len);
        pmb->pcoord->Edge1Length(k,j+1,is,ie,len_p1);
        for (int i=is; i<=ie; ++i) {
          pfld->b.x3f(k,j,i) -= (len_p1(i)*a1(k,j+1,i) - len(i)*a1(k,j,i))/area(i);
        }
      }}
    }

    if (pmb->block_size.nx3 > 1) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,area);
        pmb->pcoord->Edge2Length(k  ,j,is,ie+1,len);
        pmb->pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
        for (int i=is; i<=ie+1; ++i) {
          pfld->b.x1f(k,j,i) -= (len_p1(i)*a2(k+1,j,i) - len(i)*a2(k,j,i))/area(i);
        }
      }}

      for (int k=ks; k<=ke; ++k) {
        // reset loop limits for polar boundary
        int jl=js; int ju=je+1;
        if (pmb->block_bcs[inner_x2] == 5) jl=js+1; 
        if (pmb->block_bcs[outer_x2] == 5) ju=je;
        for (int j=jl; j<=ju; ++j) {
          pmb->pcoord->Face2Area(k,j,is,ie,area);
          pmb->pcoord->Edge1Length(k  ,j,is,ie,len);
          pmb->pcoord->Edge1Length(k+1,j,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfld->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
  

    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
  }

  //  Initialize density
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (dflag == 1) {
        for (int i=is; i<=ie; ++i) {
          phyd->u(IDN,k,j,i) = std::max(rho0, rho_floor);
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(j);
          Real r = std::max(fabs(x1*sin(x2)),pmb->pmy_mesh->mesh_size.x1min);
          Real z = fabs(x1*cos(x2));
          Real p_over_r = p0_over_r0;
          if (NON_BAROTROPIC_EOS) p_over_r = p0_over_r0*pow(r/r0, pslope);
          Real den = rho0*pow(r/r0,dslope);
          den = den*exp(gm0/p_over_r*(1./sqrt(SQR(r)+SQR(z))-1./r));
          phyd->u(IDN,k,j,i) = std::max(den, rho_floor);
        }
      }
    }
  }

  //  Initialize velocity
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (vflag == 1) {
        for (int i=is; i<=ie; ++i) {
          phyd->u(IM1,k,j,i) = 0.0;
          phyd->u(IM2,k,j,i) = 0.0;
          phyd->u(IM3,k,j,i) = 0.0;
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(j);
          Real r = std::max(fabs(x1*sin(x2)),pmb->pmy_mesh->mesh_size.x1min);
          Real z = fabs(x1*cos(x2));
          Real vel;
          if (phyd->u(IDN,k,j,i) == rho_floor) {
            vel = sqrt(gm0*SQR(r)/(SQR(r)+SQR(z))/sqrt(SQR(r)+SQR(z)));
          } else {
            if (NON_BAROTROPIC_EOS){
              Real p_over_r = p0_over_r0*pow(r/r0, pslope);
              vel = (dslope+pslope)*p_over_r/(gm0/r) + (1.+pslope) - pslope*r/x1;
              vel = sqrt(gm0/r)*sqrt(vel);
            } else {
              vel = dslope*p0_over_r0/(gm0/r)+1.0;
              vel = sqrt(gm0/r)*sqrt(vel);
            }
          }
          phyd->u(IM1,k,j,i) = 0.0;
          phyd->u(IM2,k,j,i) = 0.0;
          phyd->u(IM3,k,j,i) = vel*phyd->u(IDN,k,j,i);
        }
      }
    }
  }

  //  Initialize pressure
  if (NON_BAROTROPIC_EOS){
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(j);
          Real r = std::max(fabs(x1*sin(x2)),pmb->pmy_mesh->mesh_size.x1min);
          Real p_over_r = p0_over_r0*pow(r/r0, pslope);
          Real gamma = phyd->pf_eos->GetGamma();
          phyd->u(IEN,k,j,i) = p_over_r*phyd->u(IDN,k,j,i)/(gamma - 1.0);
          phyd->u(IEN,k,j,i) += 0.5*SQR(phyd->u(IM3,k,j,i))/phyd->u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED){
            phyd->u(IEN,k,j,i) += 
              0.5*(SQR(0.5*(pfld->b.x1f(k,j,i+1) + pfld->b.x1f(k,j,i)))
                 + SQR(0.5*(pfld->b.x2f(k,j+1,i) + pfld->b.x2f(k,j,i)))
                 + SQR(0.5*(pfld->b.x3f(k+1,j,i) + pfld->b.x3f(k,j,i))));
          }
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A3(const Real x1,const Real x2,const Real x3)
//  \brief A3: 3-component of vector potential

static Real A3(const Real x1, const Real x2, const Real x3){
  Real a3=0.0;
  if(ifield==1) {
    a3 = fabs(x1*sin(x2))*b0/2.0;
  }
  return(a3);
}

//--------------------------------------------------------------------------------------
//! \fn static Real A2(const Real x1,const Real x2,const Real x3)
//  \brief A2: 2-component of vector potential

static Real A2(const Real x1, const Real x2, const Real x3){
  Real a2=0.0;
  Real az=0.0;
  if(ifield==2) {
    Real x=x1*sin(x2)*cos(x3);
    Real y=x1*sin(x2)*sin(x3);
    Real z=x1*cos(x2);
    if(sqrt(SQR(x-r0)+SQR(y))<=0.3 && fabs(z)<0.1){
      az=1.0e-3*(0.3-sqrt(SQR(x-r0)+SQR(y)));
    }
    a2=-az*sin(x2);
  }
  return(a2);
}

//--------------------------------------------------------------------------------------
//! \fn static Real A1(const Real x1,const Real x2,const Real x3)
//  \brief A1: 1-component of vector potential

static Real A1(const Real x1, const Real x2, const Real x3){
  Real a1=0.0;
  Real az=0.0;
  if(ifield==2) {
    Real x=x1*sin(x2)*cos(x3);
    Real y=x1*sin(x2)*sin(x3);
    Real z=x1*cos(x2);
    if(sqrt(SQR(x-r0)+SQR(y))<=0.3 && fabs(z)<0.1){
      az=1.e-6*(0.3-sqrt(SQR(x-r0)+SQR(y)));
    }
    a1=az*cos(x2);
  }
  return(a1);
}
