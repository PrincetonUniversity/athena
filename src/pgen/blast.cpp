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
//! \file blast.cpp
//  \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//         cylindrical, and spherical coordinates.  Contains post-processing code
//         to check whether blast is spherical for regression tests
//
// REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for 
//   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.
//======================================================================================

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real rout = pin->GetReal("problem","radius");
  Real rin  = rout - pin->GetOrAddReal("problem","ramp",0.0);
  Real pa   = pin->GetOrAddReal("problem","pamb",1.0);
  Real da   = pin->GetOrAddReal("problem","damb",1.0);
  Real prat = pin->GetReal("problem","prat");
  Real drat = pin->GetOrAddReal("problem","drat",1.0);
  Real b0,theta;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem","b0");
    theta = (PI/180.0)*pin->GetReal("problem","angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Real x2_0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Real x3_0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Real x0,y0,z0;
  if (COORDINATE_SYSTEM == "cartesian") {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  }

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rad;
    if (COORDINATE_SYSTEM == "cartesian") {
      Real x = pcoord->x1v(i);
      Real y = pcoord->x2v(j);
      Real z = pcoord->x3v(k);
      rad = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else if (COORDINATE_SYSTEM == "cylindrical") {
      Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
      Real z = pcoord->x3v(k);
      rad = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else if (COORDINATE_SYSTEM == "spherical_polar") {
      Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
      Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
      Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      rad = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    }
    Real den = da;
    if (rad < rout) {
      if (rad < rin) {
        den = drat*da;
      } else {   // add smooth ramp in density
        Real f = (rad-rin) / (rout-rin);
        Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
        den = std::exp(log_den);
      }
    }

    phydro->u(IDN,k,j,i) = den;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
    if (NON_BAROTROPIC_EOS) {
      Real pres = pa;
      if (rad < rout) {
        if (rad < rin) {
          pres = prat*pa;
        } else {  // add smooth ramp in pressure
          Real f = (rad-rin) / (rout-rin);
          Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
          pres = std::exp(log_pres);
        }
      }
      phydro->u(IEN,k,j,i) = pres/gm1;
      if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
        phydro->u(IEN,k,j,i) += den;
    }
  }}}

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = b0*cos(theta);
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = b0*sin(theta);
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      phydro->u(IEN,k,j,i) += 0.5*b0*b0;
    }}}
  }

}

//======================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//======================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  // analysis - check shape of the spherical blast wave
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pblock->phydro->w,4,IPR,1);

  // get coordinate location of the center, convert to Cartesian
  Real x1_0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Real x2_0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Real x3_0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Real x0,y0,z0;
  if (COORDINATE_SYSTEM == "cartesian") {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  }

  // find indices of the center
  int ic, jc, kc;
  for(ic=is; ic<=ie; ic++)
    if(pblock->pcoord->x1f(ic) > x1_0) break;
  ic--;
  for(jc=pblock->js; jc<=pblock->je; jc++)
    if(pblock->pcoord->x2f(jc) > x2_0) break;
  jc--;
  for(kc=pblock->ks; kc<=pblock->ke; kc++)
    if(pblock->pcoord->x3f(kc) > x3_0) break;
  kc--;

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
      Real pmax=0.0;
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

      Real xm, ym, zm;
      Real x1m=pblock->pcoord->x1v(imax);
      Real x2m=pblock->pcoord->x2v(jmax);
      Real x3m=pblock->pcoord->x3v(kmax);
      if (COORDINATE_SYSTEM == "cartesian") {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      } else if (COORDINATE_SYSTEM == "cylindrical") {
        xm = x1m*std::cos(x2m);
        ym = x1m*std::sin(x2m);
        zm = x3m;
      } else if (COORDINATE_SYSTEM == "spherical_polar") {
        xm = x1m*std::sin(x2m)*std::cos(x3m);
        ym = x1m*std::sin(x2m)*std::sin(x3m);
        zm = x1m*std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm-x0)+SQR(ym-y0)+SQR(zm-z0));
      if(rad>rmax) rmax=rad;
      if(rad<rmin) rmin=rad;
      rave+=rad;
      nr++;
    }
  }
  rave/=(Real)nr;

  // use physical grid spacing at center of blast
  Real dr_max;
  Real  x1c = pblock->pcoord->x1v(ic);
  Real dx1c = pblock->pcoord->dx1f(ic);
  Real  x2c = pblock->pcoord->x2v(jc);
  Real dx2c = pblock->pcoord->dx2f(jc);
  Real dx3c = pblock->pcoord->dx3f(kc);
  if (COORDINATE_SYSTEM == "cartesian") {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), dx3c);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), x1c*std::sin(x2c)*dx3c);
  }
  Real deform=(rmax-rmin)/dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if((pfile = fopen(fname.c_str(),"r")) != NULL){
      if((pfile = freopen(fname.c_str(),"a",pfile)) == NULL){
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

    // The file does not exist -- open the file in write mode and add headers
    } else {
      if((pfile = fopen(fname.c_str(),"w")) == NULL){
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    fprintf(pfile,"# Offset blast wave test in %s coordinates:\n",COORDINATE_SYSTEM);
    fprintf(pfile,"# Rmax       Rmin       Rave        Deformation\n");
    fprintf(pfile,"%e  %e  %e  %e \n",rmax,rmin,rave,deform);
    fclose(pfile);
  }

  pr.DeleteAthenaArray();
  return;
}
