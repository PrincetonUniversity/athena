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
//! \file ssheet.c
//  \brief Shearing wave problem generator for 2D/3D problems.
//
//======================================================================================
//
//[JMSHI
//
// C++ headers
#include <iostream>   // cout, endl
#include <stdlib.h>   // exit
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

#if !SHEARING_BOX
#error "This problem generator requires shearing box"
#endif

static Real amp, nwx, nwy; // amplitude, Wavenumbers
static int ipert; // initial pattern
static Real gm1,iso_cs;
static Real x1size,x2size,x3size;
static Real Omega_0,qshear;

void ShearInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   int is, int ie, int js, int je, int ks, int ke);
void ShearOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   int is, int ie, int js, int je, int ks, int ke);

void LinearSlope(const int nvar, const int ny, const int nx, const AthenaArray<Real> &w,
				   AthenaArray<Real> &dw);

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // initialize global variables
  amp = pin->GetReal("problem","amp");
  nwx = pin->GetInteger("problem","nwx");
  nwy = pin->GetInteger("problem","nwy");
  ipert = pin->GetInteger("problem","ipert");
  Omega_0 = pin->GetOrAddReal("problem","Omega0",0.001);
  qshear  = pin->GetOrAddReal("problem","qshear",1.5);

  // enroll boundary value function pointers
  //EnrollUserBoundaryFunction(INNER_X1, ShearInnerX1);
  //EnrollUserBoundaryFunction(OUTER_X1, ShearOuterX1);
  return;
}



//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
//  if (pmy_mesh->mesh_size.nx2 == 1 || pmy_mesh->mesh_size.nx3 > 1) {
//	  std::cout << "[ssheet.cpp]: only works on 2D grid" << std::endl;
//	  exit(0);
//  }

//  if (NON_BAROTROPIC_EOS) {
//	  std::cout << "[ssheet.cpp]: only works for isothermal eos" << std::endl;
//	  exit(0);
//  }

  if (MAGNETIC_FIELDS_ENABLED) {
	  std::cout << "[ssheet.cpp]: only works for hydro alone" << std::endl;
	  exit(0);
  }

  Real d0 = 1.0;
  Real p0 = 1e-6;

  if (NON_BAROTROPIC_EOS) {
    gm1 = (peos->GetGamma() - 1.0);
	iso_cs = sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
	p0 = d0*SQR(iso_cs);
  }
  std::cout << "iso_cs = " << iso_cs << std::endl;
  std::cout << "d0 = " << d0 << std::endl;
  std::cout << "p0 = " << p0 << std::endl;
  std::cout << "ipert  = " << ipert  << std::endl;


  x1size = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  x2size = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  x3size = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;
  std::cout << "[ssheet.cpp]: [Lx,Ly,Lz] = [" <<x1size <<","<<x2size<<","<<x3size<<"]"<<std::endl;

  Real kx = (2.0*PI/x1size)*((double)nwx);
  Real ky = (2.0*PI/x2size)*((double)nwy);

  Real x1,x2,rd,rp,rvx,rvy;
// update the physical variables as initial conditions
  //int nx1 = (ie-is)+1 + 2*(NGHOST);
  //int nx2 = (je-js)+1 + 2*(NGHOST);
  //int nx3 = (ke-ks)+1 + 2*(NGHOST);

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      x1 = pcoord->x1v(i);
	  x2 = pcoord->x2v(j);
	  rd = d0;
	  rp = p0;
	  if (ipert == 0) {
	    // 1) pure shear bg flow:
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) -= rd*(qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = 0.0;
	  } else if (ipert == 1) {
	    // 2) initialize with shwave in velocity
        //rd = d0; //d0*(1.0+amp*cos(kx*x1 + ky*x2));
	    rvx = amp*iso_cs*sin(kx*x1 + ky*x2);
	    rvy = amp*iso_cs*(kx/ky)*sin(kx*x1 + ky*x2);
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) -= rd*(rvy + qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = 0.0;
	  } else if (ipert == 2) {
	    // 3) epicyclic oscillation
		rvx = 0.1*iso_cs;
		rvy = 0.0;
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) -= rd*(rvy + qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = 0.0;
	  } else {
          std::cout << "[ssheet.cpp] ipert = " <<ipert <<" is unrecognized " <<std::endl;
	      exit(0);
	  }
      if (NON_BAROTROPIC_EOS){
		phydro->u(IEN,k,j,i) = rp/gm1 +0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))+
											 SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
	  }
    }
  }}

  // enroll boundary value function pointers
  // in InitUserMeshProperties now.
  //
  //EnrollUserBoundaryFunction(INNER_X1, ShearInnerX1);
  //EnrollUserBoundaryFunction(OUTER_X1, ShearOuterX1);

  return;
}


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================
void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void ShearInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for ssheet problem

void ShearInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   int is, int ie, int js, int je, int ks, int ke)
{
  Real yshear, deltay, frac;
  int jplus, ji;

  //Real qshear = 1.5, Omega_0 = 1.0;
  Real qomL = qshear*Omega_0*x1size;

  yshear = qomL*pmb->pmy_mesh->time;
  deltay = fmod(yshear, x2size);
  jplus  = (int)(deltay/pco->dx2v(js));
  frac   = (fmod(deltay,pco->dx2v(js)))/pco->dx2v(js); // assuming uniform grid azimuthally
//  std::cout << "[ShearInnerX1]: [time, yshear,deltay,j+,frac] = ["<<pmb->pmy_mesh->time<<","<<yshear<<","<<","<<deltay<<","
//	  <<jplus<<","<<frac<<"]"<< std::endl;

  // Initialize boundary value arrays
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nvar = (NHYDRO);  // for now IDN, IVX, IVY and IVZ, and IEN if non_barotropic
  AthenaArray<Real> bval;
  bval.NewAthenaArray(nvar,nx2,(NGHOST));
  AthenaArray<Real> dbval;
  dbval.NewAthenaArray(nvar,nx2,(NGHOST));
  int nyzone = (je-js)+1;

  // set bval variables in inlet ghost zones
  for(int j=0; j<nx2; ++j) {
	ji = fmod(j-(jplus+(NGHOST))+nyzone, nyzone) + (NGHOST);
	if (ji < js) ji += nyzone;
	if (ji > je) ji -= nyzone;
	if ((ji < js) || (ji > je)){
	  std::cout << "[ShearInnerX1]: ji = " << ji << " < js = " << js << "or > je= " << je << std::endl;
	  exit(1);
	}
    for(int i=1; i<=(NGHOST); ++i) {
      bval(IDN,j,i-1) = a(IDN,ks,ji,ie-(NGHOST)+i);
      bval(IVX,j,i-1) = a(IVX,ks,ji,ie-(NGHOST)+i);
      bval(IVY,j,i-1) = a(IVY,ks,ji,ie-(NGHOST)+i);
      bval(IVZ,j,i-1) = a(IVZ,ks,ji,ie-(NGHOST)+i);
      if (NON_BAROTROPIC_EOS) {
	    bval(IEN,j,i-1) = a(IEN,ks,ji,ie-(NGHOST)+i);
	  }
    }
  }
  // skip b-field;

  // zone-center slope dbval for the fractional part
  LinearSlope(nvar,nx2,(NGHOST),bval,dbval);
  // piece-wise linear interpolation of the fractional part

  for(int j=js; j<=je; ++j) {
    for(int i=1; i<=(NGHOST); ++i) {
	  int ib = (NGHOST) - i;
	  a(IDN,ks,j,is-i) = (1.0-frac)*bval(IDN,j,ib) + frac*bval(IDN,j-1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IDN,j-1,ib)-dbval(IDN,j,ib));
	  a(IVX,ks,j,is-i) = (1.0-frac)*bval(IVX,j,ib) + frac*bval(IVX,j-1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IVX,j-1,ib)-dbval(IVX,j,ib));
	  a(IVY,ks,j,is-i) = (1.0-frac)*bval(IVY,j,ib) + frac*bval(IVY,j-1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IVY,j-1,ib)-dbval(IVY,j,ib)) +
					 //qshear*Omega_0*x1size*a(IDN,ks,j,is-i);
					 qshear*Omega_0*x1size;
	  a(IVZ,ks,j,is-i) = (1.0-frac)*bval(IVZ,j,ib) + frac*bval(IVZ,j-1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IVZ,j-1,ib)-dbval(IVZ,j,ib));
      if (NON_BAROTROPIC_EOS) {
	    a(IEN,ks,j,is-i) = (1.0-frac)*bval(IEN,j,ib) + frac*bval(IEN,j-1,ib) +
		               0.5*frac*(1.0-frac)*(dbval(IEN,j-1,ib)-dbval(IEN,j,ib));
	  }
	}}

}

//--------------------------------------------------------------------------------------
//! \fn void ShearOuterX1()
//  \brief Sets boundary condition on right X boundary (oib) for ssheet problem

void ShearOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   int is, int ie, int js, int je, int ks, int ke)
{
  Real yshear, deltay, frac;
  int jplus, jo;

  //Real qshear = 1.5, Omega_0 = 1.0;
  Real qomL = qshear*Omega_0*x1size;

  yshear = qomL*pmb->pmy_mesh->time;
  deltay = fmod(yshear, x2size);
  jplus  = (int)(deltay/pco->dx2v(js));
  frac   = (fmod(deltay,pco->dx2v(js)))/pco->dx2v(js);
//  std::cout << "[ShearOuterX1]: [time, yshear,deltay,j+,frac] = ["<<pmb->pmy_mesh->time<<","<<yshear<<","<<","<<deltay<<","
//	  <<jplus<<","<<frac<<"]"<< std::endl;

  // Initialize boundary value arrays
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nvar = (NHYDRO);  // for now IDN, IVX, IVY and IVZ
  AthenaArray<Real> bval;
  bval.NewAthenaArray(nvar,nx2,(NGHOST));
  AthenaArray<Real> dbval;
  dbval.NewAthenaArray(nvar,nx2,(NGHOST));
  int nyzone = (je-js)+1;

  // set primitive variables in inlet ghost zones
  for(int j=0; j<nx2; ++j) {
	jo = fmod(j+(jplus-(NGHOST)), nyzone) + (NGHOST);
	if (jo < js) jo += nyzone;
	if (jo > je) jo -= nyzone;
	if ((jo < js) || (jo > je)){
	  std::cout << "[ShearOuterX1]: jo = " << jo << " < js = " << js << "or > je= " << je << std::endl;
	  exit(1);
	}
    for(int i=1; i<=(NGHOST); ++i) {
      bval(IDN,j,i-1) = a(IDN,ks,jo,is+i-1);
      bval(IVX,j,i-1) = a(IVX,ks,jo,is+i-1);
      bval(IVY,j,i-1) = a(IVY,ks,jo,is+i-1);
      bval(IVZ,j,i-1) = a(IVZ,ks,jo,is+i-1);
      if (NON_BAROTROPIC_EOS) {
        bval(IEN,j,i-1) = a(IEN,ks,jo,is+i-1);
	  }
    }
  }
  // skip IEN and b-field;

  // zone-center slope dbval for the fractional part
  LinearSlope(nvar,nx2,(NGHOST),bval,dbval);
  // piece-wise linear interpolation of the fractional part

  for(int j=js; j<=je; ++j) {
    for(int i=1; i<=(NGHOST); ++i) {
	  int ib =  i - 1;
	  a(IDN,ks,j,ie+i) = (1.0-frac)*bval(IDN,j,ib) + frac*bval(IDN,j+1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IDN,j,ib)-dbval(IDN,j+1,ib));
	  a(IVX,ks,j,ie+i) = (1.0-frac)*bval(IVX,j,ib) + frac*bval(IVX,j+1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IVX,j,ib)-dbval(IVX,j+1,ib));
	  a(IVY,ks,j,ie+i) = (1.0-frac)*bval(IVY,j,ib) + frac*bval(IVY,j+1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IVY,j,ib)-dbval(IVY,j+1,ib)) -
					 //qshear*Omega_0*x1size*a(IDN,ks,j,ie+i);
					 qshear*Omega_0*x1size;
	  a(IVZ,ks,j,ie+i) = (1.0-frac)*bval(IVZ,j,ib) + frac*bval(IVZ,j+1,ib) +
		             0.5*frac*(1.0-frac)*(dbval(IVZ,j,ib)-dbval(IVZ,j+1,ib));
      if (NON_BAROTROPIC_EOS) {
	    a(IEN,ks,j,ie+i) = (1.0-frac)*bval(IEN,j,ib) + frac*bval(IEN,j+1,ib) +
		               0.5*frac*(1.0-frac)*(dbval(IEN,j,ib)-dbval(IEN,j+1,ib));
	  }
	}}

}

//--------------------------------------------------------------------------------------
//! \fn void LinearSlope(const int nvar, const int ny, const int nx, const AthenaArray<Real> &w, const AthenaArray<Real> &dw)
//  \brief cell center monotonic slope.
void LinearSlope(const int nvar, const int ny, const int nx, const AthenaArray<Real> &w, AthenaArray<Real> &dw) {
//
//
// This is an implementation of van Leer's harmonic mean
// monotonic slope.  Other slope monotonic slope functions
// could be substituted, if desired.
//
// Input W:  state variable
//
// Output :
//    dwc =  centered non-monotonic slope
//    dwl =  left-centered non-monotonic slope
//    dwr =  right-centered non-monotonic slope
//    dw  =  centered monotonic slope
//

  Real dwl, dwr;

  for(int k=0; k<nvar; ++k) {
    for(int j=1; j<ny-1; ++j)  {
      for(int i=0; i<nx; ++i) {
        dwl = w(k,j,i) - w(k,j-1,i);
        dwr = w(k,j+1,i) - w(k,j,i);
        if (dwl*dwr > 0.0) {
          dw(k,j,i) = 2.0*dwl*dwr/(dwl+dwr);
        } else {
          dw(k,j,i) = 0.0;
        }}
    }
    for(int i=0; i<nx; ++i) {
      dw(k,0,i) = 0.0;
	  dw(k,ny-1,i) = 0.0;
    }
  }

  return;
}
//JMSHI]
