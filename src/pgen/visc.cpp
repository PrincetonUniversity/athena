//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file disk.cpp
//  \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//  spherical polar coordinates.  Initial conditions are in vertical hydrostatic eqm.

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min
#include <cstdlib>    // srand
#include <cfloat>     // FLT_MIN

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

static void GetCartCoord(Coordinates *pco,Real &x1,Real &x2,Real &x3,int i,int j,int k);
//static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
//static Real PoverR(const Real rad, const Real phi, const Real z);
//static void VelProfileCyl(const Real rad, const Real phi, const Real z,
//  Real &v1, Real &v2, Real &v3);
//
// User-defined boundary conditions for disk simulations
//void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

// problem parameters which are useful to make global to this file
static Real v0,t0;
static Real nuiso,gm0;
static int iprob;
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Get parameters for gravitatonal potential of central point mass
  v0 = pin->GetOrAddReal("problem","v0",0.001);
  t0 = pin->GetOrAddReal("problem","t0",0.5);
  nuiso = pin->GetOrAddReal("problem","nuiso",0.0);
  iprob = pin->GetOrAddInteger("problem","iprob",0);
  gm0 = pin->GetOrAddReal("problem","GM", 0.0);

//  // enroll user-defined boundary condition
//  if(mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
//    EnrollUserBoundaryFunction(INNER_X1, DiskInnerX1);
//  }
//  if(mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
//    EnrollUserBoundaryFunction(OUTER_X1, DiskOuterX1);
//  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes viscous shear flow.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  //Real rad, phi, z;
  Real v1=0.0, v2=0.0, v3=0.0;
  Real d0 = 1.0, p0=1.0, x0=1.0;
  Real x1,x2,x3;
  Real rad,z,phi,theta;

  //  Initialize density and momenta in Cartesian grids
  if (iprob == 0) { //visc column
	if (nuiso == 0) nuiso = 0.03;
    for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //GetCartCoord(pcoord,x1,x2,x3,i,j,k); // convert to cylindrical coordinates
        if(COORDINATE_SYSTEM == "cartesian"){
          x1=pcoord->x1v(i);
          x2=pcoord->x2v(j);
          x3=pcoord->x3v(k);
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
          v2 = v0/sqrt(4.0*PI*nuiso*t0)*exp(-SQR(x1-x0)/(4.0*nuiso*t0));
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
        }
        if(COORDINATE_SYSTEM == "cylindrical"){
          rad=pcoord->x1v(i);
          phi=pcoord->x2v(j);
          z=pcoord->x3v(k);
          x1=rad*cos(phi);
          x2=rad*sin(phi);
          x3=z;
          v2 = v0/sqrt(4.0*PI*nuiso*t0)*exp(-SQR(x1-x0)/(4.0*nuiso*t0));
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(v1*cos(phi)+v2*sin(phi));
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*(-v1*sin(phi)+v2*cos(phi));
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
        }
        if(COORDINATE_SYSTEM == "spherical_polar"){
          rad=fabs(pcoord->x1v(i)*sin(pcoord->x2v(j)));
          phi=pcoord->x3v(k);
          theta=pcoord->x2v(j);
          z=pcoord->x1v(i)*cos(pcoord->x2v(j));
          x1=rad*cos(phi);
          x2=rad*sin(phi);
          x3=z;
          v2 = v0/sqrt(4.0*PI*nuiso*t0)*exp(-SQR(x1-x0)/(4.0*nuiso*t0));
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(v1*cos(phi)*sin(theta)+v2*sin(phi)*sin(theta)+v3*cos(theta));
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*(v1*cos(phi)*cos(theta)+v2*sin(phi)*cos(theta)+v3*sin(theta));
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*(-v1*sin(phi)+v2*cos(phi));
        }
      }
    }}
  } else if (iprob == 1) { // visc ring gaussian profile at t=0
     if(COORDINATE_SYSTEM != "cylindrical" && gm0 == 0.0) {
       std::cout << "[visc]: viscous ring test only compatible with cylindrical coord with point mass in center" << std::endl;
       exit(0);
     }
	 Real width = 0.1;
     for(int k=ks; k<=ke; ++k) {
     for (int j=js; j<=je; ++j) {
       for (int i=is; i<=ie; ++i) {
         rad=pcoord->x1v(i);
         d0 = exp(-SQR(rad-x0)/2.0/SQR(width));
         if(d0 < 1e-6) d0 += 1e-6;
		 v1 = -3.0*nuiso*(0.5/rad - (rad-x0)/SQR(width));
         v2 = sqrt(gm0/rad); // set it to be Mach=100
         phydro->u(IDN,k,j,i) = d0;
         phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
         phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
         phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
    }}}
  } else {
    std::cout << "[visc]: iprob has to be either 0 or 1" << std::endl;
    exit(0);
  }

  return;
}

//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

static void GetCartCoord(Coordinates *pco,Real &x1,Real &x2,Real &x3,int i,int j,int k)
{
  Real rad,phi,z;
  if(COORDINATE_SYSTEM == "cylindrical"){
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
    x1=rad*cos(phi);
    x2=rad*sin(phi);
    x3=z;
  } else if(COORDINATE_SYSTEM == "spherical_polar"){
    rad=fabs(pco->x1v(i)*sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*cos(pco->x2v(j));
    x1=rad*cos(phi);
    x2=rad*sin(phi);
    x3=z;
  } else if(COORDINATE_SYSTEM == "cartesian"){
    x1=pco->x1v(i);
    x2=pco->x2v(j);
    x3=pco->x3v(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
////
//void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
//{
//  Real rad,phi,z;
//  Real v1=, v2, v3;
//  for (int k=ks; k<=ke; ++k) {
//    for (int j=js; j<=je; ++j) {
//      for (int i=1; i<=(NGHOST); ++i) {
//        GetCylCoord(pco,rad,phi,z,is-i,j,k);
//        prim(IDN,k,j,is-i) = DenProfileCyl(rad,phi,z);
//        VelProfileCyl(rad,phi,z,v1,v2,v3);
//        prim(IM1,k,j,is-i) = v1;
//        prim(IM2,k,j,is-i) = v2;
//        prim(IM3,k,j,is-i) = v3;
//        if (NON_BAROTROPIC_EOS)
//          prim(IEN,k,j,is-i) = PoverR(rad, phi, z)*prim(IDN,k,j,is-i);
//      }
//    }
//  }
//}
//
//void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
//{
//  Real rad,phi,z;
//  Real v1, v2, v3;
//  for (int k=ks; k<=ke; ++k) {
//    for (int j=js; j<=je; ++j) {
//      for (int i=1; i<=(NGHOST); ++i) {
//        GetCylCoord(pco,rad,phi,z,ie+i,j,k);
//        prim(IDN,k,j,ie+i) = DenProfileCyl(rad,phi,z);
//        VelProfileCyl(rad,phi,z,v1,v2,v3);
//        prim(IM1,k,j,ie+i) = v1;
//        prim(IM2,k,j,ie+i) = v2;
//        prim(IM3,k,j,ie+i) = v3;
//        if (NON_BAROTROPIC_EOS)
//          prim(IEN,k,j,ie+i) = PoverR(rad, phi, z)*prim(IDN,k,j,ie+i);
//      }
//    }
//  }
//}

