//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cart_disk.cpp
//  \brief Initializes Keplerian accretion disk in cartesian cooordinates

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

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
static Real PoverR(const Real rad, const Real phi, const Real z);
static void VelProfileCyl(const Real rad, const Real phi, const Real z,
  Real &v1, Real &v2, Real &v3);
static Real RampFunc(const Real rad, const Real phi, const Real z,
					 const Real v1, const Real v2, const Real v3);
//User source terms
void Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
//
// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

// problem parameters which are useful to make global to this file
static Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
static Real dfloor;
static Real rsoft,rin,rout;
static Real tsink; // mass removing time scale
static Real rsink;
static Real beta_th;//dimensionless cooling time
static Real rbuf1,rbuf2;
static Real alpha;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);
  rsoft = pin->GetOrAddReal("problem","rsoft",0.0);
  rsink = pin->GetOrAddReal("problem","rsink",0.0);
  rin = pin->GetOrAddReal("problem","rin",0.0);
  rout = pin->GetOrAddReal("problem","rout",5.0);
  tsink = pin->GetOrAddReal("problem","tsink",100.0);
  beta_th = pin->GetOrAddReal("problem","beta_th",0.0);
  alpha = pin->GetOrAddReal("problem","alpha",0.0);
  rbuf1 = pin->GetOrAddReal("problem","rbuf1",8.0);
  rbuf2 = pin->GetOrAddReal("problem","rbuf2",10.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure and cooling parameters
  if(NON_BAROTROPIC_EOS){
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  }else{
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));

  // Enroll user-defined physical source terms
  if (beta_th != 0.0 && NON_BAROTROPIC_EOS) EnrollUserExplicitSourceFunction(Cooling);


  // enroll user-defined boundary condition
  //if(mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
  //  EnrollUserBoundaryFunction(INNER_X1, DiskInnerX1);
  //}
  //if(mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
  //  EnrollUserBoundaryFunction(OUTER_X1, DiskOuterX1);
  //}
  //if(mesh_bcs[INNER_X2] == GetBoundaryFlag("user")) {
  //  EnrollUserBoundaryFunction(INNER_X2, DiskInnerX2);
  //}
  //if(mesh_bcs[OUTER_X2] == GetBoundaryFlag("user")) {
  //  EnrollUserBoundaryFunction(OUTER_X2, DiskOuterX2);
  //}
//  if(mesh_bcs[INNER_X3] == GetBoundaryFlag("user")) {
//    EnrollUserBoundaryFunction(INNER_X3, DiskInnerX3);
//  }
//  if(mesh_bcs[OUTER_X3] == GetBoundaryFlag("user")) {
//    EnrollUserBoundaryFunction(OUTER_X3, DiskOuterX3);
//  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real rad, phi, z;
  Real v1, v2, v3;

  //  Initialize density and momenta
  for(int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
      // compute initial conditions in cylindrical coordinates
	  // 1) pure power law with truncation radii rin rout
      //if (rad >= rin)
      //  phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
      //else
      //  phydro->u(IDN,k,j,i) = dfloor;
	  // 2) power law profile with exponential taper interior to rin
	  if (rad >= rsink)
	    phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z); //*exp(-pow((rad/rin),-3.0));
	  else
		phydro->u(IDN,k,j,i) = dfloor;
      VelProfileCyl(rad,phi,z,v1,v2,v3);
      phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
      phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
      phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
      if (NON_BAROTROPIC_EOS){
        Real p_over_r = PoverR(rad,phi,z);
        phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                   + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k)
{
  rad=sqrt(SQR(pco->x1v(i))+SQR(pco->x2v(j)));
  phi=atan2(pco->x2v(j),pco->x1v(i));
  z=pco->x3v(k);
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

static Real DenProfileCyl(const Real rad, const Real phi, const Real z)
{
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*pow(rad/r0,dslope);
  Real dentem = denmid*exp(gm0/p_over_r*(1./sqrt(SQR(rad)+SQR(rsoft)+SQR(z))-1./sqrt(SQR(rad)+SQR(rsoft))));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates

static Real PoverR(const Real rad, const Real phi, const Real z)
{
  Real poverr;
  poverr = p0_over_r0*pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

static void VelProfileCyl(const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3)
{
  Real p_over_r = PoverR(rad, phi, z);
  Real rad1 = sqrt(SQR(rad)+SQR(rsoft));
  Real rad2 = sqrt(SQR(rad)+SQR(rsoft)+SQR(z));
  Real vel = sqrt((dslope+pslope)*p_over_r + gm0*SQR(rad/rad1)/rad1
            +pslope*gm0*(1.0/rad1 -1.0/rad2));

  v1 = -vel*sin(phi);
  v2 = vel*cos(phi);
  v3 = 0.0;

  // correct with radial drift due to accretion
  if (NON_BAROTROPIC_EOS) {
	vel = 1.5*p_over_r*gamma_gas*alpha/sqrt(gm0/rad);
    v1 += -vel*cos(phi);
	v2 += -vel*sin(phi);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes ramp function R/\tau
// R(x) = a*x^2 + b*x  + c; where a=1/(r2^2-r1r2), and b=-a*r1
// tau  = 2PI/Omega(r1)

static Real RampFunc(const Real rad, const Real phi, const Real z,
					 const Real v1, const Real v2, const Real v3)
{
  Real ramp,tau;
  if (rad >= rbuf2) {
	ramp = 100.0;
  } else if (rad >= rbuf1) {
    Real fac = 1.0/(SQR(rbuf2)-rbuf1*rbuf2);
    ramp = fac*SQR(rad)-fac*rbuf1*rad;
  } else {
	  ramp = 0.0;
  }
  tau = 2.0*PI/sqrt(gm0/rbuf1/SQR(rbuf1));
  return ramp/tau;
}

void MeshBlock::UserWorkInLoop(void)
{
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real rad, phi, z;
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // apply sink cells within the rsink( ususally >= rsoft)
        if (rad <= rsink) { // sink cells within r<=rsink
          phydro->u(IDN,k,j,i) -= pmy_mesh->dt*phydro->u(IDN,k,j,i)/tsink;
        }
		// apply buffer zones within [rbuf1,rbuf2] to quench m=4 mode
		if (rad >= rbuf1) {
		  Real v1, v2, v3;
		  Real den0 = DenProfileCyl(rad,phi,z);
          VelProfileCyl(rad,phi,z,v1,v2,v3);
		  Real ramp = RampFunc(rad,phi,z,v1,v2,v3);
          phydro->u(IDN,k,j,i) -= pmy_mesh->dt*(phydro->u(IDN,k,j,i)-den0)*ramp;
          phydro->u(IM1,k,j,i) -= pmy_mesh->dt*(phydro->u(IM1,k,j,i)-den0*v1)*ramp;
          phydro->u(IM2,k,j,i) -= pmy_mesh->dt*(phydro->u(IM2,k,j,i)-den0*v2)*ramp;
          phydro->u(IM3,k,j,i) -= pmy_mesh->dt*(phydro->u(IM3,k,j,i)-den0*v3)*ramp;
          if (NON_BAROTROPIC_EOS){
            Real p_over_r = PoverR(rad,phi,z);
            Real eng0 = p_over_r*den0/(gamma_gas - 1.0);
			eng0 += 0.5*den0*(SQR(v1)+SQR(v2)+SQR(v3));
            phydro->u(IEN,k,j,i) -= pmy_mesh->dt*(phydro->u(IEN,k,j,i)-eng0)*ramp;
          }
		}

      }
    }
  }
}

void Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  // apply extremely short cooling
  //Real gam = pmb->peos->GetGamma();
  //Real gam1 = pmb->peos->GetGamma()-1.0;
  Real gam1 = gamma_gas-1.0;
  Real tau = 0.01;
  Real rad, phi, z;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
		Real eth = cons(IEN,k,j,i)-0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i))+
										SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);
		Real eth0 = cons(IDN,k,j,i)*PoverR(rad,phi,z)/gam1;
		Real tcool= std::max(beta_th*2.0*PI/sqrt(gm0/SQR(rad)/rad),dt);
        cons(IEN,k,j,i) -= dt*(eth-eth0)*dt/tcool;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
//

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        GetCylCoord(pco,rad,phi,z,is-i,j,k);
        prim(IDN,k,j,is-i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,j,is-i) = v1;
        prim(IM2,k,j,is-i) = v2;
        prim(IM3,k,j,is-i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,is-i) = PoverR(rad, phi, z)*prim(IDN,k,j,is-i);
      }
    }
  }
}

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        GetCylCoord(pco,rad,phi,z,ie+i,j,k);
        prim(IDN,k,j,ie+i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,j,ie+i) = v1;
        prim(IM2,k,j,ie+i) = v2;
        prim(IM3,k,j,ie+i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,ie+i) = PoverR(rad, phi, z)*prim(IDN,k,j,ie+i);
      }
    }
  }
}

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,js-j,k);
        prim(IDN,k,js-j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,js-j,i) = v1;
        prim(IM2,k,js-j,i) = v2;
        prim(IM3,k,js-j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,js-j,i) = PoverR(rad, phi, z)*prim(IDN,k,js-j,i);
      }
    }
  }
}

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,je+j,k);
        prim(IDN,k,je+j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,k,je+j,i) = v1;
        prim(IM2,k,je+j,i) = v2;
        prim(IM3,k,je+j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,je+j,i) = PoverR(rad, phi, z)*prim(IDN,k,je+j,i);
      }
    }
  }
}

//void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
//{
//  Real rad,phi,z;
//  Real v1, v2, v3;
//  for (int k=1; k<=(NGHOST); ++k) {
//    for (int j=js; j<=je; ++j) {
//      for (int i=is; i<=ie; ++i) {
//        GetCylCoord(pco,rad,phi,z,i,j,ks-k);
//        prim(IDN,ks-k,j,i) = DenProfileCyl(rad,phi,z);
//        VelProfileCyl(rad,phi,z,v1,v2,v3);
//        prim(IM1,ks-k,j,i) = v1;
//        prim(IM2,ks-k,j,i) = v2;
//        prim(IM3,ks-k,j,i) = v3;
//        if (NON_BAROTROPIC_EOS)
//          prim(IEN,ks-k,j,i) = PoverR(rad, phi, z)*prim(IDN,ks-k,j,i);
//      }
//    }
//  }
//}
//
//void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
//{
//  Real rad,phi,z;
//  Real v1, v2, v3;
//  for (int k=1; k<=(NGHOST); ++k) {
//    for (int j=js; j<=je; ++j) {
//      for (int i=is; i<=ie; ++i) {
//        GetCylCoord(pco,rad,phi,z,i,j,ke+k);
//        prim(IDN,ke+k,j,i) = DenProfileCyl(rad,phi,z);
//        VelProfileCyl(rad,phi,z,v1,v2,v3);
//        prim(IM1,ke+k,j,i) = v1;
//        prim(IM2,ke+k,j,i) = v2;
//        prim(IM3,ke+k,j,i) = v3;
//        if (NON_BAROTROPIC_EOS)
//          prim(IEN,ke+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ke+k,j,i);
//      }
//    }
//  }
//}
