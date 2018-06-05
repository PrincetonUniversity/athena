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

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
static Real PoverR(const Real rad, const Real phi, const Real z);
static void VelProfileCyl(const Real rad, const Real phi, const Real z,
  Real &v1, Real &v2, Real &v3);

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

// problem parameters which are useful to make global to this file
static Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
static Real dfloor;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));

  // enroll user-defined boundary condition
  if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, DiskInnerX1);
  }
  if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X1, DiskOuterX1);
  }
  if (mesh_bcs[INNER_X2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X2, DiskInnerX2);
  }
  if (mesh_bcs[OUTER_X2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X2, DiskOuterX2);
  }
  if (mesh_bcs[INNER_X3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X3, DiskInnerX3);
  }
  if (mesh_bcs[OUTER_X3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X3, DiskOuterX3);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad, phi, z;
  Real v1, v2, v3;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
      // compute initial conditions in cylindrical coordinates
      phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
      VelProfileCyl(rad,phi,z,v1,v2,v3);

      phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
      phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
      phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
      if (NON_BAROTROPIC_EOS) {
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

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (COORDINATE_SYSTEM == "cylindrical") {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    rad=fabs(pco->x1v(i)*sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

static Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*pow(rad/r0,dslope);
  Real dentem = denmid*exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates

static Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

static void VelProfileCyl(const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(gm0/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel);
  if (COORDINATE_SYSTEM == "cylindrical") {
    v1=0.0;
    v2=vel;
    v3=0.0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    v1=0.0;
    v2=0.0;
    v3=vel;
  }
  return;
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
//

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
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
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
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
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
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
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
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

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,ks-k);
        prim(IDN,ks-k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,ks-k,j,i) = v1;
        prim(IM2,ks-k,j,i) = v2;
        prim(IM3,ks-k,j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,ks-k,j,i) = PoverR(rad, phi, z)*prim(IDN,ks-k,j,i);
      }
    }
  }
}

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,ke+k);
        prim(IDN,ke+k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(rad,phi,z,v1,v2,v3);
        prim(IM1,ke+k,j,i) = v1;
        prim(IM2,ke+k,j,i) = v2;
        prim(IM3,ke+k,j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,ke+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ke+k,j,i);
      }
    }
  }
}
