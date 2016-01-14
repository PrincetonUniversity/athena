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
//! \file disk_cyl.cpp
//  \brief Problem generator for accretion disk problems.  
//
// Problem generator for disk problems in cylindrical coordinate system.
//======================================================================================

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

#error "disk_cyl.cpp is outdated and must be rewritten."

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// File scope variables
static Real x1Max, x1Min;
static Real rho0, rho_floor, rho_MIN;
static int ICd, ICv;
static Real rstart=0.0,rtrunc=0.0;
static Real GM=0.0, q_m=0.0, R0 = 1.0;

// Function Declarations
static Real ICden(const Real x1);
static Real ICvel(const Real x1);
static Real KeplerVel(const Real x1);


void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  x1Max = mesh_size.x1max;
  x1Min = mesh_size.x1min;

// Get parameters for grav terms
  GM = pin->GetOrAddReal("problem","GM",0.0);

// Get initial density
  rho0 = pin->GetReal("problem","rho0");
  rho_floor = pin->GetReal("problem","rho_floor"); 
  rho_MIN = 0.1*rho_floor;
  ICd = pin->GetInteger("problem","ICd");
  if(ICd == 4) {
    rstart = pin->GetReal("problem","rstart");
    if (rstart == 0) rstart = x1Min;
    rtrunc = pin->GetReal("problem","rtrunc");
  }
  ICv = pin->GetInteger("problem","ICv");

  return;
}


void Mesh::TerminateUserMeshProperties(void)
{
  return;
}


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;

// Get initial pressure
  Real pressure = 0.0;
  if(NON_BAROTROPIC_EOS){
    pressure = pin->GetReal("problem","pressure");
  }

// Initialize the disk
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = ICden(x1);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = ICvel(x1)*phydro->u(IDN,k,j,i);
        phydro->u(IM3,k,j,i) = 0.0;
	    if(NON_BAROTROPIC_EOS)
          phydro->u(IEN,k,j,i)= pressure/(phydro->peos->GetGamma()-1.0)
          + 0.5*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) +
                 SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }
  }


  return;
}


void MeshBlock::UserWorkInLoop(void)
{
  return;
}

// Initial Density Profile
Real ICden(const Real x1){
  Real den = 0.0;
  std::stringstream msg;
  
  //Initial Condition 1: for no inflow
  if(ICd == 1) {
    den = rho0;
  } 
  //Initial Condition 2: for L1-inflow
  if(ICd == 2) {
    if(x1<=0.2*x1Max) den = 0.2*rho0;
    else den = rho_floor;
  } 
  //Initial Condition 3: empty initial disk
  if(ICd == 3) {
    den = rho_floor;
  } 
  //Initial Condition 4: truncated initial disk
  if(ICd == 4) {
    if(x1>=rstart && x1<=rtrunc) den = rho0;
    else den = rho_floor;
  } 
  //Initial Condition 5: truncated initial disk
  if(ICd == 5) {
    if(x1>=rstart && x1<=rtrunc) den = rho0;
    else if(x1 > rtrunc)
      den = rho0*exp(-SQR(x1-rtrunc)/2./SQR(0.1*(x1Max-rtrunc)));
    else if(x1 < rstart)
      den = rho0*exp(-SQR(x1-rstart)/2./SQR(0.1*(x1Min-rstart)));
  } //ICd 5
  
  if((ICd!=1)&&(ICd!=2)&&(ICd!=3)&&(ICd!=4)&&(ICd!=5)){
    msg << "### FATAL ERROR in Problem Generator" << std::endl 
        << "The ICd flag value "<<ICd <<" is not defined"<<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return den;
}

// Initial Velocity Profile
Real ICvel(const Real x1){
  Real vel = 0.0;
  std::stringstream msg;

  //Initial Condition 1: static gas
  if(ICv == 1) {
    vel = 0.0;
  } 
  //Initial Condition 2: moving gas
  else if(ICv == 2) {
    Real rtrans=x1Min+0.1*(x1Max-x1Min);
    if(x1 < rtrans)
      vel = 0;
    else 
      vel = sqrt(x1-rtrans);
  }
  //Initial Condition 2: moving gas
  else if(ICv == 3) {
    vel = KeplerVel(x1); 
  }
  else {
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "The ICv flag value "<<ICv <<" is not defined"<<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return vel;
}


// Kepler Velocity
static Real KeplerVel(const Real x1){
  Real v = GM/sqrt(x1); 
  return v;
}


