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
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "../fluid/fluid.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/eos/eos.hpp"    // ParameterInput

// inline functions
inline Real SQR(Real x){
  return ((x)*(x));
}

//======================================================================================
/*! \file disk_cyl.cpp
 *  \brief Problem generator for accretion disk problems.  
 *
 * Problem generator for disk problems in disk_cyl system.
 *====================================================================================*/

// File scope variables
static Real x1Max, x1Min;
static Real rho0, rho_floor, rho_MIN;
static int ICd, ICv;
static Real rstart=0.0,rtrunc=0.0;
static Real q_m, R0 = 1.0;

// Function Declarations
static Real ICden(const Real x1);
static Real ICvel(const Real x1);
static Real KeplerVel(const Real x1);

void Fluid::InitFluid(ParameterInput *pin)
{
  MeshBlock *pb = pmy_block;
  std::stringstream msg;

  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;
  x1Max = pb->block_size.x1max;
  x1Min = pb->block_size.x1min;

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

// Get initial pressure
  Real pressure = 0.0;
  if(NON_BAROTROPIC_EOS){
    pressure = pin->GetReal("problem","pressure");
  }
std::cout<<"[Problem Generator]: pressure="<<pressure<<std::endl;

// Get options for boundary conditions
  int Hbc = pin->GetInteger("problem","Hbc");
  int ix1_bc = pin->GetInteger("mesh","ix1_bc");
  int ox1_bc = pin->GetInteger("mesh","ox1_bc");

// Initialize the disk
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pb->x1v(i);
        u(IDN,k,j,i) = ICden(x1);
        u(IM1,k,j,i) = 0.0;
        u(IM2,k,j,i) = ICvel(x1)*u(IDN,k,j,i);
        u(IM3,k,j,i) = 0.0;
	if(NON_BAROTROPIC_EOS) u(IEN,k,j,i) = pressure/(pf_eos->GetGamma() - 1.0) +
				0.5*(SQR(u(IM1,k,j,i))+SQR(u(IM2,k,j,i))+SQR(u(IM3,k,j,i)))/u(IDN,k,j,i);
      }
    }
  }


  return;
}

/* Initial Density Profile*/
Real ICden(const Real x1){
  Real den = 0.0;
  std::stringstream msg;
  
  /*Initial Condition 1: for no inflow*/
  if(ICd == 1) {
    den = rho0;
  } 
  /*Initial Condition 2: for L1-inflow*/
  if(ICd == 2) {
    if(x1<=0.2*x1Max) den = 0.2*rho0;
    else den = rho_floor;
  } 
  /*Initial Condition 3: empty initial disk*/
  if(ICd == 3) {
    den = rho_floor;
  } 
  /*Initial Condition 4: truncated initial disk*/
  if(ICd == 4) {
    if(x1>=rstart && x1<=rtrunc) den = rho0;
    else den = rho_floor;
  } 
  /*Initial Condition 5: truncated initial disk*/
  if(ICd == 5) {
    if(x1>=rstart && x1<=rtrunc) den = rho0;
    else if(x1 > rtrunc)
      den = rho0*exp(-SQR(x1-rtrunc)/2./SQR(0.1*(x1Max-rtrunc)));
    else if(x1 < rstart)
      den = rho0*exp(-SQR(x1-rstart)/2./SQR(0.1*(x1Min-rstart)));
  } /*ICd 5*/
  
  if((ICd!=1)&&(ICd!=2)&&(ICd!=3)&&(ICd!=4)&&(ICd!=5)){
    msg << "### FATAL ERROR in Problem Generator" << std::endl 
        << "The ICd flag value "<<ICd <<" is not defined"<<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return den;
}

/* Initial Velocity Profile*/
Real ICvel(const Real x1){
  Real vel = 0.0;
  std::stringstream msg;

  /*Initial Condition 1: static gas*/
  if(ICv == 1) {
    vel = 0.0;
  } 
  /*Initial Condition 2: moving gas*/
  else if(ICv == 2) {
    Real rtrans=x1Min+0.1*(x1Max-x1Min);
    if(x1 < rtrans)
      vel = 0;
    else 
      vel = sqrt(x1-rtrans);
  }
  /*Initial Condition 2: moving gas*/
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
  Real v = 1.0/sqrt(x1); 
  return v;
}


