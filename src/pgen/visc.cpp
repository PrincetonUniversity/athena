//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file visc.cpp
//  iprob = 0 - test viscous shear flow density column in various coordinate systems
//  iprob = 1 - test viscous spreading of Keplerain ring

// C++ headers
#include <algorithm>  // min
#include <cstdlib>    // srand
#include <cfloat>     // FLT_MIN
#include <cmath>      // std::sqrt()
#include <fstream>
#include <iostream>   // endl
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <sstream>    // stringstream

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

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  v0 = pin->GetOrAddReal("problem","v0",0.001);
  t0 = pin->GetOrAddReal("problem","t0",0.5);
  nuiso = pin->GetOrAddReal("problem","nuiso",0.0);
  iprob = pin->GetOrAddInteger("problem","iprob",0);
  gm0 = pin->GetOrAddReal("problem","GM", 0.0);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes viscous shear flow.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //Real rad, phi, z;
  Real v1=0.0, v2=0.0, v3=0.0;
  Real d0 = 1.0, p0=1.0, x0=1.0;
  Real x1,x2,x3;
  Real rad,z,phi,theta;

  //  Initialize density and momenta in Cartesian grids
  if (iprob == 0) { //visc column
    if (nuiso == 0) {
      nuiso = 0.03;
    }
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {
            x1=pcoord->x1v(i);
            x2=pcoord->x2v(j);
            x3=pcoord->x3v(k);
            phydro->u(IDN,k,j,i) = d0;
            phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
            v2 = v0/std::sqrt(4.0*PI*nuiso*t0)*exp(-SQR(x1-x0)/(4.0*nuiso*t0));
            phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
            phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            rad=pcoord->x1v(i);
            phi=pcoord->x2v(j);
            z=pcoord->x3v(k);
            x1=rad*cos(phi);
            x2=rad*sin(phi);
            x3=z;
            v2 = v0/std::sqrt(4.0*PI*nuiso*t0)*exp(-SQR(x1-x0)/(4.0*nuiso*t0));
            phydro->u(IDN,k,j,i) = d0;
            phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(v1*cos(phi)+v2*sin(phi));
            phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*(-v1*sin(phi)+v2*cos(phi));
            phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
          } else { // (COORDINATE_SYSTEM == "spherical_polar") {
            rad=fabs(pcoord->x1v(i)*sin(pcoord->x2v(j)));
            phi=pcoord->x3v(k);
            theta=pcoord->x2v(j);
            z=pcoord->x1v(i)*cos(pcoord->x2v(j));
            x1=rad*cos(phi);
            x2=rad*sin(phi);
            x3=z;
            v2 = v0/std::sqrt(4.0*PI*nuiso*t0)*exp(-SQR(x1-x0)/(4.0*nuiso*t0));
            phydro->u(IDN,k,j,i) = d0;
            phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(v1*cos(phi)*sin(theta)
                                                         +v2*sin(phi)*sin(theta)
                                                         +v3*cos(theta));
            phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*(v1*cos(phi)*cos(theta)
                                                         +v2*sin(phi)*cos(theta)
                                                         +v3*sin(theta));
            phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*(-v1*sin(phi)+v2*cos(phi));
          }
        }
      }
    }
  } else if (iprob == 1) { // visc ring gaussian profile at t=0
    if (COORDINATE_SYSTEM != "cylindrical" && gm0 == 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in visc.cpp ProblemGenerator" << std::endl
          << "viscous ring test only compatible with cylindrical coord"
          << std::endl << "with point mass in center" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    Real width = 0.1;
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          rad=pcoord->x1v(i);
          d0 = exp(-SQR(rad-x0)/2.0/SQR(width));
          if (d0 < 1e-6) d0 += 1e-6;
          v1 = -3.0*nuiso*(0.5/rad - (rad-x0)/SQR(width));
          v2 = std::sqrt(gm0/rad); // set it to be Mach=100
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
        }
      }
    }
  } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in visc.cpp ProblemGenerator" << std::endl
          << "viscous iprob has to be either 0 or 1" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  return;
}
