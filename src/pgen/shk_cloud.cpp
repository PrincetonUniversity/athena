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
//! \file shk_cloud.c
//  \brief Problem generator for shock-cloud problem
//
// The shock-cloud problem consists of a planar shock impacting a single spherical cloud
// Input parameters are:
//    - problem/Mach   = Mach number of incident shock
//    - problem/drat   = density ratio of cloud to ambient
//    - problem/beta   = ratio of Pgas/Pmag
//
// The cloud radius is fixed at 1.0.  The center of the coordinate system defines the
// center of the cloud, and should be in the middle of the cloud. The shock is initially
// at x1=-2.0.  A typical grid domain should span x1 in [-3.0,7.0] , y and z in 
//[-2.5,2.5] (see input file in /tst).
//======================================================================================

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

// postshock flow variables are shared with IIB function
static Real gmma1,dl,pl,ul;
static Real bxl,byl,bzl;

// fixes BCs on L-x1 (left edge) of grid to postshock flow.
void ShockCloudInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                       FaceField &b, int is, int ie, int js, int je, int ks, int ke);

//======================================================================================
//! \fn void Mesh::InitUserMeshProperties(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
// Set IIB value function pointer
  EnrollUserBoundaryFunction(INNER_X1, ShockCloudInnerX1);
  return;
}


//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================

void Mesh::TerminateUserMeshProperties(void)
{
  // nothing to do
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the shock-cloud interaction test
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real gmma  = phydro->peos->GetGamma();
  gmma1 = gmma - 1.0;

// Read input parameters

  Real xshock = -2.0;
  Real rad    = 1.0;
  Real mach = pin->GetReal("problem","Mach");
  Real drat = pin->GetReal("problem","drat");
  Real beta;
  if (MAGNETIC_FIELDS_ENABLED) beta = pin->GetReal("problem","beta");
  
// Set paramters in ambient medium ("R-state" for shock)

  Real dr = 1.0;
  Real pr = 1.0/(phydro->peos->GetGamma());
  Real ur = 0.0;

// Uses Rankine Hugoniot relations for adiabatic gas to initialize problem

  Real jump1 = (gmma + 1.0)/(gmma1 + 2.0/(mach*mach));
  Real jump2 = (2.0*gmma*mach*mach - gmma1)/(gmma + 1.0);
  Real jump3 = 2.0*(1.0 - 1.0/(mach*mach))/(gmma + 1.0);

  dl = dr*jump1;
  pl = pr*jump2;
  ul = ur + jump3*mach*sqrt(gmma*pr/dr);

// Initialize the grid

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    // postshock flow
    if(pcoord->x1v(i) < xshock) {
      phydro->u(IDN,k,j,i) = dl;
      phydro->u(IM1,k,j,i) = ul*dl;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = pl/gmma1 + 0.5*dl*(ul*ul);

    // preshock ambient gas
    } else {
      phydro->u(IDN,k,j,i) = dr;
      phydro->u(IM1,k,j,i) = ur*dr;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = pr/gmma1 + 0.5*dr*(ur*ur);
    }

    // cloud interior
    Real diag = sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)) + SQR(pcoord->x3v(k)));
    if (diag < rad) {
      phydro->u(IDN,k,j,i) = dr*drat;
      phydro->u(IM1,k,j,i) = ur*dr*drat;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = pr/gmma1 + 0.5*dr*drat*(ur*ur);
    }
  }}}

// initialize interface B, assuming longitudinal field only B=(1,0,0)

  if (MAGNETIC_FIELDS_ENABLED) {
    Real bxr = sqrt(2.0/beta);
    Real byr = 0.0;
    Real bzr = 0.0;
    bxl = sqrt(2.0/beta);
    byl = 0.0;
    bzl = 0.0;

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      if(pcoord->x1v(i) < xshock) {
        pfield->b.x1f(k,j,i) = bxl;
      } else {
        pfield->b.x1f(k,j,i) = bxr;
      }
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      if(pcoord->x1v(i) < xshock) {
        pfield->b.x2f(k,j,i) = byl;
      } else {
        pfield->b.x2f(k,j,i) = byr;
      }
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      if(pcoord->x1v(i) < xshock) {
        pfield->b.x3f(k,j,i) = bzl;
      } else {
        pfield->b.x3f(k,j,i) = bzr;
      }
    }}}

    // initialize total energy

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      if(pcoord->x1v(i) < xshock) {
        phydro->u(IEN,k,j,i) += 0.5*(bxl*bxl + byl*byl + bzl*bzl);
      } else {
        phydro->u(IEN,k,j,i) += 0.5*(bxr*bxr + byr*byr + bxr*bzr);
      }
    }}}
  }
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
//! \fn void ShockCloudInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) 
// Note quantities at this boundary are held fixed at the downstream state

void ShockCloudInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                       FaceField &b, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1; i<=(NGHOST); ++i) {
      a(IDN,k,j,is-i) = dl;
      a(IVX,k,j,is-i) = ul;
      a(IVY,k,j,is-i) = 0.0;
      a(IVZ,k,j,is-i) = 0.0;
      a(IEN,k,j,is-i) = pl;
    }
  }}
}
