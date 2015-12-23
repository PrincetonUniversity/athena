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
//! \file orszag-tang.c
//  \brief Problem generator for Orszag-Tang vortex problem.
//
// REFERENCE: For example, see: G. Toth,  "The div(B)=0 constraint in shock capturing
//   MHD codes", JCP, 161, 605 (2000)				      */
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

void Mesh::ProblemGenerator(Hydro *phyd, Field *pfld, ParameterInput *pin)
{
  MeshBlock *pmb = phyd->pmy_block;
  Coordinates *pco = pmb->pcoord;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real gm1 = (phyd->peos->GetGamma() - 1.0);

  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx2,nx1);

  Real B0 = 1.0/sqrt(4.0*PI);
  Real d0 = 25.0/(36.0*PI);
  Real v0 = 1.0;
  Real p0 = 5.0/(12*PI);

// Initialize vector potential

  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    az(j,i) = B0/(4.0*PI)*cos(4.0*PI*pco->x1f(i)) + B0/(2.0*PI)*cos(2.0*PI*pco->x2f(j));
  }}

// Initialize density, momentum, face-centered fields

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    phyd->u(IDN,k,j,i) = d0;
    phyd->u(IM1,k,j,i) = -d0*v0*sin(2.0*PI*pco->x2v(j));
    phyd->u(IM2,k,j,i) =  d0*v0*sin(2.0*PI*pco->x1v(i));
    phyd->u(IM3,k,j,i) = 0.0;
  }}}

// initialize interface B

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfld->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pco->dx2f(j);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfld->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pco->dx1f(i);
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfld->b.x3f(k,j,i) = 0.0;
  }}}

// initialize total energy

  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phyd->u(IEN,k,j,i) = p0/gm1 +
          0.5*(SQR(0.5*(pfld->b.x1f(k,j,i) + pfld->b.x1f(k,j,i+1))) +
               SQR(0.5*(pfld->b.x2f(k,j,i) + pfld->b.x2f(k,j+1,i))) +
               SQR(0.5*(pfld->b.x3f(k,j,i) + pfld->b.x3f(k+1,j,i)))) + (0.5)*
          (SQR(phyd->u(IM1,k,j,i)) + SQR(phyd->u(IM2,k,j,i)) + SQR(phyd->u(IM3,k,j,i)))
          /phyd->u(IDN,k,j,i);
    }}}
  }

  az.DeleteAthenaArray();
  return;
}
