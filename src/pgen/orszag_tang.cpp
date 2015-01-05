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

// Primary header
#include "../mesh.hpp" 

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/eos/eos.hpp"    // eos
#include "../fluid/fluid.hpp"      // Fluid
#include "../field/field.hpp"      // Field

//======================================================================================
//! \file orszag-tang.c
//  \brief Problem generator for Orszag-Tang vortex problem.
//
// REFERENCE: For example, see: G. Toth,  "The div(B)=0 constraint in shock capturing
//   MHD codes", JCP, 161, 605 (2000)				      */
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real gm1 = (pfl->pf_eos->GetGamma() - 1.0);

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
    az(j,i) = B0/(4.0*PI)*cos(4.0*PI*pmb->x1f(i)) + B0/(2.0*PI)*cos(2.0*PI*pmb->x2f(j));
  }}

// Initialize density, momentum, face-centered fields

  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfl->u(IDN,ks,j,i) = d0;
    pfl->u(IM1,ks,j,i) = -d0*v0*sin(2.0*PI*pmb->x2v(j));
    pfl->u(IM2,ks,j,i) =  d0*v0*sin(2.0*PI*pmb->x1v(i));
    pfl->u(IM3,ks,j,i) = 0.0;
  }}

// initialize interface B

  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfd->b.x1f(ks,j,i) = (az(j+1,i) - az(j,i))/pmb->dx2f(j);
  }}
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfd->b.x2f(ks,j,i) = (az(j,i) - az(j,i+1))/pmb->dx1f(i);
  }}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfd->b.x3f(k,j,i) = 0.0;
  }}}

// initialize total energy

  if (NON_BAROTROPIC_EOS) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfl->u(IEN,ks,j,i) = p0/gm1 +
          0.5*(SQR(0.5*(pfd->b.x1f(ks,j,i) + pfd->b.x1f(ks,j,i+1))) +
               SQR(0.5*(pfd->b.x2f(ks,j,i) + pfd->b.x2f(ks,j+1,i))) +
               SQR(0.5*(pfd->b.x3f(ks,j,i) + pfd->b.x3f(ks+1,j,i)))) + (0.5)*
          (SQR(pfl->u(IM1,ks,j,i)) + SQR(pfl->u(IM2,ks,j,i)) + SQR(pfl->u(IM3,ks,j,i)))
          /pfl->u(IDN,ks,j,i);
    }}
  }

  az.DeleteAthenaArray();
  return;
}
