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
//! \file rotor.c
//  \brief Sets up 2D rotor test problem.
// The center of the grid is assumed to have coordinates (x1,x2) = [0,0]; the grid
// initialization must be consistent with this
//
// REFERENCE: G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes", JCP, 161,
//   605 (2000)
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  Coordinates *pco = pmb->pcoord;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real gm1 = (pfl->pf_eos->GetGamma() - 1.0);

// Read initial conditions from 'athinput'

  Real v0  = pin->GetReal("problem","v0");
  Real p0  = pin->GetReal("problem","p0");
  Real bx0 = pin->GetReal("problem","bx0");
  Real r0  = pin->GetReal("problem","r0");
  Real r1  = pin->GetReal("problem","r1");

// Initialize the grid.  Note the center is always assumed to have coordinates
// x1=0, x2=0; the grid range in the input file must be consistent with this

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfl->u(IDN,k,j,i) = 1.0;
      pfl->u(IM1,k,j,i) = 0.0;
      pfl->u(IM2,k,j,i) = 0.0;
      pfl->u(IM3,k,j,i) = 0.0;

// reset density, velocity if cell is inside rotor

      Real rad = sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(j)));
      if (rad <= r0) {
        pfl->u(IDN,k,j,i) = 10.0;
        pfl->u(IM1,k,j,i) = -100.0*v0*pco->x2v(j);
        pfl->u(IM2,k,j,i) = 100.0*v0*pco->x1v(i);
      } else {

// smooth solution between r0 and r1.  For no smoothing, set r1<0 in input

        if (rad <= r1) {
          Real frac = (0.115 - rad)/(0.015);
          pfl->u(IDN,k,j,i) = 1.0 + 9.0*frac;
          pfl->u(IM1,k,j,i) = -frac*100.0*v0*pco->x2v(j);
          pfl->u(IM2,k,j,i) =  frac*100.0*v0*pco->x1v(i);
        }
      }

    }
  }}

// initialize interface B

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfd->b.x1f(k,j,i) = bx0;
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfd->b.x2f(k,j,i) = 0.0;
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfd->b.x3f(k,j,i) = 0.0;
  }}}

// initialize total energy

  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfl->u(IEN,k,j,i) = p0/gm1 + 0.5*bx0*bx0 +
          (SQR(pfl->u(IM1,k,j,i)) + SQR(pfl->u(IM2,k,j,i)))/pfl->u(IDN,k,j,i);
      }
    }}
  }

  return;
}
