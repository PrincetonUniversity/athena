//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rotor.c
//  \brief Sets up 2D rotor test problem.
// The center of the grid is assumed to have coordinates (x1,x2) = [0,0]; the grid
// initialization must be consistent with this
//
// REFERENCE: G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes", JCP, 161,
//   605 (2000)
//========================================================================================

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Rotor test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gm1 = peos->GetGamma() - 1.0;

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
      phydro->u(IDN,k,j,i) = 1.0;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;

      // reset density, velocity if cell is inside rotor
      Real rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)));
      if (rad <= r0) {
        phydro->u(IDN,k,j,i) = 10.0;
        phydro->u(IM1,k,j,i) = -100.0*v0*pcoord->x2v(j);
        phydro->u(IM2,k,j,i) = 100.0*v0*pcoord->x1v(i);
      } else {

        // smooth solution between r0 and r1.  For no smoothing, set r1<0 in input
        if (rad <= r1) {
          Real frac = (0.115 - rad)/(0.015);
          phydro->u(IDN,k,j,i) = 1.0 + 9.0*frac;
          phydro->u(IM1,k,j,i) = -frac*100.0*v0*pcoord->x2v(j);
          phydro->u(IM2,k,j,i) =  frac*100.0*v0*pcoord->x1v(i);
        }
      }

    }
  }}

  // initialize interface B
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfield->b.x1f(k,j,i) = bx0;
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = 0.0;
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x3f(k,j,i) = 0.0;
  }}}

  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IEN,k,j,i) = p0/gm1 + 0.5*bx0*bx0 +
          (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }}
  }

  return;
}
