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
//! \file blast.cpp
//  \brief Problem generator for spherical blast wave problem.
//
// REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for 
//   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

#include <cmath>

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real rout = pin->GetReal("problem","radius");
  Real rin = rout - pin->GetOrAddReal("problem","ramp",0.0);
  Real pa  = pin->GetReal("problem","pamb");
  Real da  = pin->GetOrAddReal("problem","damb",1.0);
  Real prat = pin->GetReal("problem","prat");
  Real drat = pin->GetOrAddReal("problem","drat",1.0);
  Real b0,theta;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem","b0");
    theta = (PI/180.0)*pin->GetReal("problem","angle");
  }
  Real gamma = phydro->peos->GetGamma();
  Real gm1 = gamma - 1.0;

// setup uniform ambient medium with spherical over-pressured region

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rad = sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)) + SQR(pcoord->x3v(k)));
    Real den = da;
    if (rad < rout) {
      if (rad < rin) {
        den = drat*da;
      } else {
        Real f = (rad-rin) / (rout-rin);
        Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
        den = std::exp(log_den);
      }
    }

    phydro->u(IDN,k,j,i) = den;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
    if (NON_BAROTROPIC_EOS) {
      Real pres = pa;
      if (rad < rout) {
        if (rad < rin) {
          pres = prat*pa;
        } else {
          Real f = (rad-rin) / (rout-rin);
          Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
          pres = std::exp(log_pres);
        }
      }
      phydro->u(IEN,k,j,i) = pres/gm1;
      if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
        phydro->u(IEN,k,j,i) += den;
    }
  }}}

// initialize interface B and total energy

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = b0*cos(theta);
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = b0*sin(theta);
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      phydro->u(IEN,k,j,i) += 0.5*b0*b0;
    }}}
  }

}
