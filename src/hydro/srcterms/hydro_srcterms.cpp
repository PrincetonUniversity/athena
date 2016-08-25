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
//======================================================================================
//! \file hydro_srcterms.cpp
//  \brief Class to implement source terms in the hydro equations
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../parameter_input.hpp"

// this class header
#include "hydro_srcterms.hpp"

#include <iostream>

// HydroSourceTerms constructor

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin)
{
  pmy_hydro_ = phyd;
  hydro_sourceterms_defined = false;

  // read point mass or constant acceleration parameters from input block
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  if (gm_ != 0.0) hydro_sourceterms_defined = true;

  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  if (g1_ != 0.0) hydro_sourceterms_defined = true;

  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  if (g2_ != 0.0) hydro_sourceterms_defined = true;

  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  if (g3_ != 0.0) hydro_sourceterms_defined = true;
//[JMSHI
  // read shearing box parameters from input block
  Omega_0_ = pin->GetOrAddReal("problem","Omega0",0.0);
  qshear_  = pin->GetOrAddReal("problem","qshear",0.0);
  ShBoxCoord_ = pin->GetOrAddInteger("problem","shboxcoord",1);
  if ((Omega_0_ !=0.0) && (qshear_ != 0.0)) hydro_sourceterms_defined = true;
//JMSHI]

  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
  if(UserSourceTerm != NULL) hydro_sourceterms_defined = true;
}

// destructor

HydroSourceTerms::~HydroSourceTerms()
{
}

//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::AddHydroSourceTerms
//  \brief Adds source terms to conserved variables

void HydroSourceTerms::AddHydroSourceTerms(const Real time, const Real dt,
     const AthenaArray<Real> *flux, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // accleration due to point mass (MUST BE AT ORIGIN)
  if (gm_ != 0.0) PointMass(dt, flux, prim, cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0) ConstantAcceleration(dt, flux, prim,cons);

  // Add new source terms here
  // MyNewSourceTerms()
//[JMSHI
  // shearing box source terms: tidal and Coriolis forces
  if ((Omega_0_ !=0.0) && (qshear_ != 0.0)) ShearingBoxSourceTerms(dt, flux, prim, cons);
//JMSHI]

  //  user-defined source terms
  if (UserSourceTerm != NULL)
    UserSourceTerm(pmb, time,dt,prim,bcc,cons);

  return;
}

