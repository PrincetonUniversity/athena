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
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>

#include "athena.hpp"
#include "athena_arrays.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "fluid.hpp"

//======================================================================================
//! \file fluid.cpp
//  \brief implementation of functions in class Fluid
//======================================================================================

// constructor, initializes data structures and parameters, calls problem generator

Fluid::Fluid(ParameterInput *pin, Mesh *pm)
{
// Read adiabatic index from input file

  gamma_ = pin->GetReal("fluid","gamma");
  gamma_m1_ = gamma_ - 1.0;

// Allocate FluidData structure on root Domain, initialize times

  proot = new FluidData;
  proot->time = 0.0;
  proot->dt   = 1.0;

// Allocate memory for primitive/conserved variables on root Domain

  int ncells1 = pm->root.pblock->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pm->root.pblock->block_size.nx2 + 2*(NGHOST);
  int ncells3 = pm->root.pblock->block_size.nx3 + 2*(NGHOST);

  proot->u.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  proot->w.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

// Allocate memory for primitive/conserved variables at half-time step, and scratch

  u1_.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);
  w1_.NewAthenaArray(NVAR,ncells3,ncells2,ncells1);

  wl_.NewAthenaArray(NVAR,ncells1);
  wr_.NewAthenaArray(NVAR,ncells1);
  flx_.NewAthenaArray(NVAR,ncells1);

// call problem generator (set by file in pgen/ compiled by configure) on root domain

  ProblemGenerator(pin, &(pm->root));

// initialize integration algorithms, source terms, boundary conditions

}

// destructor

Fluid::~Fluid()
{
}

//--------------------------------------------------------------------------------------
/*! \fn  void Fluid::Predict()
 *  \brief integrates fluid over the predict step */

void Fluid::Predict(Mesh *pm)
{
  PredictVanLeer2(pm);

//  PredictSourceTerms
//  BoundaryValues
//  ConservedToPrimitive

  CorrectVanLeer2(pm);

//  CorrectSourceTerms
//  BoundaryValues
//  ConservedToPrimitive
}
