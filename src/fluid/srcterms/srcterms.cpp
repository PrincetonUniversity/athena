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

 *====================================================================================*/

// Primary header
#include "srcterms.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../../athena.hpp"          // Real
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../../mesh.hpp"            // MeshBlock
#include "../../coordinates/coordinates.hpp"  // VectorBetweenPoints()
#include "../fluid.hpp"              // Fluid
#include "../../parameter_input.hpp" // ParameterInput

//======================================================================================
/*! \file srcterms.cpp
 *  \brief implements functions that compute physical source terms in the fluid
 *====================================================================================*/

// FluidSourceTerms constructor - sets function pointers for each of the physical source
// terms to be included in the calculation.

FluidSourceTerms::FluidSourceTerms(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;

  pfirst_mass = NULL;
  PointMass *pnew_mass;
  PointMass *plast = pfirst_mass;

  InputBlock *pib = pin->pfirst_block;
  while (pib != NULL) {
    if (pib->block_name.compare(0,9,"pointmass") == 0) {
      pnew_mass = new PointMass;
      pnew_mass->gm = pin->GetReal(pib->block_name,"GM");
      pnew_mass->position.x1 = pin->GetReal(pib->block_name,"x1_0");
      pnew_mass->position.x2 = pin->GetReal(pib->block_name,"x2_0");
      pnew_mass->position.x3 = pin->GetReal(pib->block_name,"x3_0");
      pnew_mass->velocity.x1 = pin->GetOrAddReal(pib->block_name,"v1_0",0.0);
      pnew_mass->velocity.x2 = pin->GetOrAddReal(pib->block_name,"v2_0",0.0);
      pnew_mass->velocity.x3 = pin->GetOrAddReal(pib->block_name,"v3_0",0.0);
      pnew_mass->pnext = NULL;
      if (pfirst_mass == NULL) {
        pfirst_mass = pnew_mass;
      } else {
        plast->pnext = pnew_mass;
      }
      plast = pnew_mass;
    }
    pib = pib->pnext;
  }

// Allocate memory for scratch arrays used in integrator, and internal scratch arrays
// Only allocate arrays needed for forces in x1 direction for now 

  MeshBlock *pmb = pmy_fluid_->pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  volume_i_.NewAthenaArray(ncells1);
  src_terms_i_.NewAthenaArray(ncells1);

// compute constant factors used to compute face-areas and cell volumes and store in
// local scratch arrays.  This helps improve performance.

#pragma simd
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
    volume_i_(i)    = 0.5*(pmb->x1f(i+1)*pmb->x1f(i+1) - pmb->x1f(i)*pmb->x1f(i));
    src_terms_i_(i) = pmb->dx1f(i)/volume_i_(i);
  }
}

// destructor

FluidSourceTerms::~FluidSourceTerms()
{
// iterate through linked list of point masses and delete each node
  PointMass *ppm = pfirst_mass;
  while (ppm != NULL) {
    PointMass *pold_mass = ppm;
    ppm = ppm->pnext;
    delete pold_mass;
  }

  volume_i_.DeleteAthenaArray();
  src_terms_i_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief
 */

void FluidSourceTerms::PhysicalSourceTerms(const Real dt, const AthenaArray<Real> &prim,
  AthenaArray<Real> &cons)
{
  if (pfirst_mass == NULL) return;

  MeshBlock *pmb = pmy_fluid_->pmy_block;
  PointMass *ppm = pfirst_mass;
  ThreeVector r,p1;

// Source terms due to point mass gravity

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
    while (ppm != NULL) {
      p1.x2 = pmb->x2v(j);
      p1.x3 = pmb->x3v(k);
#pragma simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        p1.x1 = pmb->x1v(i);
        r = pmb->pcoord->VectorBetweenPoints(p1,ppm->position);
        Real d = sqrt(r.x1*r.x1 + r.x2*r.x2 + r.x3*r.x3);
        Real force = (ppm->gm)*prim(IDN,k,j,i)/(d*d);

        Real src = dt*src_terms_i_(i)*pmb->x1v(i)*force*(r.x1/d);
        cons(IM1,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVX,k,j,i);

if (i==2 && j==2 && k==0) {
std::cout << r.x1 << "  " << r.x2 << "  " << r.x3 << std::endl;
std::cout << src/dt << "  " << d << std::endl; 
}
 

        if (pmb->block_size.nx2 > 1) {
          src = dt*force*(r.x2/d);
          cons(IM2,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);
if (i==2 && j==2 && k==0) {
std::cout << src/dt << "  " << d << std::endl; 
}
 
        }

        if (pmb->block_size.nx3 > 1) {
          src = dt*force*(r.x3/d);
          cons(IM3,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVZ,k,j,i);
if (i==2 && j==2 && k==0) {
std::cout << src/dt << "  " << d << std::endl; 
}
 
        }
      }
      ppm = ppm->pnext;
    }
  }}

  return;
}
