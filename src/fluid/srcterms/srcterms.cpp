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
 * distribution.  If not see <http://www.gnu.org/licenses/>.
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
  pt_mass_ = pin->GetOrAddReal("problem","GM",0.0);

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
  pmy_fluid_ = NULL; // Fluid destructor will free this memory
  volume_i_.DeleteAthenaArray();
  src_terms_i_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief
 */

void FluidSourceTerms::PhysicalSourceTerms(Real dt, AthenaArray<Real> &prim,
  AthenaArray<Real> &cons)
{
  Real src[NVAR];
  MeshBlock *pmb = pmy_fluid_->pmy_block;
  if (pt_mass_ == 0.0) return;

// src = GM/R_{cell_center}

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma simd
    for (int i=pmb->is; i<=pmb->ie; ++i) {
      src[IM1] = src_terms_i_(i)*(pt_mass_*prim(IDN,k,j,i)/(pmb->x1v(i)));

      Real& uim1 = cons(IM1,k,j,i);
      uim1 -= dt*src[IM1];
    }
  }}

  return;
}
