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
#include "field.hpp"

// C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs(), sqrt()

// Athena headers
#include "../athena.hpp"                  // array access, macros, Real
#include "../athena_arrays.hpp"           // AthenaArray
#include "../mesh.hpp"                    // MeshBlock, Mesh
#include "integrators/field_integrator.hpp"  // FieldIntegrator

//======================================================================================
//! \file field.cpp
//  \brief implementation of functions in class Field
//======================================================================================

// constructor, initializes data structures and parameters

Field::Field(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_mblock = pmb;

// Allocate memory for interface fields, but only when needed.

  if (MAGNETIC_FIELDS_ENABLED) {
    int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
    if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

//  Note the extra cell in each longitudinal dirn for interface fields

    b.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    b1.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b1.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b1.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    bcc.NewAthenaArray (NFIELD,ncells3,ncells2,ncells1);
    bcc1.NewAthenaArray(NFIELD,ncells3,ncells2,ncells1);

    e.x1f.NewAthenaArray((NFIELDM1), ncells3   , ncells2   ,(ncells1+1));
    e.x2f.NewAthenaArray((NFIELDM1), ncells3   ,(ncells2+1), ncells1   );
    e.x3f.NewAthenaArray((NFIELDM1),(ncells3+1), ncells2   , ncells1   );

    wght.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    wght.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    wght.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    e1.NewAthenaArray((ncells3+1),(ncells2+1), ncells1   );
    e2.NewAthenaArray((ncells3+1), ncells2   ,(ncells1+1));
    e3.NewAthenaArray( ncells3   ,(ncells2+1),(ncells1+1));

// Construct ptrs to objects of various classes needed to integrate B-field

    pint = new FieldIntegrator(this, pin);

  }
}

// destructor

Field::~Field()
{
  b.x1f.DeleteAthenaArray();
  b.x2f.DeleteAthenaArray();
  b.x3f.DeleteAthenaArray();
  b1.x1f.DeleteAthenaArray();
  b1.x2f.DeleteAthenaArray();
  b1.x3f.DeleteAthenaArray();
  bcc.DeleteAthenaArray();
  bcc1.DeleteAthenaArray();

  e.x1f.DeleteAthenaArray();
  e.x2f.DeleteAthenaArray();
  e.x3f.DeleteAthenaArray();
  wght.x1f.DeleteAthenaArray();
  wght.x2f.DeleteAthenaArray();
  wght.x3f.DeleteAthenaArray();
  e1.DeleteAthenaArray();
  e2.DeleteAthenaArray();
  e3.DeleteAthenaArray();
}
