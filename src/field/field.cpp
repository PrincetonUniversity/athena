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
#include "../coordinates/coordinates.hpp" // Coordinates

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

    e.x1e.NewAthenaArray((ncells3+1),(ncells2+1), ncells1   );
    e.x2e.NewAthenaArray((ncells3+1), ncells2   ,(ncells1+1));
    e.x3e.NewAthenaArray( ncells3   ,(ncells2+1),(ncells1+1));

    ei.x1f.NewAthenaArray(((NFIELD)-1), ncells3   , ncells2   ,(ncells1+1));
    ei.x2f.NewAthenaArray(((NFIELD)-1), ncells3   ,(ncells2+1), ncells1   );
    ei.x3f.NewAthenaArray(((NFIELD)-1),(ncells3+1), ncells2   , ncells1   );
    wght.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    wght.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    wght.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

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

  e.x1e.DeleteAthenaArray();
  e.x2e.DeleteAthenaArray();
  e.x3e.DeleteAthenaArray();
  ei.x1f.DeleteAthenaArray();
  ei.x2f.DeleteAthenaArray();
  ei.x3f.DeleteAthenaArray();
  wght.x1f.DeleteAthenaArray();
  wght.x2f.DeleteAthenaArray();
  wght.x3f.DeleteAthenaArray();
}


void Field::CalculateCellCenteredField(const InterfaceField &bf, AthenaArray<Real> &bc,
            Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  int nthreads = pmy_mblock->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
    for (int j=js; j<=je; ++j){
    // calc cell centered fields first
#pragma simd
      for (int i=is; i<=ie; ++i){
        const Real& b1_i   = bf.x1f(k,j,i  );
        const Real& b1_ip1 = bf.x1f(k,j,i+1);
        const Real& b2_j   = bf.x2f(k,j  ,i);
        const Real& b2_jp1 = bf.x2f(k,j+1,i);
        const Real& b3_k   = bf.x3f(k  ,j,i);
        const Real& b3_kp1 = bf.x3f(k+1,j,i);

        Real& bcc1 = bc(IB1,k,j,i);
        Real& bcc2 = bc(IB2,k,j,i);
        Real& bcc3 = bc(IB3,k,j,i);

        // cell center B-fields are defined as spatial interpolation at the volume center
        const Real& x1f_i  = pco->x1f(i);
        const Real& x1f_ip = pco->x1f(i+1);
        const Real& x1v_i  = pco->x1v(i);
        const Real& dx1_i  = pco->dx1f(i);
        Real lw=(x1f_ip-x1v_i)/dx1_i;
        Real rw=(x1v_i -x1f_i)/dx1_i;
        bcc1 = lw*b1_i + rw*b1_ip1;
        const Real& x2f_j  = pco->x2f(j);
        const Real& x2f_jp = pco->x2f(j+1);
        const Real& x2v_j  = pco->x2v(j);
        const Real& dx2_j  = pco->dx2f(j);
        lw=(x2f_jp-x2v_j)/dx2_j;
        rw=(x2v_j -x2f_j)/dx2_j;
        bcc2 = lw*b2_j + rw*b2_jp1;
        const Real& x3f_k  = pco->x3f(k);
        const Real& x3f_kp = pco->x3f(k+1);
        const Real& x3v_k  = pco->x3v(k);
        const Real& dx3_k  = pco->dx3f(k);
        lw=(x3f_kp-x3v_k)/dx3_k;
        rw=(x3v_k -x3f_k)/dx3_k;
        bcc3 = lw*b3_k + rw*b3_kp1;
      }
    }
  }
}
  return;
}

