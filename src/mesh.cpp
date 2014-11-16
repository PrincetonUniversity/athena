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
#include "mesh.hpp"

// C++ headers
#include <cfloat>     // FLT_MAX
#include <cmath>      // std::abs(), pow()
#include <iostream>   // cout, endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "athena.hpp"                   // enums, macros, Real
#include "athena_arrays.hpp"            // AthenaArray
#include "coordinates/coordinates.hpp"  // Coordinates
#include "fluid/fluid.hpp"              // Fluid
#include "field/field.hpp"              // Field
#include "bvals/bvals.hpp"              // BoundaryValues
#include "fluid/eos/eos.hpp"              // FluidEqnOfState
#include "fluid/integrators/fluid_integrator.hpp"  // FluidIntegrator
#include "field/integrators/field_integrator.hpp"  // FieldIntegrator
#include "parameter_input.hpp"          // ParameterInput

//======================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in classes Mesh, MeshDomain, and MeshBlock
//======================================================================================

//--------------------------------------------------------------------------------------
// Mesh constructor, builds mesh at start of calculation using parameters in input file

Mesh::Mesh(ParameterInput *pin)
{
  std::stringstream msg;

// read time and cycle limits from input file

  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  cfl_number = pin->GetReal("time","cfl_number");
  time = start_time;
  dt   = (FLT_MAX);

  nlim = pin->GetOrAddInteger("time","nlim",-1);
  ncycle = 0;

// read number of OpenMP threads for mesh

  nthreads_mesh = pin->GetOrAddReal("mesh","max_num_threads",1);
  if (nthreads_mesh < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but max_num_threads=" 
        << nthreads_mesh << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// read number of grid cells in mesh (root domain) from input file.  

  mesh_size.nx1 = pin->GetInteger("mesh","nx1");
  if (mesh_size.nx1 < 4) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx1 must be >= 4, but nx1=" 
        << mesh_size.nx1 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  mesh_size.nx2 = pin->GetInteger("mesh","nx2");
  if (mesh_size.nx2 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx2 must be >= 1, but nx2=" 
        << mesh_size.nx2 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  mesh_size.nx3 = pin->GetInteger("mesh","nx3");
  if (mesh_size.nx3 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx3 must be >= 1, but nx3=" 
        << mesh_size.nx3 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.nx2 == 1 && mesh_size.nx3 > 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file: nx2=1, nx3=" << mesh_size.nx3 
        << ", 2D problems in x1-x3 plane not supported" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// read physical size of mesh (root domain) from input file.  

  mesh_size.x1min = pin->GetReal("mesh","x1min");
  mesh_size.x2min = pin->GetReal("mesh","x2min");
  mesh_size.x3min = pin->GetReal("mesh","x3min");

  mesh_size.x1max = pin->GetReal("mesh","x1max");
  mesh_size.x2max = pin->GetReal("mesh","x2max");
  mesh_size.x3max = pin->GetReal("mesh","x3max");

  if (mesh_size.x1max <= mesh_size.x1min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x1max must be larger than x1min: x1min=" << mesh_size.x1min 
        << " x1max=" << mesh_size.x1max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.x2max <= mesh_size.x2min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x2max must be larger than x2min: x2min=" << mesh_size.x2min 
        << " x2max=" << mesh_size.x2max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.x3max <= mesh_size.x3min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x3max must be larger than x3min: x3min=" << mesh_size.x3min 
        << " x3max=" << mesh_size.x3max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// read ratios of grid cell size in each direction

  mesh_size.x1rat = pin->GetOrAddReal("mesh","x1rat",1.0);
  mesh_size.x2rat = pin->GetOrAddReal("mesh","x2rat",1.0);
  mesh_size.x3rat = pin->GetOrAddReal("mesh","x3rat",1.0);

  if (std::abs(mesh_size.x1rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x1rat <= 1.1, x1rat=" 
        << mesh_size.x1rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (std::abs(mesh_size.x2rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x2rat <= 1.1, x2rat=" 
        << mesh_size.x2rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (std::abs(mesh_size.x3rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x3rat <= 1.1, x3rat=" 
        << mesh_size.x3rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// read BC flags for each of the 6 boundaries in turn.  Error tests performed in
// BoundaryValues constructor

  mesh_bcs.ix1_bc = pin->GetOrAddInteger("mesh","ix1_bc",0);
  mesh_bcs.ox1_bc = pin->GetOrAddInteger("mesh","ox1_bc",0);
  mesh_bcs.ix2_bc = pin->GetOrAddInteger("mesh","ix2_bc",0);
  mesh_bcs.ox2_bc = pin->GetOrAddInteger("mesh","ox2_bc",0);
  mesh_bcs.ix3_bc = pin->GetOrAddInteger("mesh","ix3_bc",0);
  mesh_bcs.ox3_bc = pin->GetOrAddInteger("mesh","ox3_bc",0);

// allocate root domain

  pdomain = new MeshDomain(mesh_size, mesh_bcs, this, pin);
}

// destructor

Mesh::~Mesh()
{
  delete pdomain;
}

//--------------------------------------------------------------------------------------
// MeshDomain constructor: builds array of MeshBlocks based on input arguments.

MeshDomain::MeshDomain(RegionSize in_size, RegionBCs in_bcs, Mesh* pm, ParameterInput *pin)
{
  pmy_mesh = pm;
  domain_size = in_size;
  domain_bcs  = in_bcs;

// allocate block on this domain
// In future w MPI: calculate array of blocks, set their region sizes, and initialize

  pblock = new MeshBlock(domain_size, domain_bcs, this, pin);

  return;
}

// destructor

MeshDomain::~MeshDomain()
{
  delete pblock;
}

//--------------------------------------------------------------------------------------
// MeshBlock constructor: builds 1D vectors of cell positions and spacings, and
// constructs coordinate, boundary condition, fluid and field objects.

MeshBlock::MeshBlock(RegionSize in_size, RegionBCs in_bcs, MeshDomain *pd,
  ParameterInput *pin)
{
  pmy_domain = pd;
  block_size = in_size;
  block_bcs  = in_bcs;

// initialize grid indices

  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  if (block_size.nx2 > 1) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
  } else {
    js = je = 0;
  }

  if (block_size.nx3 > 1) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
  } else {
    ks = ke = 0;
  }

  std::cout << "is=" << is << " ie=" << ie << std::endl;
  std::cout << "js=" << js << " je=" << je << std::endl;
  std::cout << "ks=" << ks << " ke=" << ke << std::endl;

// allocate arrays for sizes and positions of cells

  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

// cell sizes
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

// cell positions. Note the extra element for cell face positions
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

// X1-DIRECTION: initialize sizes and positions of cell FACES (dx1f,x1f)

  Real dx = block_size.x1max - block_size.x1min;
  if (block_size.x1rat == 1.0) {
    dx1f(is) = dx/(Real)block_size.nx1;
  } else {
    dx1f(is) = dx*(block_size.x1rat - 1.0)/(pow(block_size.x1rat,block_size.nx1) - 1.0);
  }

  x1f(is) = block_size.x1min;
  for (int i=is+1; i<=ie; ++i) {
    x1f(i) = x1f(i-1) + dx1f(i-1);
    dx1f(i) = dx1f(i-1)*block_size.x1rat;
  }
  x1f(ie+1) = block_size.x1max;

// cell face positions and spacing in ghost zones

  if (block_bcs.ix1_bc == 1) {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  } else {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(is-i) = dx1f(is-i+1)/block_size.x1rat;
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }

  if (block_bcs.ox1_bc == 1) {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  } else {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(ie+i  ) = dx1f(ie+i-1)*block_size.x1rat;
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

// X2-DIRECTION: initialize spacing and positions of cell FACES (dx2f,x2f)

  dx = block_size.x2max - block_size.x2min;
  if (block_size.nx2 == 1) {
    dx2f(js) = dx;
    x2f(js  ) = block_size.x2min;
    x2f(je+1) = block_size.x2max;
  } else {
    if (block_size.x2rat == 1.0) {
      dx2f(js) = dx/(Real)block_size.nx2;
    } else {
      dx2f(js)=dx*(block_size.x2rat - 1.0)/(pow(block_size.x2rat,block_size.nx2) - 1.0);
    }

    x2f(js) = block_size.x2min;
    for (int j=js+1; j<=je; ++j) {
      x2f(j) = x2f(j-1) + dx2f(j-1);
      dx2f(j) = dx2f(j-1)*block_size.x2rat;
    }
    x2f(je+1) = block_size.x2max;

// cell face positions and spacing in ghost zones

    if (block_bcs.ix2_bc == 1) {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    } else {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(js-j) = dx2f(js-j+1)/block_size.x2rat;
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }

    if (block_bcs.ox2_bc == 1) {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(je+j  ) = dx2f(je-j+1);
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    } else {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(je+j  ) = dx2f(je+j-1)*block_size.x2rat;
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    }
  }

// X3-DIRECTION: initialize spacing and positions of cell FACES (dx3f,x3f)

  dx = block_size.x3max - block_size.x3min;
  if (block_size.nx3 == 1) {
    dx3f(ks) = dx;
    x3f(ks  ) = block_size.x3min;
    x3f(ke+1) = block_size.x3max;
  } else {
    if (block_size.x3rat == 1.0) {
      dx3f(ks) = dx/(Real)block_size.nx3;
    } else {
      dx3f(ks)=dx*(block_size.x3rat - 1.0)/(pow(block_size.x3rat,block_size.nx3) - 1.0);
    }

    x3f(ks) = block_size.x3min;
    for (int k=ks+1; k<=ke; ++k) {
      x3f(k) = x3f(k-1) + dx3f(k-1);
      dx3f(k) = dx3f(k-1)*block_size.x3rat;
    }
    x3f(ke+1) = block_size.x3max;

// cell face positions and spacing in ghost zones

    if (block_bcs.ix3_bc == 1) {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    } else {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ks-k) = dx3f(ks-k+1)/block_size.x3rat;
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }

    if (block_bcs.ox3_bc == 1) {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ke+k  ) = dx3f(ke-k+1);
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    } else {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ke+k  ) = dx3f(ke+k-1)*block_size.x3rat;
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    }
  }

// construct Coordinates and Fluid objects stored in MeshBlock class.  Note that the
// initial conditions for the fluid are set in problem generator called from main, not
// in the Fluid constructor
 
  pcoord = new Coordinates(this, pin);
  pfluid = new Fluid(this, pin);
  pfield = new Field(this, pin);
  pbval  = new BoundaryValues(this, pin);

  return;
}

// destructor

MeshBlock::~MeshBlock()
{
  dx1f.DeleteAthenaArray();  
  dx2f.DeleteAthenaArray();  
  dx3f.DeleteAthenaArray();  
  dx1v.DeleteAthenaArray();  
  dx2v.DeleteAthenaArray();  
  dx3v.DeleteAthenaArray();  
  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();

  delete pcoord;
  delete pfluid;
  delete pfield;
  delete pbval;
}

//--------------------------------------------------------------------------------------
// \!fn void Mesh::ForAllDomains(enum ActionOnDomain action, ParameterInput *pin)
// \brief function that loops over all MeshDomains and calls MeshBlock functions

void Mesh::ForAllDomains(enum ActionOnDomain action, ParameterInput *pin)
{

// Eventually this will be a loop over all domains

  MeshBlock* pmb = pdomain->pblock;
  if (pmb != NULL)  {
    Fluid* pfluid = pdomain->pblock->pfluid;
    Field* pfield = pdomain->pblock->pfield;

    switch (action) {

      case pgen: // call problem generator
        ProblemGenerator(pfluid,pfield,pin);
        break;

      case fluid_bcs_n: // set fluid BCs at t^n
        pmb->pbval->ApplyBVals(pfluid->u);
        break;

      case fluid_bcs_nhalf: // set fluid BCs at t^{intermediate}
        pmb->pbval->ApplyBVals(pfluid->u1);
        break;

      case bfield_bcs_n: // set bfield BCs at t^n
        pmb->pbval->ApplyBVals(pfield->b);
        break;

      case bfield_bcs_nhalf: // set bfield BCs at t^{intermediate}
        pmb->pbval->ApplyBVals(pfield->b1);
        break;

      case fluid_predict: // integrate fluid to intermediate step 
        pfluid->pf_integrator->Predict(pdomain->pblock);
        break;

      case fluid_correct: // integrate fluid for full timestep, t^n --> t^{n+1}
        pfluid->pf_integrator->Correct(pdomain->pblock);
        break;

      case bfield_predict: // integrate fluid to intermediate step 
        pfield->pint->CT(pmb, pfield->b, pfield->b1, 0.5*dt);
        break;

      case bfield_correct: // integrate fluid for full timestep, t^n --> t^{n+1}
        pfield->pint->CT(pmb, pfield->b, pfield->b, dt);
        break;

      case primitives_n: // compute primitives from conseerved at t^n
        pfluid->pf_eos->ConservedToPrimitive(pfluid->u,pfluid->w,pfluid->w1,
           pfield->b, pfield->bcc);
        break;

      case primitives_nhalf: // compute primitives from conseerved at t^{intermediate}
        pfluid->pf_eos->ConservedToPrimitive(pfluid->u1,pfluid->w1,pfluid->w,
           pfield->b1, pfield->bcc1);
        break;

      case new_timestep: // calculate new time step
        pfluid->NewTimeStep(pdomain->pblock);
        break;

    }
  }
}
