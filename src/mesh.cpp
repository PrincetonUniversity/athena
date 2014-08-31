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
#include "fluid/fluid.hpp"                    // Fluid
#include "fluid/bvals/bvals.hpp"              // FluidBoundaryConditions
#include "fluid/eos/eos.hpp"              // FluidEqnOfState
#include "fluid/integrators/integrators.hpp"  // FluidIntegrator
#include "outputs/outputs.hpp"          // Outputs
#include "parameter_input.hpp"          // ParameterInput

//======================================================================================
/*! \file mesh.cpp
 *  \brief implementation of functions in classes Mesh, MeshDomain, and MeshBlock
 *====================================================================================*/

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
// function FluidBoundaryConditions::InitBoundaryConditions

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
// constructs coordinate, boundary condition, and fluid objects.

MeshBlock::MeshBlock(RegionSize in_size, RegionBCs in_bcs, MeshDomain *pd, ParameterInput *pin)
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

// allocate arrays for spacing and positions of cell faces and volume centers

  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

  x1f.NewAthenaArray(ncells1 + 1);
  x2f.NewAthenaArray(ncells2 + 1);
  x3f.NewAthenaArray(ncells3 + 1);
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

// X1-DIRECTION: initialize spacing and positions of cell FACES (dx1f,x1f)

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

// construct various objects stored in MeshBlock class.  Constructors do the following:
//   Coordinates: initializes volume-centered coordinates (x1v,dx1v,...)
//   Fluid: allocates memory for arrays, and constructs FluidIntegrator.
//   Outputs:
// initial conditions for fluid set in problem generator called from main
 
  pcoord   = new Coordinates(this, pin);
  pfluid   = new Fluid(this, pin);
  poutputs = new Outputs(this, pin);

/*******************
  for (int i=0; i<((ie-is+1)+2*(NGHOST)); ++i) {
    printf("i=%i  x1f=%e  dx1f=%e x1v=%e dx1v=%e \n",i,x1f(i),dx1f(i),x1v(i),dx1v(i));
  }
  printf("i=%i  x1f=%e  \n",((ie-is+1)+2*NGHOST),x1f(((ie-is+1)+2*NGHOST)));

  if (block_size.nx2 == 1) {
    printf("j=%i  x2f=%e  dx2f=%e x2v=%e dx2v=%e\n",js,x2f(js),dx2f(js),x2v(js),dx2v(js));
  } else {
    for (int j=0; j<((je-js+1)+2*(NGHOST)); ++j) {
      printf("j=%i  x2f=%e  dx2f=%e x2v=%e dx2v=%e \n",j,x2f(j),dx2f(j),x2v(j),dx2v(j));
    }
    printf("j=%i  x2f=%e  \n",((je-js+1)+2*NGHOST),x2f(((je-js+1)+2*NGHOST)));
  }

  if (block_size.nx3 == 1) {
    printf("k=%i  x3f=%e  dx3f=%e x3v=%e dx3v=%e\n",ks,x3f(ks),dx3f(ks),x3v(ks),dx3v(ks));
  } else {
    for (int k=0; k<((ke-ks+1)+2*(NGHOST)); ++k) {
      printf("k=%i  x3f=%e  dx3f=%e x3v=%e dx3v=%e \n",k,x3f(k),dx3f(k),x3v(k),dx3v(k));
    }
  printf("k=%i  x3f= %e  \n",((ke-ks+1)+2*NGHOST),x3f(((ke-ks+1)+2*NGHOST)));
  }
*******************/
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
}

//--------------------------------------------------------------------------------------
// \!fn void Mesh::ForAllDomains(enum ActionOnDomain action, ParameterInput *pin)
// \brief function that loops over all MeshDomains and calls MeshBlock functions

void Mesh::ForAllDomains(enum ActionOnDomain action, ParameterInput *pin)
{

// Eventually this will be a loop over all domains

  if (pdomain->pblock != NULL)  {
    Fluid* pf = pdomain->pblock->pfluid;

    switch (action) {

      case init_fluid:
        pdomain->pblock->pfluid->InitFluid(pin);
        break;

      case init_outputs:
        pdomain->pblock->poutputs->InitOutputTypes(pin);
        break;

      case fluid_bcs_n:
        pdomain->pblock->pfluid->pf_bcs->ApplyBoundaryConditions(pf->u);
        break;

      case fluid_bcs_nhalf:
        pdomain->pblock->pfluid->pf_bcs->ApplyBoundaryConditions(pf->u1);
        break;

      case fluid_predict:
        pdomain->pblock->pfluid->pf_integrator->Predict(pdomain->pblock);
        break;

      case fluid_correct:
        pdomain->pblock->pfluid->pf_integrator->Correct(pdomain->pblock);
        break;

      case primitives_n:
        pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(pf->u,pf->w1,pf->w);
        break;

      case primitives_nhalf:
        pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(pf->u1,pf->w,pf->w1);
        break;

      case new_timestep:
        pdomain->pblock->pfluid->NewTimeStep(pdomain->pblock);
        break;

      case make_outputs:
        pdomain->pblock->poutputs->MakeOutputs();
        break;
    }
  }
}
