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

#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "athena.hpp"
#include "athena_arrays.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "fluid.hpp"
#include "bvals/bvals.hpp"
#include "convert_var/convert_var.hpp"
#include "integrators/integrators.hpp"
#include "geometry/geometry.hpp"

//======================================================================================
/*! \file mesh.cpp
 *  \brief implementation of functions in classes Mesh, Domain, and Block
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

  if (abs(mesh_size.x1rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x1rat <= 1.1, x1rat=" 
        << mesh_size.x1rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (abs(mesh_size.x2rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x2rat <= 1.1, x2rat=" 
        << mesh_size.x2rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (abs(mesh_size.x3rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x3rat <= 1.1, x3rat=" 
        << mesh_size.x3rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// allocate root domain

  pdomain = new Domain(mesh_size, this);
}

// Mesh destructor

Mesh::~Mesh()
{
  delete pdomain;
}

//--------------------------------------------------------------------------------------
// Domain constructor: builds array of Blocks based on input arguments.  May be called
// at any time in a simulation, whenever AMR creates a new Domain

Domain::Domain(RegionSize dom_size, Mesh* pm)
{
  pparent_mesh = pm;
  domain_size = dom_size;

// allocate block on this domain
// In future w MPI: calculate array of blocks, set their region sizes, and initialize

  pblock = new Block(domain_size, this);

  return;
}

// Domain destructor

Domain::~Domain()
{
  delete pblock;
}

//--------------------------------------------------------------------------------------
// Block constructor: builds 1D vectors of cell sizes and spacing based on input
// arguments. Constructs Fluid and Geometry objects (but these are initialized with
// separate init functions).  May be called at any time in an AMR simulations. 

Block::Block(RegionSize blk_size, Domain *pd)
{
  pparent_domain = pd;
  block_size = blk_size;

// initilize grid indices

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

// allocate arrays for positions and spacing of cell faces and volume centers

  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

  x1f.NewAthenaArray(ncells1 + 1);
  x2f.NewAthenaArray(ncells2 + 1);
  x3f.NewAthenaArray(ncells3 + 1);

  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);

  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

// initialize positions and spacing of x1-cell faces.

  Real dx = block_size.x1max - block_size.x1min;
  x1f(is  ) = block_size.x1min;
  x1f(ie+1) = block_size.x1max;
  if (block_size.x1rat == 1.0) {
    dx1f(is) = dx/(Real)block_size.nx1;
  } else {
    dx1f(is) = dx*(block_size.x1rat - 1.0)/(pow(block_size.x1rat,block_size.nx1) - 1.0);
  }
  for (int i=is+1; i<=ie; ++i) {
    dx1f(i) = block_size.x1rat*dx1f(i-1);
     x1f(i) = x1f(i-1) + dx1f(i-1);
  }
// compute x1-positions starting from x1max, then average with above to prevent
// round-off from accumulating in last cell.  Recompute dx1f for consistency.
  for (int i=ie; i>is; --i) {
    x1v(i) = x1f(i+1) - dx1f(i);  // x1v is being used as a temporary variable!!
  }
  for (int i=is+1; i<=ie; ++i) {
     x1f(i) = 0.5*(x1f(i) + x1v(i));
    dx1f(i) = x1f(i) - x1f(i-1);
  }
// cell face face positions and spacing in ghost zones
  for (int i=1; i<=(NGHOST); ++i) {
    dx1f(is-i) = dx1f(is-i+1)/block_size.x1rat;
     x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
  }
  for (int i=1; i<=(NGHOST); ++i) {
    dx1f(ie+i  ) = dx1f(ie+i-1)*block_size.x1rat;
     x1f(ie+i+1) =  x1f(ie+i) - dx1f(ie+i);
  }

/********************/
  for (int i=0; i<((ie-is+1)+2*(NGHOST)); ++i) {
    printf("i=%i  x1f= %e  dx1f=%e \n",i,x1f(i),dx1f(i));
  }
/********************/



  for (int i=0; i<ncells1; ++i) {
    dx1v(i) = dx;
//    dx1f(i) = dx;
  }

  dx = (block_size.x2max - block_size.x2min)/(Real)block_size.nx2;
  for (int j=0; j<ncells2; ++j) {
    dx2v(j) = dx;
    dx2f(j) = dx;
  }

  dx = (block_size.x3max - block_size.x3min)/(Real)block_size.nx3;
  for (int k=0; k<ncells3; ++k) {
    dx3v(k) = dx;
    dx3f(k) = dx;
  }

// Grid positions are calculated starting from both xmin and xmax and averaged
// to reduce round-off error
// assumes uniform Cartesian mesh for now

  x1f(is) = block_size.x1min;
  x2f(js) = block_size.x2min;
  x3f(ks) = block_size.x3min;

  x1f(ie+1) = block_size.x1max;
  x2f(je+1) = block_size.x2max;
  x3f(ke+1) = block_size.x3max;

// Construct Fluid and Geometry objects.  Fluid is initialized in seperate functions
// depending on whether this block is created at the start of a new simulation, or as
// part of a new domain in an AMR simulation.  Positions and spacing of cell volume
// centers depend on geometry and so are initialized in the Geometry class.

  pfluid = new Fluid(this);
  pgeometry = new Geometry(this);

  return;
}

// Block destructor

Block::~Block()
{
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();

  dx1v.DeleteAthenaArray();  
  dx2v.DeleteAthenaArray();  
  dx3v.DeleteAthenaArray();  
  dx1f.DeleteAthenaArray();  
  dx2f.DeleteAthenaArray();  
  dx3f.DeleteAthenaArray();  

  delete pfluid;
  delete pgeometry;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Mesh::InitializeAcrossDomains(enum QuantityToBeInit qnty, ParameterInput *pin)
{

// Eventually this will be a loop over all domains

  if (pdomain->pblock != NULL)  {

    switch (qnty) {
      case fluid:
        Fluid *pf = pdomain->pblock->pfluid;

// call problem generator to set initial conditions
        pf->InitProblem(pin);

// set function pointers for BCs, then apply BCs for u
        pf->pf_bcs->InitBoundaryConditions(pin);
        pf->pf_bcs->ApplyBoundaryConditions(pf->u);

// compute w everywhere (including ghost zones)
        pf->pcons_to_prim->ComputePrimitives(pf->u, pf->w);

        break;
    }

  }
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Mesh::UpdateAcrossDomains(enum UpdateAction action)
{
// Eventually this will be a loop over all domains

  if (pdomain->pblock != NULL)  {
    Fluid *pf = pdomain->pblock->pfluid;

    switch (action) {
      case fluid_bvals_n:
        pf->pf_bcs->ApplyBoundaryConditions(pf->u);
        break;
      case fluid_bvals_nhalf:
        pf->pf_bcs->ApplyBoundaryConditions(pf->u1);
        break;
      case fluid_predict:
        pf->pf_integrator->Predict(pdomain->pblock);
        break;
      case fluid_correct:
        pf->pf_integrator->Correct(pdomain->pblock);
        break;
      case convert_vars_n:
        pf->pcons_to_prim->ComputePrimitives(pf->u,pf->w);
        break;
      case convert_vars_nhalf:
        pf->pcons_to_prim->ComputePrimitives(pf->u1,pf->w1);
        break;
      case new_timestep:
        pf->NewTimeStep(pdomain->pblock);
        break;
      case data_output:
        break;
    }

  }
}
