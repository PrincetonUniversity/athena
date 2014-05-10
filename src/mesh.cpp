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
#include "bvals/boundary_conditions.hpp"
#include "convert_var/convert_var.hpp"

//======================================================================================
/*! \file mesh.cpp
 *  \brief implementation of functions in classes Mesh, Domain, and Block
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// Mesh constructor, builds mesh based on parameters in input file

Mesh::Mesh(ParameterInput *pin)
{
// read mesh (root domain) size from input file.  
// need to add error checking

  mesh_size.nx1 = pin->GetInteger("mesh","nx1");
  mesh_size.nx2 = pin->GetInteger("mesh","nx2");
  mesh_size.nx3 = pin->GetInteger("mesh","nx3");

  mesh_size.x1min = pin->GetReal("mesh","x1min");
  mesh_size.x2min = pin->GetReal("mesh","x2min");
  mesh_size.x3min = pin->GetReal("mesh","x3min");

  mesh_size.x1max = pin->GetReal("mesh","x1max");
  mesh_size.x2max = pin->GetReal("mesh","x2max");
  mesh_size.x3max = pin->GetReal("mesh","x3max");

// allocate root domain

  pdomain = new Domain(mesh_size);

}

// Mesh destructor

Mesh::~Mesh()
{
  delete[] pdomain;
}

//--------------------------------------------------------------------------------------
// Domain constructor

Domain::Domain(RegionSize region)
{
  domain_size = region;

// calculate array of blocks in domain, set their region sizes, and initialize

  pblock = new Block(domain_size);

  return;
}

// Domain destructor

Domain::~Domain()
{
  delete[] pblock;
}

//--------------------------------------------------------------------------------------
// Block constructor

Block::Block(RegionSize region)
{
  block_size = region;

// initilize grid indices
// assumes 3D for now, need to add error checks

  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = block_size.nx2 + 2*(NGHOST);
  int ncells3 = block_size.nx3 + 2*(NGHOST);

  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  js = NGHOST;
  je = js + block_size.nx2 - 1;

  ks = NGHOST;
  ke = ks + block_size.nx3 - 1;

  std::cout << "is=" << is << " ie=" << ie << std::endl;
  std::cout << "js=" << js << " je=" << je << std::endl;
  std::cout << "ks=" << ks << " ke=" << ke << std::endl;

// allocate arrays for grid spacing and positions
// assumes 3D for now

  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

  x1f.NewAthenaArray(ncells1 + 1);
  x2f.NewAthenaArray(ncells2 + 1);
  x3f.NewAthenaArray(ncells3 + 1);

  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);

// initialize grid spacing
// assumes uniform Cartesian mesh for now

  Real dx = (block_size.x1max - block_size.x1min)/(Real)block_size.nx1;
  for (int i=0; i<ncells1; ++i) {
    dx1v(i) = dx;
    dx1f(i) = dx;
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

// Create Fluid

//  pfluid = new Fluid(pin, this);

  return;
}

// Block destructor

Block::~Block()
{
// delete vectors storing cell positions and spacing

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
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Mesh::InitializeOnDomains(enum QuantityToBeInitialized qnty, ParameterInput *pin)
{

// Eventually this will be a loop over all domains

  if (pdomain->pblock != NULL)  {

    switch (qnty) {
      case fluid:
// construct new fluid, call problem generator
        pdomain->pblock->pfluid = new Fluid(pin,pdomain->pblock);
        Fluid *pf = pdomain->pblock->pfluid;
        pf->Problem(pin);

// construct new BCs, and set them for u
        pf->pbvals = new BoundaryConditions(pin,pf);
        pf->pbvals->SetBoundaryValues(pf->u);

// construction new variable conversion class, and compute w
        pf->pcons_to_prim = new ConvertVariables(pf);
        pf->pcons_to_prim->ComputePrimitives(pf->u, pf->w);

        break;
    }

  }
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Mesh::StepThroughDomains(enum AlgorithmSteps action)
{
// Eventually this will be a loop over all domains

  if (pdomain->pblock != NULL)  {
    Fluid *pf = pdomain->pblock->pfluid;

    switch (action) {
      case fluid_bvals_n:
        pf->pbvals->SetBoundaryValues(pf->u);
        break;
      case fluid_bvals_nhalf:
        pf->pbvals->SetBoundaryValues(pf->u1_);
        break;
      case fluid_predict:
        break;
      case fluid_correct:
        break;
      case convert_vars_n:
        break;
      case convert_vars_nhalf:
        break;
      case data_output:
        break;
    }

  }
}

