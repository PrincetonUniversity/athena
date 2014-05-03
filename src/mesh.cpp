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

//======================================================================================
/*! \file mesh.cpp
 *  \brief implementation of functions in class Mesh
 *====================================================================================*/

// constructor, builds mesh based on parameters in input file

Mesh::Mesh(ParameterInput *pin)
{
// read size of root domain from input file.  
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

// set region size of root domain, and initialize

  root.domain_size = mesh_size;
  InitDomain(&root);
}

// destructor

Mesh::~Mesh()
{
  delete[] root.pblock;
}

//--------------------------------------------------------------------------------------
/*! \fn  void Mesh::InitDomain()
 *  \brief Initializes Domain struct */  

void Mesh::InitDomain(Domain *pd)
{

// allocate blocks in domain, set their region sizes, and initialize

  pd->pblock = new Block;
  pd->pblock->block_size = pd->domain_size;
  InitBlocks(pd->pblock); 

  return;
}

//--------------------------------------------------------------------------------------
/*! \fn  void Mesh::InitBlocks()
 *  \brief Initializes data in Blocks struct */  

void Mesh::InitBlocks(Block *pb)
{

// initilize grid indices
// assumes 3D for now, need to add error checks

  int ncells1 = pb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = pb->block_size.nx2 + 2*(NGHOST);
  int ncells3 = pb->block_size.nx3 + 2*(NGHOST);

  pb->is = NGHOST;
  pb->ie = pb->is + pb->block_size.nx1 - 1;

  pb->js = NGHOST;
  pb->je = pb->js + pb->block_size.nx2 - 1;

  pb->ks = NGHOST;
  pb->ke = pb->ks + pb->block_size.nx3 - 1;

  std::cout << "is=" << pb->is << " ie=" << pb->ie << std::endl;
  std::cout << "js=" << pb->js << " je=" << pb->je << std::endl;
  std::cout << "ks=" << pb->ks << " ke=" << pb->ke << std::endl;

// allocate arrays for grid spacing and positions
// assumes 3D for now

  pb->x1v.NewAthenaArray(ncells1);
  pb->x2v.NewAthenaArray(ncells2);
  pb->x3v.NewAthenaArray(ncells3);

  pb->dx1v.NewAthenaArray(ncells1);
  pb->dx2v.NewAthenaArray(ncells2);
  pb->dx3v.NewAthenaArray(ncells3);

  pb->x1f.NewAthenaArray(ncells1 + 1);
  pb->x2f.NewAthenaArray(ncells2 + 1);
  pb->x3f.NewAthenaArray(ncells3 + 1);

  pb->dx1f.NewAthenaArray(ncells1);
  pb->dx2f.NewAthenaArray(ncells2);
  pb->dx3f.NewAthenaArray(ncells3);

// initialize grid spacing
// assumes uniform Cartesian mesh for now

  Real dx = (pb->block_size.x1max - pb->block_size.x1min)/(Real)pb->block_size.nx1;
  for (int i=0; i<ncells1; ++i) {
    pb->dx1v(i) = dx;
    pb->dx1f(i) = dx;
  }

  dx = (pb->block_size.x2max - pb->block_size.x2min)/(Real)pb->block_size.nx2;
  for (int j=0; j<ncells2; ++j) {
    pb->dx2v(j) = dx;
    pb->dx2f(j) = dx;
  }

  dx = (pb->block_size.x3max - pb->block_size.x3min)/(Real)pb->block_size.nx3;
  for (int k=0; k<ncells3; ++k) {
    pb->dx3v(k) = dx;
    pb->dx3f(k) = dx;
  }

// Grid positions are calculated starting from both xmin and xmax and averaged
// to reduce round-off error
// assumes uniform Cartesian mesh for now

  pb->x1f(pb->is) = pb->block_size.x1min;
  pb->x2f(pb->js) = pb->block_size.x2min;
  pb->x3f(pb->ks) = pb->block_size.x3min;

  pb->x1f(pb->ie+1) = pb->block_size.x1max;
  pb->x2f(pb->je+1) = pb->block_size.x2max;
  pb->x3f(pb->ke+1) = pb->block_size.x3max;

  return;
}
