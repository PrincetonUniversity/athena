#ifndef MESH_HPP
#define MESH_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file mesh.hpp
 *  \brief defines class Mesh
 *  Contains data structures and functions related to the computational mesh
 *====================================================================================*/

class ParameterInput;

//! \struct RegionSize
//  \brief physical size and number of cells in a region (Domain or Block)

typedef struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  int nx1, nx2, nx3;
} RegionSize;

//! \struct Block
//  \brief data associated with a single block inside a domain

struct FluidData;

typedef struct Block {
  AthenaArray<Real> x1v, x2v, x3v, dx1v, dx2v, dx3v;
  AthenaArray<Real> x1f, x2f, x3f, dx1f, dx2f, dx3f; 
  int is,ie,js,je,ks,ke;
  RegionSize block_size;
} Block;

//! \struct Domain
//  \brief data associated with a domain

typedef struct Domain {
  Block *pblock;
  RegionSize domain_size;
} Domain;

//! \class Mesh
//  \brief mesh data and functions

class Mesh {
public:
  Mesh(ParameterInput *pin);
  ~Mesh();

  RegionSize mesh_size;
  Domain root;

// public functions implemented in mesh.cpp

  void InitDomain(Domain *pd);
  void InitBlocks(Block *pb);

private:

// private functions implemented in mesh.cpp

};
#endif
