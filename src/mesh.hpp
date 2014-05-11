#ifndef MESH_HPP
#define MESH_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file mesh.hpp
 *  \brief defines classes Mesh, Domain, and Block
 *  These classes contain data and functions related to the computational mesh
 *====================================================================================*/

class ParameterInput;
class Fluid;
class Mesh;
class Domain;

//! \struct RegionSize
//  \brief physical size and number of cells in a region (Domain or Block)

typedef struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  int nx1, nx2, nx3;
} RegionSize;

//! \class Block
//  \brief data/functions associated with a single block inside a domain

class Block {
public:
  Block(RegionSize region, Domain *pd);
  ~Block();

  Fluid *pfluid;

  AthenaArray<Real> x1v, x2v, x3v, dx1v, dx2v, dx3v;
  AthenaArray<Real> x1f, x2f, x3f, dx1f, dx2f, dx3f; 
  int is,ie,js,je,ks,ke;
  RegionSize block_size;

  Domain *pmy_domain;
};

//! \class Domain
//  \brief data/functions associated with a domain

class Domain {
public:
  Domain(RegionSize region, Mesh *pm);
  ~Domain();

  Block *pblock;
  RegionSize domain_size;
  Mesh *pmy_mesh;
};

//! \class Mesh
//  \brief data/functions associated with the mesh

class Mesh {
public:
  Mesh(ParameterInput *pin);
  ~Mesh();

  Domain *pdomain;
  RegionSize mesh_size;

  Real time, start_time, dt, tlim, cfl_number;
  int ncycle, nlim;

  void InitializeOnDomains(enum QuantityToBeInitialized qnty, ParameterInput *pin);
  void StepThroughDomains(enum AlgorithmSteps action);

};
#endif
