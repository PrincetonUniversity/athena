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
class Geometry;

//! \struct RegionSize
//  \brief physical size and number of cells in a region (Domain or Block)

typedef struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real x1rat, x2rat, x3rat; // ratio of x(i)/x(i-1)
  int nx1, nx2, nx3;        // number of active cells (not including ghost zones)
} RegionSize;

//! \class Block
//  \brief data/functions associated with a single block inside a domain

class Block {
public:
  Block(RegionSize blk_size, Domain *pd);
  ~Block();

  Domain *pparent_domain;                            // ptr to parent Domain
  AthenaArray<Real> x1f, x2f, x3f, dx1f, dx2f, dx3f; // cell face   centers and spacing
  AthenaArray<Real> x1v, x2v, x3v, dx1v, dx2v, dx3v; // cell volume centers and spacing
  int is,ie,js,je,ks,ke;
  RegionSize block_size;

  Fluid *pfluid;
  Geometry *pgeometry;
};

//! \class Domain
//  \brief data/functions associated with a domain inside the mesh

class Domain {
public:
  Domain(RegionSize dom_size, Mesh *pm);
  ~Domain();

  Mesh *pparent_mesh;  // ptr to parent Mesh

  Block *pblock;
  RegionSize domain_size;
};

//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh {
public:
  Mesh(ParameterInput *pin);
  ~Mesh();

  Domain *pdomain;
  RegionSize mesh_size;

  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle;

  void InitializeAcrossDomains(enum QuantityToBeInit qnty, ParameterInput *pin);
  void UpdateAcrossDomains(enum UpdateAction action);
};
#endif
