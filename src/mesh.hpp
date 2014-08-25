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

// Athena headers
#include "athena.hpp"         // macros, Real
#include "athena_arrays.hpp"  // AthenaArray

class ParameterInput;
class Mesh;
class Domain;
class Coordinates;
class FluidBoundaryConditions;
class Fluid;
class OutputList;

//! \struct RegionSize
//  \brief physical size and number of cells in a Mesh, Domain or Block

typedef struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real x1rat, x2rat, x3rat; // ratio of x(i)/x(i-1)
  int nx1, nx2, nx3;        // number of active cells (not including ghost zones)
} RegionSize;

//! \struct RegionBoundaryFlags
//  \brief boundary condition flags for a Mesh, Domain or Block

typedef struct RegionBoundaryFlags {
  int ix1_bc, ix2_bc, ix3_bc;  // inner-x (left edge) BC flags
  int ox1_bc, ox2_bc, ox3_bc;  // outer-x (right edge) BC flags
} RegionBoundaryFlags;

//--------------------------------------------------------------------------------------
//! \class Block
//  \brief data/functions associated with a single block inside a domain

class Block {
public:
  Block(RegionSize input_size, RegionBoundaryFlags input_bndry, Domain *pd);
  ~Block();

  Domain *pmy_domain;  // ptr to Domain containing this Block
  RegionSize block_size;
  RegionBoundaryFlags block_bndry;

  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  int is,ie,js,je,ks,ke;

  FluidBoundaryConditions *pf_bcs;
  Coordinates *pcoord;
  Fluid *pfluid;
  OutputList *poutputs;
};

//--------------------------------------------------------------------------------------
//! \class Domain
//  \brief data/functions associated with a domain inside the mesh

class Domain {
public:
  Domain(RegionSize input_size, RegionBoundaryFlags input_bndry, Mesh *pm);
  ~Domain();

  Mesh *pmy_mesh;  // ptr to Mesh containing this Domain

  Block *pblock;
  RegionSize domain_size;
  RegionBoundaryFlags domain_bndry;
};

//--------------------------------------------------------------------------------------
//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh {
public:
  Mesh(ParameterInput *pin);
  ~Mesh();

  Domain *pdomain;
  RegionSize mesh_size;
  RegionBoundaryFlags mesh_bndry;

  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle;

  void ForAllDomains(enum ActionOnDomain action, ParameterInput *pin);
};
#endif
