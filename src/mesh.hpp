#ifndef MESH_HPP
#define MESH_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file mesh.hpp
//  \brief defines classes Mesh, and MeshBlock
//  These classes contain data and functions related to the computational mesh
//======================================================================================

// Athena headers
#include "athena.hpp"         // macros, Real
#include "athena_arrays.hpp"  // AthenaArray

class ParameterInput;
class Mesh;
class Coordinates;
class Fluid;
class Field;
class BoundaryValues;

//! \struct RegionSize
//  \brief physical size and number of cells in a Mesh or MeshBlock

typedef struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real x1rat, x2rat, x3rat; // ratio of x(i)/x(i-1)
  int nx1, nx2, nx3;        // number of active cells (not including ghost zones)
} RegionSize;

//! \struct RegionBCs
//  \brief boundary condition flags for a Mesh or MeshBlock

typedef struct RegionBCs {
  int ix1_bc, ix2_bc, ix3_bc;  // inner-x (left edge) BC flags
  int ox1_bc, ox2_bc, ox3_bc;  // outer-x (right edge) BC flags
} RegionBCs;

//! \class MeshBlock
//  \brief data/functions associated with a single block

class MeshBlock {
public:
  MeshBlock(RegionSize input_size, RegionBCs input_bcs, Mesh *pm, ParameterInput *pin);
  ~MeshBlock();
  RegionSize block_size;
  RegionBCs  block_bcs;
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock

  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  int is,ie,js,je,ks,ke;

  Coordinates *pcoord;
  Fluid *pfluid;
  Field *pfield;
  BoundaryValues *pbval;

};

//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh {
public:
  Mesh(ParameterInput *pin);
  ~Mesh();
  RegionSize mesh_size;
  RegionBCs  mesh_bcs;

  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle;
  int nthreads_mesh;

  MeshBlock *pblock;

  void ForAllMeshBlocks(enum ActionOnBlock action, ParameterInput *pin);
  void ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin); // files in /pgen
};
#endif
