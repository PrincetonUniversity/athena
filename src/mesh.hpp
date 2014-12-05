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
#include "blockuid/blockuid.hpp"
#include "outputs/wrapper.hpp"

class ParameterInput;
class Mesh;
class Coordinates;
class Fluid;
class Field;
class BoundaryValues;

//! \struct RegionSize
//  \brief physical size and number of cells in a Mesh

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


//! \struct NeighborBlock
//  \brief neighbor rank, level, and ids

typedef struct NeighborBlock {
  int rank, level, gid;
  NeighborBlock() : rank(-1), level(-1), gid(-1) {};
} NeighborBlock;

//! \class MeshBlock
//  \brief data/functions associated with a single block

class MeshBlock {
private:
  BlockUID uid;
  NeighborBlock neighbor[6][2][2];
  Real cost;
  friend class RestartOutput;
public:
  MeshBlock(int igid, BlockUID iuid, RegionSize input_size,
            RegionBCs input_bcs, Mesh *pm, ParameterInput *pin);
  MeshBlock(int igid, Mesh *pm, ParameterInput *pin, BlockUID *list,
                     int *nslist, ResFile& resfile, ResSize_t offset, Real icost);
  ~MeshBlock();
  size_t GetBlockSize(void);

  RegionSize block_size;
  RegionBCs  block_bcs;
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock
  void SetNeighbor(enum direction, int nrank, int nid, int nlevel);
  void SetNeighbor(enum direction, int nrank, int nid, int nlevel, int fb1, int fb2);

  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  int is,ie,js,je,ks,ke;
  int gid;

  Coordinates *pcoord;
  Fluid *pfluid;
  Field *pfield;
  BoundaryValues *pbval;
  MeshBlock *prev, *next;
};


//! \class DefaultMeshFunctions
//  \brief default functions which can be overridden by user-defined ones
class DefaultMeshFunctions {
public:
  Real MeshGeneratorX1() {};
  Real MeshGeneratorX2() {};
  Real MeshGeneratorX3() {};
  Real MeshGeneratorX1(Real rat1) {};
  Real MeshGeneratorX2(Real rat2) {};
  Real MeshGeneratorX3(Real rat3) {};
  void MeshRestart() {};
};

//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh : public DefaultMeshFunctions {
private:
  int root_level, max_level;
  int nbtotal, nblocal, nbstart, nbend;
  friend class RestartOutput;
public:
  Mesh(ParameterInput *pin, int test_flag=0);
  Mesh(ParameterInput *pin, const char *prestart_file, int test_flag=0);
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

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nid, int nlevel)
// \brief set neighbor information, for the same or a coarser level
inline void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nid, int nlevel)
{
  neighbor[dir][0][0].rank=nrank;
  neighbor[dir][0][0].gid=nid;
  neighbor[dir][0][0].level=nlevel;
}

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nid, int nlevel,
//                                  int fb1, int fb2)
// \brief set neighbor information, for a finer level
inline void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nid, int nlevel,
                                  int fb1, int fb2)
{
  neighbor[dir][fb2][fb1].rank=nrank;
  neighbor[dir][fb2][fb1].gid=nid;
  neighbor[dir][fb2][fb1].level=nlevel;
}

#endif
