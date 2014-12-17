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
#include "blockuid.hpp"
#include "wrapio.hpp"

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


//! \class MeshBlock
//  \brief data/functions associated with a single block

class MeshBlock {
private:
  BlockUID uid;
  NeighborBlock neighbor[6][2][2];
  Real cost;
  friend class RestartOutput;
  friend class BoundaryValues;
public:
  MeshBlock(int igid, BlockUID iuid, RegionSize input_size,
            RegionBCs input_bcs, Mesh *pm, ParameterInput *pin);
  MeshBlock(int igid, Mesh *pm, ParameterInput *pin, BlockUID *list,
            WrapIO& resfile, WrapIOSize_t offset, Real icost, int *ranklist);
  ~MeshBlock();
  size_t GetBlockSizeInBytes(void);

  void SetNeighbor(enum direction, int nrank, int nlevel, int nid);
  void SetNeighbor(enum direction, int nrank, int nlevel, int nid, int fb1, int fb2);

  RegionSize block_size;
  RegionBCs  block_bcs;
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock

  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  int is,ie,js,je,ks,ke;
  int gid;
  Real block_dt;

  Coordinates *pcoord;
  Fluid *pfluid;
  Field *pfield;
  BoundaryValues *pbval;
  MeshBlock *prev, *next;
};

//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh {
private:
  int root_level, max_level;
  int nbtotal, nbstart, nbend;
  Real MeshGeneratorX1(Real x, RegionSize rs);
  Real MeshGeneratorX2(Real x, RegionSize rs);
  Real MeshGeneratorX3(Real x, RegionSize rs);
  bool adaptive;

  friend class RestartOutput;
  friend class MeshBlock;
public:
  Mesh(ParameterInput *pin, int test_flag=0);
  Mesh(ParameterInput *pin, WrapIO &resfile, int test_flag=0);
  ~Mesh();

  RegionSize mesh_size;
  RegionBCs  mesh_bcs;

  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle;
  int nthreads_mesh;

  MeshBlock *pblock;

  void ForAllMeshBlocks(enum ActionOnBlock action, ParameterInput *pin);
  void ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin); // files in /pgen
  void NewTimeStep(void);
};

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel, int nid)
// \brief set neighbor information, for the same or a coarser level
inline void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel, int nid)
{
  neighbor[dir][0][0].rank=nrank;
  neighbor[dir][0][0].level=nlevel;
  neighbor[dir][0][0].gid=nid;
}

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nid, int nlevel,
//                                  int fb1, int fb2)
// \brief set neighbor information, for a finer level
inline void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel, int nid, 
                                  int fb1, int fb2)
{
  neighbor[dir][fb2][fb1].rank=nrank;
  neighbor[dir][fb2][fb1].level=nlevel;
  neighbor[dir][fb2][fb1].gid=nid;
}


//--------------------------------------------------------------------------------------
// \!fn Real Mesh::MeshGeneratorX1(Real x, RegionSize rs)
// \brief x1 mesh generator function, x is the logical location; x=i/nx1
inline Real Mesh::MeshGeneratorX1(Real x, RegionSize rs)
{
  Real lw, rw;
  if(rs.x1rat==1.0)
    rw=x;
  else
  {
    Real ratn=pow(rs.x1rat,rs.nx1);
    Real rnx=pow(rs.x1rat,x*rs.nx1);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  lw=1.0-rw;
  return rs.x1min*lw+rs.x1max*rw;
}

//--------------------------------------------------------------------------------------
// \!fn Real Mesh::MeshGeneratorX2(Real x, RegionSize rs)
// \brief x2 mesh generator function, x is the logical location; x=j/nx2
inline Real Mesh::MeshGeneratorX2(Real x, RegionSize rs)
{
  Real lw, rw;
  if(rs.x2rat==1.0)
    rw=x;
  else
  {
    Real ratn=pow(rs.x2rat,rs.nx2);
    Real rnx=pow(rs.x2rat,x*rs.nx2);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  lw=1.0-rw;
  return rs.x2min*lw+rs.x2max*rw;
}

//--------------------------------------------------------------------------------------
// \!fn Real Mesh::MeshGeneratorX3(Real x, RegionSize rs)
// \brief x3 mesh generator function, x is the logical location; x=k/nx3
inline Real Mesh::MeshGeneratorX3(Real x, RegionSize rs)
{
  Real lw, rw;
  if(rs.x3rat==1.0)
    rw=x;
  else
  {
    Real ratn=pow(rs.x3rat,rs.nx3);
    Real rnx=pow(rs.x3rat,x*rs.nx3);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  lw=1.0-rw;
  return rs.x3min*lw+rs.x3max*rw;
}

#endif
