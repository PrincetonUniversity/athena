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


#include <stdint.h>  // int64_t

// Athena headers
#include "athena.hpp"         // macros, Real
#include "athena_arrays.hpp"  // AthenaArray
#include "blockuid.hpp"
#include "wrapio.hpp"
#include "tasklist.hpp"

class ParameterInput;
class Mesh;
class Coordinates;
class Fluid;
class Field;
class BoundaryValues;
struct Task;
class TaskList;

//! \struct RegionSize
//  \brief physical size and number of cells in a Mesh

typedef struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real x1rat, x2rat, x3rat; // ratio of x(i)/x(i-1)
  int nx1, nx2, nx3;        // number of active cells (not including ghost zones)
} RegionSize;


//! \class MeshBlock
//  \brief data/functions associated with a single block

class MeshBlock {
private:
  BlockUID uid;
  NeighborBlock neighbor[6][2][2];
  Real cost;
  Real new_block_dt;
  Task *task;
  long int task_flag;
  int ntask, firsttask, ntodo;

  friend class RestartOutput;
  friend class BoundaryValues;
  friend class Mesh;
  friend class Fluid;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif
public:
  MeshBlock(int igid, int ilid, BlockUID iuid, RegionSize input_size,
            int *input_bcs, Mesh *pm, ParameterInput *pin);
  MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin, BlockUID *list,
  WrapIO& resfile, WrapIOSize_t offset, Real icost, int *ranklist, int *nslist);
  ~MeshBlock();
  size_t GetBlockSizeInBytes(void);

  void SetNeighbor(enum direction, int nrank, int nlevel, int ngid, int nlid);
  void SetNeighbor(enum direction, int nrank, int nlevel, int ngid, int nlid,
                   int fb1, int fb2);

  void SetTaskList(TaskList& tl);
  enum tlstatus DoOneTask(void);

  RegionSize block_size;
  int block_bcs[6];
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock

  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  int is,ie,js,je,ks,ke;
  int gid, lid;

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
  int num_mesh_threads_;
  int *nslist, *nblist;
  Real MeshGeneratorX1(Real x, RegionSize rs);
  Real MeshGeneratorX2(Real x, RegionSize rs);
  Real MeshGeneratorX3(Real x, RegionSize rs);
  bool adaptive;

  void MeshTest(BlockUID *buid, int *ranklist, Real *costlist);

  friend class RestartOutput;
  friend class MeshBlock;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif
public:
  Mesh(ParameterInput *pin, int test_flag=0);
  Mesh(ParameterInput *pin, WrapIO &resfile, int test_flag=0);
  ~Mesh();

  RegionSize mesh_size;
  int mesh_bcs[6];

  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle;

  MeshBlock *pblock;

  int64_t GetTotalCells(void);
  int GetNumMeshThreads() const {return num_mesh_threads_;}
  void Initialize(int res_flag, ParameterInput *pin);
  void UpdateOneStep(void);
  void SetTaskList(TaskList& tl);
  void ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin); // files in /pgen
  void NewTimeStep(void);
};

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel,
//                                  int ngid, int nlid)
// \brief set neighbor information, for the same or a coarser level
inline void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel,
                                   int ngid, int nlid)
{
  neighbor[dir][0][0].rank=nrank;
  neighbor[dir][0][0].level=nlevel;
  neighbor[dir][0][0].gid=ngid;
  neighbor[dir][0][0].lid=nlid;
}

//--------------------------------------------------------------------------------------
// \!fn void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel,
//                                  int ngid, int nlid, int fb1, int fb2)
// \brief set neighbor information, for a finer level
inline void MeshBlock::SetNeighbor(enum direction dir, int nrank, int nlevel,
                                  int ngid, int nlid, int fb1, int fb2)
{
  neighbor[dir][fb2][fb1].rank=nrank;
  neighbor[dir][fb2][fb1].level=nlevel;
  neighbor[dir][fb2][fb1].gid=ngid;
  neighbor[dir][fb2][fb1].lid=nlid;
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
