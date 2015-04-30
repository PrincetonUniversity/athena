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
class BlockUID;
class BlockTree;


//! \struct NeighborBlock
//  \brief neighbor rank, level, and ids

typedef struct NeighborBlock {
  int rank, level, gid, lid, ox1, ox2, ox3, bufid, fi1, fi2;
  enum neighbor_type type;
  NeighborBlock() : rank(-1), level(-1), gid(-1), lid(-1),
    ox1(-1), ox2(-1), ox3(-1), bufid(-1), fi1(-1), fi2(-1), type(neighbor_none) {};
  void SetNeighbor(int irank, int ilevel, int igid, int ilid, int iox1, int iox2, int iox3,
                   enum neighbor_type itype, int ibid, int ifi1, int ifi2);
} NeighborBlock;



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
  NeighborBlock neighbor[56];
  Real cost;
  Real new_block_dt;
  Task *task;
  long int task_flag;
  int ntask, firsttask, ntodo, nneighbor;

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
  void SearchAndSetNeighbors(BlockTree &tree, int *ranklist, int *nslist);
  void SetTaskList(TaskList& tl);
  enum tasklist_status DoOneTask(void);
  void SetCoarserCoordinates(void);

  RegionSize block_size;
  int block_bcs[6];
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock

  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  int is,ie,js,je,ks,ke;
  int gid, lid;

  AthenaArray<Real> coarse_dx1f, coarse_dx2f, coarse_dx3f;
  AthenaArray<Real> coarse_x1f, coarse_x2f, coarse_x3f;
  AthenaArray<Real> coarse_dx1v, coarse_dx2v, coarse_dx3v;
  AthenaArray<Real> coarse_x1v, coarse_x2v, coarse_x3v;
  AthenaArray<Real> coarse_data;
  int cis,cie,cjs,cje,cks,cke,cnghost;

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
  int maxneighbor_;
  int num_mesh_threads_;
  int *nslist, *nblist;
  Real MeshGeneratorX1(Real x, RegionSize rs);
  Real MeshGeneratorX2(Real x, RegionSize rs);
  Real MeshGeneratorX3(Real x, RegionSize rs);
  bool adaptive, multilevel, face_only;
  BlockUID *buid;
  BlockTree tree;
  long int nrbx1, nrbx2, nrbx3;

  void MeshTest(int *ranklist, Real *costlist, int dim);

  friend class RestartOutput;
  friend class MeshBlock;
  friend class BoundaryValues;
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
  MeshBlock* FindMeshBlock(int tgid);
};


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
