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

// C/C++ headers
#include <stdint.h>  // int64_t

// Athena++ classes headers
#include "athena.hpp"
#include "athena_arrays.hpp"
#include "meshblocktree.hpp"
#include "outputs/wrapper.hpp"
#include "task_list.hpp"
#include "mesh_refinement.hpp"
#include "bvals/bvals.hpp"

// Forward declarations
class ParameterInput;
class Mesh;
class MeshRefinement;
class Coordinates;
class Hydro;
class Field;
class BoundaryValues;
class TaskList;
class MeshBlockTree;

//! \struct NeighborBlock
//  \brief neighbor rank, level, and ids

typedef struct NeighborBlock {
  int rank, level, gid, lid, ox1, ox2, ox3, fi1, fi2, bufid, targetid;
  enum neighbor_type type;
  enum BoundaryFace fid;
  enum edgeid eid;
  bool polar;
  bool self_neighbor;
  NeighborBlock() : rank(-1), level(-1), gid(-1), lid(-1), ox1(-1), ox2(-1), ox3(-1),
    bufid(-1), targetid(-1), fi1(-1), fi2(-1), type(neighbor_none),
    fid(FACE_UNDEF), eid (edgeid_undefined) {};
  void SetNeighbor(int irank, int ilevel, int igid, int ilid, int iox1, int iox2,
                   int iox3, enum neighbor_type itype, int ibid, int itargetid,
                   bool iself_neighbor, int ifi1, int ifi2, bool ipolar);
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
  LogicalLocation loc;
  NeighborBlock neighbor[56];
  Real cost;
  Real new_block_dt;
  unsigned long int finished_tasks[4];
  int first_task, num_tasks_todo, nneighbor;

  friend class RestartOutput;
  friend class BoundaryValues;
  friend class Mesh;
  friend class Hydro;
  friend class Coordinates;
  friend class TaskList;
  friend class MeshRefinement;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif
public:
  MeshBlock(int igid, int ilid, LogicalLocation iloc, RegionSize input_size,
            enum BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin);
  MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin, LogicalLocation *llist,
  IOWrapper& resfile, IOWrapperSize_t offset, Real icost, int *ranklist, int *nslist);
  ~MeshBlock();
  size_t GetBlockSizeInBytes(void);
  void SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist);
  int FindNeighborGID(int ox1, int ox2, int ox3);
  void IntegrateConservative(Real *tcons);

  RegionSize block_size;
  enum BoundaryFlag block_bcs[6];
  int nblevel[3][3][3];
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock

  int is,ie,js,je,ks,ke;
  int gid, lid;

  int cis,cie,cjs,cje,cks,cke,cnghost;

  Coordinates *pcoord, *pcoarsec;
  Hydro *phydro;
  Field *pfield;
  BoundaryValues *pbval;
  MeshRefinement *pmr;

  MeshBlock *prev, *next;
};

//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh {
private:
  int root_level, max_level, current_level;
  int nbtotal, nbstart, nbend;
  int maxneighbor_;
  int num_mesh_threads_;
  int *nslist, *nblist, *ranklist;
  Real *costlist;
  Real MeshGeneratorX1(Real x, RegionSize rs);
  Real MeshGeneratorX2(Real x, RegionSize rs);
  Real MeshGeneratorX3(Real x, RegionSize rs);
  LogicalLocation *loclist;
  MeshBlockTree tree;
  long int nrbx1, nrbx2, nrbx3;

  void MeshTest(int dim);

  friend class RestartOutput;
  friend class MeshBlock;
  friend class BoundaryValues;
  friend class Coordinates;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif
public:
  Mesh(ParameterInput *pin, int test_flag=0);
  Mesh(ParameterInput *pin, IOWrapper &resfile, int test_flag=0);
  ~Mesh();

  RegionSize mesh_size;
  enum BoundaryFlag mesh_bcs[6];

  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle;
  bool adaptive, multilevel, face_only;

  TaskList *ptlist;
  MeshBlock *pblock;

  int64_t GetTotalCells(void);
  int GetNumMeshThreads() const {return num_mesh_threads_;}
  void Initialize(int res_flag, ParameterInput *pin);
  void UpdateOneStep(void);
  void ProblemGenerator(Hydro *phyd, Field *pfld, ParameterInput *pin); // in /pgen
  void NewTimeStep(void);
  MeshBlock* FindMeshBlock(int tgid);
  void TestConservation(void);
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
