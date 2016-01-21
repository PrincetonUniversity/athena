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
#include "bvals/bvals.hpp"
#include "mesh_refinement/mesh_refinement.hpp"

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
  bool polar; // flag indicating boundary is across a pole
  NeighborBlock() : rank(-1), level(-1), gid(-1), lid(-1), ox1(-1), ox2(-1), ox3(-1),
    bufid(-1), targetid(-1), fi1(-1), fi2(-1), type(neighbor_none),
    fid(FACE_UNDEF), eid (edgeid_undefined), polar(false) {};
  void SetNeighbor(int irank, int ilevel, int igid, int ilid, int iox1, int iox2,
                   int iox3, enum neighbor_type itype, int ibid, int itargetid,
                   bool ipolar, int ifi1, int ifi2);
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
  NeighborBlock neighbor[56];
  Real cost;
  Real new_block_dt;
  unsigned long int finished_tasks[4];
  int first_task, num_tasks_todo, nneighbor;

  void ProblemGenerator(ParameterInput *pin); // in /pgen

  friend class RestartOutput;
  friend class BoundaryValues;
  friend class Mesh;
  friend class Hydro;
  friend class TaskList;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif
public:
  LogicalLocation loc;
  MeshBlock(int igid, int ilid, LogicalLocation iloc, RegionSize input_size,
            enum BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin);
  MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin, LogicalLocation iloc,
  IOWrapper& resfile, IOWrapperSize_t offset, Real icost, int *ranklist, int *nslist);
  ~MeshBlock();
  size_t GetBlockSizeInBytes(void);
  void SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist);
  int FindNeighborGID(int ox1, int ox2, int ox3);
  void IntegrateConservative(Real *tcons);
  void UserWorkInLoop(void); // in /pgen

  RegionSize block_size;
  enum BoundaryFlag block_bcs[6];
  int nblevel[3][3][3];
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock

  int is,ie,js,je,ks,ke;
  int gid, lid;

  int cis,cie,cjs,cje,cks,cke,cnghost;

  Coordinates *pcoord;
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
  int nbtotal;
  int maxneighbor_;
  int num_mesh_threads_;
  int *nslist, *ranklist, *nblist;
  Real *costlist;
  int *nref, *nderef, *bnref, *bnderef, *rdisp, *brdisp, *ddisp, *bddisp;
  LogicalLocation *loclist;
  MeshBlockTree tree;
  long int nrbx1, nrbx2, nrbx3;

  bool user_meshgen_[3];
  MeshGenFunc_t MeshGenerator_[3];
  SrcTermFunc_t UserSourceTerm_;
  BValFunc_t BoundaryFunction_[6];
  AMRFlag_t AMRFlag_;

  void MeshTest(int dim);

  // methods in /pgen
  void InitUserMeshProperties(ParameterInput *pin);
  void TerminateUserMeshProperties(void);

  void EnrollUserBoundaryFunction (enum BoundaryFace face, BValFunc_t my_func);
  void EnrollUserRefinementCondition(AMRFlag_t amrflag);
  void EnrollUserMeshGenerator(enum direction dir, MeshGenFunc_t my_mg);
  void EnrollUserSourceTermFunction(SrcTermFunc_t my_func);

  void LoadBalancing(Real *clist, int *rlist, int *slist, int *nlist, int nb);

  friend class RestartOutput;
  friend class MeshBlock;
  friend class BoundaryValues;
  friend class Coordinates;
  friend class MeshRefinement;
  friend class HydroSourceTerms;
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
  void SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                 enum BoundaryFlag *block_bcs);
  void UpdateOneStep(void);
  void NewTimeStep(void);
  void AdaptiveMeshRefinement(ParameterInput *pin);
  MeshBlock* FindMeshBlock(int tgid);
  void TestConservation(void);
};


//--------------------------------------------------------------------------------------
// \!fn Real DefaultMeshGeneratorX1(Real x, RegionSize rs)
// \brief x1 mesh generator function, x is the logical location; x=i/nx1
inline Real DefaultMeshGeneratorX1(Real x, RegionSize rs)
{
  Real lw, rw;
  if(rs.x1rat==1.0)
    rw=x, lw=1.0-x;
  else
  {
    Real ratn=pow(rs.x1rat,rs.nx1);
    Real rnx=pow(rs.x1rat,x*rs.nx1);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  return rs.x1min*lw+rs.x1max*rw;
}

//--------------------------------------------------------------------------------------
// \!fn Real DefaultMeshGeneratorX2(Real x, RegionSize rs)
// \brief x2 mesh generator function, x is the logical location; x=j/nx2
inline Real DefaultMeshGeneratorX2(Real x, RegionSize rs)
{
  Real lw, rw;
  if(rs.x2rat==1.0)
    rw=x, lw=1.0-x;
  else
  {
    Real ratn=pow(rs.x2rat,rs.nx2);
    Real rnx=pow(rs.x2rat,x*rs.nx2);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  return rs.x2min*lw+rs.x2max*rw;
}

//--------------------------------------------------------------------------------------
// \!fn Real DefaultMeshGeneratorX3(Real x, RegionSize rs)
// \brief x3 mesh generator function, x is the logical location; x=k/nx3
inline Real DefaultMeshGeneratorX3(Real x, RegionSize rs)
{
  Real lw, rw;
  if(rs.x3rat==1.0)
    rw=x, lw=1.0-x;
  else
  {
    Real ratn=pow(rs.x3rat,rs.nx3);
    Real rnx=pow(rs.x3rat,x*rs.nx3);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  return rs.x3min*lw+rs.x3max*rw;
}

#endif
