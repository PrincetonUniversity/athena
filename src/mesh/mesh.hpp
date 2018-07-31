#ifndef MESH_MESH_HPP_
#define MESH_MESH_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.hpp
//  \brief defines Mesh and MeshBlock classes, and various structs used in them
//  The Mesh is the overall grid structure, and MeshBlocks are local patches of data
//  (potentially on different levels) that tile the entire domain.

// C/C++ headers
#include <stdint.h>  // int64_t
#include <string>

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../outputs/io_wrapper.hpp"
#include "../task_list/task_list.hpp"
#include "../bvals/bvals.hpp"
#include "meshblock_tree.hpp"
#include "mesh_refinement.hpp"

// Forward declarations
class ParameterInput;
class Mesh;
class MeshRefinement;
class MeshBlockTree;
class BoundaryValues;
class GravityBoundaryValues;
class TaskList;
class TaskState;
class Coordinates;
class Reconstruction;
class Hydro;
class Field;
class Gravity;
class MGGravityDriver;
class EquationOfState;
class FFTDriver;
class FFTGravityDriver;
class TurbulenceDriver;

//----------------------------------------------------------------------------------------
//! \class MeshBlock
//  \brief data/functions associated with a single block

class MeshBlock {
  friend class RestartOutput;
  friend class BoundaryValues;
  friend class GravityBoundaryValues;
  friend class Mesh;
  friend class Hydro;
  friend class TaskList;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif

public:
  MeshBlock(int igid, int ilid, LogicalLocation iloc, RegionSize input_size,
            enum BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin, int igflag,
            bool ref_flag = false);
  MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin, LogicalLocation iloc,
            RegionSize input_block, enum BoundaryFlag *input_bcs, Real icost,
            char *mbdata, int igflag);
  ~MeshBlock();

  // data
  Mesh *pmy_mesh;  // ptr to Mesh containing this MeshBlock
  LogicalLocation loc;
  RegionSize block_size;
  int is,ie,js,je,ks,ke;
  int gid, lid;
  int cis,cie,cjs,cje,cks,cke,cnghost;
  int gflag;
  // At every cycle n, hydro and field registers (u, b) are advanced from t^n -> t^{n+1},
  // the time-integration scheme may partially substep several storage register pairs
  // (u,b), (u1,b1), (u2, b2), ..., (umn, bm) through the dt interval.
  // Track their time abscissae at the end of each stage (l) as (dt_m^l) relative to t^n
  Real stage_abscissae[MAX_NSTAGE][MAX_NREGISTER];

  // user output variables for analysis
  int nuser_out_var;
  AthenaArray<Real> user_out_var;
  std::string *user_out_var_names_;

  // user MeshBlock data that can be stored in restart files
  AthenaArray<Real> *ruser_meshblock_data;
  AthenaArray<int> *iuser_meshblock_data;

  // mesh-related objects
  Coordinates *pcoord;
  BoundaryValues *pbval;
  GravityBoundaryValues *pgbval;
  Reconstruction *precon;
  MeshRefinement *pmr;

  // physics-related objects
  Hydro *phydro;
  Field *pfield;
  Gravity *pgrav;
  EquationOfState *peos;

  MeshBlock *prev, *next;

  // functions
  size_t GetBlockSizeInBytes(void);
  void SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist);
  void UserWorkInLoop(void); // in ../pgen
  void InitUserMeshBlockData(ParameterInput *pin); // in ../pgen
  void UserWorkBeforeOutput(ParameterInput *pin); // in ../pgen

private:
  // data
  Real cost;
  Real new_block_dt;
  TaskState tasks;
  int nreal_user_meshblock_data_, nint_user_meshblock_data_;

  // functions
  void AllocateRealUserMeshBlockDataField(int n);
  void AllocateIntUserMeshBlockDataField(int n);
  void AllocateUserOutputVariables(int n);
  void SetUserOutputVariableName(int n, const char *name);

  void ProblemGenerator(ParameterInput *pin); // in ../pgen
};

//----------------------------------------------------------------------------------------
//! \class Mesh
//  \brief data/functions associated with the overall mesh

class Mesh {
  friend class RestartOutput;
  friend class HistoryOutput;
  friend class MeshBlock;
  friend class BoundaryBase;
  friend class BoundaryValues;
  friend class MGBoundaryValues;
  friend class GravityBoundaryValues;
  friend class Coordinates;
  friend class MeshRefinement;
  friend class HydroSourceTerms;
  friend class Hydro;
  friend class FFTDriver;
  friend class FFTGravityDriver;
  friend class TurbulenceDriver;
  friend class MultigridDriver;
  friend class MGGravityDriver;
  friend class Gravity;
  friend class HydroDiffusion;
  friend class FieldDiffusion;
#ifdef HDF5OUTPUT
  friend class ATHDF5Output;
#endif

public:
  explicit Mesh(ParameterInput *pin, int test_flag=0);
  Mesh(ParameterInput *pin, IOWrapper &resfile, int test_flag=0);
  ~Mesh();

  // accessors
  int GetNumMeshBlocksThisRank(int my_rank) {return nblist[my_rank];}
  int GetNumMeshThreads() const {return num_mesh_threads_;}
  int64_t GetTotalCells() {return static_cast<int64_t> (nbtotal)*
     pblock->block_size.nx1*pblock->block_size.nx2*pblock->block_size.nx3;}

  // data
  RegionSize mesh_size;
  enum BoundaryFlag mesh_bcs[6];
  Real start_time, tlim, cfl_number, time, dt;
  int nlim, ncycle, ncycle_out;
  int nbtotal, nbnew, nbdel;
  bool adaptive, multilevel;
  int gflag;
  int turb_flag; // turbulence flag

  MeshBlock *pblock;

  TurbulenceDriver *ptrbd;
  FFTGravityDriver *pfgrd;
  MGGravityDriver *pmgrd;

  AthenaArray<Real> *ruser_mesh_data;
  AthenaArray<int> *iuser_mesh_data;

  // functions
  void Initialize(int res_flag, ParameterInput *pin);
  void SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                 enum BoundaryFlag *block_bcs);
  void NewTimeStep(void);
  void AdaptiveMeshRefinement(ParameterInput *pin);
  unsigned int CreateAMRMPITag(int lid, int ox1, int ox2, int ox3);
  MeshBlock* FindMeshBlock(int tgid);
  void ApplyUserWorkBeforeOutput(ParameterInput *pin);
  void UserWorkAfterLoop(ParameterInput *pin); // method in ../pgen

private:
  // data
  int root_level, max_level, current_level;
  int num_mesh_threads_;
  int *nslist, *ranklist, *nblist;
  Real *costlist;
  int *nref, *nderef, *bnref, *bnderef, *rdisp, *brdisp, *ddisp, *bddisp;
  LogicalLocation *loclist;
  MeshBlockTree tree;
  int64_t nrbx1, nrbx2, nrbx3;
  // flags are false if using non-uniform or user meshgen function
  bool use_uniform_meshgen_fn_[3];
  int nreal_user_mesh_data_, nint_user_mesh_data_;

  int nuser_history_output_;
  std::string *user_history_output_names_;

  // global constants
  Real four_pi_G_, grav_eps_, grav_mean_rho_;

  // functions
  MeshGenFunc_t MeshGenerator_[3];
  SrcTermFunc_t UserSourceTerm_;
  BValFunc_t BoundaryFunction_[6];
  AMRFlagFunc_t AMRFlag_;
  TimeStepFunc_t UserTimeStep_;
  HistoryOutputFunc_t *user_history_func_;
  MetricFunc_t UserMetric_;
  ViscosityCoeff_t ViscosityCoeff_;
  ConductionCoeff_t ConductionCoeff_;
  FieldDiffusionCoeff_t FieldDiffusivity_;
  MGBoundaryFunc_t MGBoundaryFunction_[6];

  void AllocateRealUserMeshDataField(int n);
  void AllocateIntUserMeshDataField(int n);
  void OutputMeshStructure(int dim);
  void LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb);

  // methods in ../pgen
  void InitUserMeshData(ParameterInput *pin);
  void EnrollUserBoundaryFunction (enum BoundaryFace face, BValFunc_t my_func);
  void EnrollUserRefinementCondition(AMRFlagFunc_t amrflag);
  void EnrollUserMeshGenerator(enum CoordinateDirection dir, MeshGenFunc_t my_mg);
  void EnrollUserExplicitSourceFunction(SrcTermFunc_t my_func);
  void EnrollUserTimeStepFunction(TimeStepFunc_t my_func);
  void AllocateUserHistoryOutput(int n);
  void EnrollUserHistoryOutput(int i, HistoryOutputFunc_t my_func, const char *name);
  void EnrollUserMetric(MetricFunc_t my_func);
  void EnrollUserMGBoundaryFunction(enum BoundaryFace dir, MGBoundaryFunc_t my_bc);
  void EnrollViscosityCoefficient(ViscosityCoeff_t my_func);
  void EnrollConductionCoefficient(ConductionCoeff_t my_func);
  void EnrollFieldDiffusivity(FieldDiffusionCoeff_t my_func);
  void SetGravitationalConstant(Real g) { four_pi_G_=4.0*PI*g; }
  void SetFourPiG(Real fpg) { four_pi_G_=fpg; }
  void SetGravityThreshold(Real eps) { grav_eps_=eps; }
  void SetMeanDensity(Real d0) { grav_mean_rho_=d0; }
};


//----------------------------------------------------------------------------------------
// \!fn Real ComputeMeshGeneratorX(int64_t index, int64_t nrange, bool sym_interval)
// \brief wrapper fn to compute Real x logical location for either [0, 1] or [-0.5, 0.5]
//        real cell ranges for MeshGenerator_[] functions (default/user vs. uniform)

inline Real ComputeMeshGeneratorX(int64_t index, int64_t nrange, bool sym_interval) {
  // index is typically 0, ... nrange for non-ghost boundaries
  if (sym_interval == false) {
    // to map to fractional logical position [0.0, 1.0], simply divide by # of faces
    return static_cast<Real>(index)/static_cast<Real>(nrange);
  } else {
    // to map to a [-0.5, 0.5] range, rescale int indices around 0 before FP conversion
    // if nrange is even, there is an index at center x=0.0; map it to (int) 0
    // if nrange is odd, the center x=0.0 is between two indices; map them to -1, 1
    int64_t noffset = index - (nrange)/2;
    int64_t noffset_ceil = index - (nrange+1)/2; // = noffset if nrange is even
    //std::cout << "noffset, noffset_ceil = " << noffset << ", " << noffset_ceil << "\n";
    // average the (possibly) biased integer indexing
    return static_cast<Real>(noffset + noffset_ceil)/(2.0*nrange);
  }
}

//----------------------------------------------------------------------------------------
// \!fn Real DefaultMeshGeneratorX1(Real x, RegionSize rs)
// \brief x1 mesh generator function, x is the logical location; x=i/nx1, real in [0, 1]

inline Real DefaultMeshGeneratorX1(Real x, RegionSize rs) {
  Real lw, rw;
  if (rs.x1rat==1.0) {
    rw=x, lw=1.0-x;
  } else {
    Real ratn=pow(rs.x1rat,rs.nx1);
    Real rnx=pow(rs.x1rat,x*rs.nx1);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  // linear interp, equally weighted from left (x(xmin)=0.0) and right (x(xmax)=1.0)
  return rs.x1min*lw+rs.x1max*rw;
}

//----------------------------------------------------------------------------------------
// \!fn Real DefaultMeshGeneratorX2(Real x, RegionSize rs)
// \brief x2 mesh generator function, x is the logical location; x=j/nx2, real in [0, 1]

inline Real DefaultMeshGeneratorX2(Real x, RegionSize rs) {
  Real lw, rw;
  if (rs.x2rat==1.0) {
    rw=x, lw=1.0-x;
  } else {
    Real ratn=pow(rs.x2rat,rs.nx2);
    Real rnx=pow(rs.x2rat,x*rs.nx2);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  return rs.x2min*lw+rs.x2max*rw;
}

//----------------------------------------------------------------------------------------
// \!fn Real DefaultMeshGeneratorX3(Real x, RegionSize rs)
// \brief x3 mesh generator function, x is the logical location; x=k/nx3, real in [0, 1]

inline Real DefaultMeshGeneratorX3(Real x, RegionSize rs) {
  Real lw, rw;
  if (rs.x3rat==1.0) {
    rw=x, lw=1.0-x;
  } else {
    Real ratn=pow(rs.x3rat,rs.nx3);
    Real rnx=pow(rs.x3rat,x*rs.nx3);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
  }
  return rs.x3min*lw+rs.x3max*rw;
}

//----------------------------------------------------------------------------------------
// \!fn Real UniformMeshGeneratorX1(Real x, RegionSize rs)
// \brief x1 mesh generator function, x is the logical location; real cells in [-0.5, 0.5]

inline Real UniformMeshGeneratorX1(Real x, RegionSize rs) {
  // linear interp, equally weighted from left (x(xmin)=-0.5) and right (x(xmax)=0.5)
  return static_cast<Real>(0.5)*(rs.x1min+rs.x1max) + (x*rs.x1max - x*rs.x1min);
}

//----------------------------------------------------------------------------------------
// \!fn Real UniformMeshGeneratorX2(Real x, RegionSize rs)
// \brief x2 mesh generator function, x is the logical location; real cells in [-0.5, 0.5]

inline Real UniformMeshGeneratorX2(Real x, RegionSize rs) {
  return static_cast<Real>(0.5)*(rs.x2min+rs.x2max) + (x*rs.x2max - x*rs.x2min);
}

//----------------------------------------------------------------------------------------
// \!fn Real UniformMeshGeneratorX3(Real x, RegionSize rs)
// \brief x3 mesh generator function, x is the logical location; real cells in [-0.5, 0.5]

inline Real UniformMeshGeneratorX3(Real x, RegionSize rs) {
  return static_cast<Real>(0.5)*(rs.x3min+rs.x3max) + (x*rs.x3max - x*rs.x3min);
}

#endif  // MESH_MESH_HPP_
