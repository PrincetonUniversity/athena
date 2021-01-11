#ifndef BVALS_BVALS_HPP_
#define BVALS_BVALS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.hpp
//! \brief defines BoundaryBase, BoundaryValues classes used for setting BCs on all data

// C headers

// C++ headers
#include <string>   // string
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "bvals_interfaces.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
//! \todo (felker):
//! - how many of these foward declarations are actually needed now?
//! - Can #include "./bvals_interfaces.hpp" suffice?
class Mesh;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;
struct RegionSize;

//! free functions to return boundary flag given input string
BoundaryFlag GetBoundaryFlag(const std::string& input_string);
//! free functions to return string given boundary flag
std::string GetBoundaryString(BoundaryFlag input_flag);
//! confirming that the MeshBlock's boundaries are all valid selections
void CheckBoundaryFlag(BoundaryFlag block_flag, CoordinateDirection dir);

//----------------------------------------------------------------------------------------
//! \class BoundaryBase
//! \brief Base class for all BoundaryValues classes (BoundaryValues and MGBoundaryValues)
//!
//! \note
//! - 1x pair (neighbor index, buffer ID) per entire SET of separate variable buffers
//!   (Hydro, Field, Passive Scalar, Gravity, etc.). Greedy allocation for worst-case
//!   of refined 3D; only 26 entries needed/initialized if unrefined 3D, e.g.
//!   - 1D refined or unrefined=2
//!   - 2D refined=12, unrefined=8
//!   - 3D refined=56, unrefined=26.
//!     Refinement adds: 3*6 faces + 1*12 edges = +30 neighbors

class BoundaryBase {
 public:
  BoundaryBase(Mesh *pm, LogicalLocation iloc, RegionSize isize,
               BoundaryFlag *input_bcs);
  virtual ~BoundaryBase();
  // 1x pair (neighbor index, buffer ID) per entire SET of separate variable buffers
  // (Hydro, Field, Passive Scalar, Gravity, etc.). Greedy allocation for worst-case
  // of refined 3D; only 26 entries needed/initialized if unrefined 3D, e.g.
  static NeighborIndexes ni[56];
  static int bufid[56];

  NeighborBlock neighbor[56];
  int nneighbor;
  int nblevel[3][3][3];
  LogicalLocation loc;
  BoundaryFlag block_bcs[6];

  static int CreateBvalsMPITag(int lid, int bufid, int phys);
  static int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
  static int BufferID(int dim, bool multilevel);
  static int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);

  void SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist);

 protected:
  // 1D refined or unrefined=2
  // 2D refined=12, unrefined=8
  // 3D refined=56, unrefined=26. Refinement adds: 3*6 faces + 1*12 edges = +30 neighbors
  static int maxneighbor_;

  // used only in fc/
  int num_north_polar_blocks_, num_south_polar_blocks_;
  SimpleNeighborBlock *polar_neighbor_north_, *polar_neighbor_south_;

  Mesh *pmy_mesh_;
  RegionSize block_size_;
  AthenaArray<Real> sarea_[2];

 private:
  // calculate 3x shared static data members when constructing only the 1st class instance
  // int maxneighbor_=BufferID() computes ni[] and then calls bufid[]=CreateBufferID()
  static bool called_;
};

//----------------------------------------------------------------------------------------
//! \class BoundaryValues
//! \brief centralized class for interacting with each individual variable boundary data
//!        (design pattern ~ mediator)

class BoundaryValues : public BoundaryBase, //public BoundaryPhysics,
                       public BoundaryCommunication {
 public:
  BoundaryValues(MeshBlock *pmb, BoundaryFlag *input_bcs, ParameterInput *pin);
  ~BoundaryValues();

  //! variable-length arrays of references to BoundaryVariable instances
  //! containing all BoundaryVariable instances:
  std::vector<BoundaryVariable *> bvars;
  //! subset of bvars that are exchanged in the main TimeIntegratorTaskList
  std::vector<BoundaryVariable *> bvars_main_int;
  //! subset of bvars that are exchanged in the SuperTimeStepTaskList
  std::vector<BoundaryVariable *> bvars_sts;

  // inherited functions (interface shared with BoundaryVariable objects):
  // ------
  // called before time-stepper:
  void SetupPersistentMPI() final; // setup MPI requests

  // called before and during time-stepper:
  void StartReceiving(BoundaryCommSubset phase) final {return;};
  void ClearBoundary(BoundaryCommSubset phase) final {return;};
  void StartReceivingShear(BoundaryCommSubset phase) final;
  void ComputeShear(const Real time_fc, const Real time_int);

  // non-inhertied / unique functions (do not exist in BoundaryVariable objects):
  // (these typically involve a coupled interaction of boundary variable/quantities)
  // ------
  void StartReceivingSubset(BoundaryCommSubset phase,
                            std::vector<BoundaryVariable *> bvars_subset);
  void ClearBoundarySubset(BoundaryCommSubset phase,
                           std::vector<BoundaryVariable *> bvars_subset);
  void ApplyPhysicalBoundaries(const Real time, const Real dt,
                               std::vector<BoundaryVariable *> bvars_subset);
  void ProlongateBoundaries(const Real time, const Real dt,
                            std::vector<BoundaryVariable *> bvars_subset);

  // compute the shear at each integrator stage
  //! \todo (felker):
  //! - consider making this fn private again if calling within StartRecv()
  void FindShearBlock(const Real time);

  //! safety check of user's boundary fns in Mesh::Initialize before SetupPersistentMPI()
  void CheckUserBoundaries();

  int AdvanceCounterPhysID(int num_phys);

 private:
  MeshBlock *pmy_block_;      // ptr to MeshBlock containing this BoundaryValues
  int nface_, nedge_;         // used only in fc/flux_correction_fc.cpp calculations

  //! if a BoundaryPhysics or user fn should be applied at each MeshBlock boundary
  //! false --> e.g. block, polar, (shear-) periodic boundaries
  bool apply_bndry_fn_[6]{};   // C++11: in-class initializer of non-static member
  // C++11: direct-list-initialization -> value init of array -> zero init of each scalar

  //! For spherical polar coordinates edge-case: if one MeshBlock wraps entirely around
  //! the pole (azimuthally), shift the k-axis by nx3/2 for cell- and face-centered
  //! variables, & emf, using this temporary 1D array.
  AthenaArray<Real> azimuthal_shift_;

  //! local counter for generating unique MPI tags for per-MeshBlock BoundaryVariable
  //! communication (subset of Mesh::next_phys_id_)
  int bvars_next_phys_id_;

  // Shearing box (shared with Field and Hydro)
  // KGF: remove the redundancies in these variables:
  int shearing_box; // flag for shearing box: 0 = none, 1: xy, 2: xz
  int joverlap_, joverlap_flux_; // # of cells the shear runs over one block
  Real ssize_;                   // # of ghost cells in x-z plane
  Real eps_, eps_flux_;          // fraction part of the shear
  Real qomL_;
  int xorder_, xgh_;
  AthenaArray<Real> pflux_;    // pencil buffer for remapping

  std::int64_t nblx2;
  //! it is possible for a MeshBlock to have is_shear={true, true},
  //! if it is the only block along x1
  bool is_shear[2]; // inner_x1=0, outer_x1=1
  SimpleNeighborBlock *shbb_[2];
  std::int64_t loc_shear[2];  // x1 LogicalLocation of block(s) on inner/outer shear bndry

  // tomo-ono: 3x arrays and 4x arrays are required for int and fc, respectively
  ShearNeighborData<4> sb_data_[2];
  ShearNeighborData<3> sb_flux_data_[2];

  // ProlongateBoundaries() wraps the following S/AMR-operations (within nneighbor loop):
  // (the next function is also called within 3x nested loops over nk,nj,ni)
  void RestrictGhostCellsOnSameLevel(const NeighborBlock& nb, int nk, int nj, int ni);
  void ApplyPhysicalBoundariesOnCoarseLevel(
      const NeighborBlock& nb, const Real time, const Real dt,
      int si, int ei, int sj, int ej, int sk, int ek,
      std::vector<BoundaryVariable *> bvars_subset);
  void ProlongateGhostCells(const NeighborBlock& nb,
                            int si, int ei, int sj, int ej, int sk, int ek);

  void DispatchBoundaryFunctions(
      MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
      int il, int iu, int jl, int ju, int kl, int ku, int ngh,
      AthenaArray<Real> &prim, FaceField &b, BoundaryFace face,
      std::vector<BoundaryVariable *> bvars_subset);

  void CheckPolarBoundaries();  // called in BoundaryValues() ctor

  // temporary--- Added by @tomidakn on 2015-11-27 in f0f989f85f
  //! \todo (felker):
  //! - consider removing this friendship designation
  friend class Mesh;
  //! currently, this class friendship is required for copying send/recv buffers between
  //! BoundaryVariable objects within different MeshBlocks on the same MPI rank:
  friend class BoundaryVariable;
  friend class FaceCenteredBoundaryVariable;  // needs nface_, nedge_, num_north/south_...
  //! \todo (felker):
  //! - consider removing these friendship designations:
  friend class CellCenteredBoundaryVariable;
  friend class HydroBoundaryVariable;  // needed for shearing box quantities
};
#endif // BVALS_BVALS_HPP_
