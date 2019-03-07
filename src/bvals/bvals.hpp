#ifndef BVALS_BVALS_HPP_
#define BVALS_BVALS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.hpp
//  \brief defines BoundaryBase, BoundaryValues classes used for setting BCs on all data

// C headers

// C++ headers
#include <string>   // string
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "./bvals_interfaces.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
// TODO(felker): how many of these foward declarations are actually needed now?
// Can #include "./bvals_interfaces.hpp" suffice?
class Mesh;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;
struct RegionSize;

// functions to return boundary flag given input string, and vice versa
enum BoundaryFlag GetBoundaryFlag(std::string input_string);
std::string GetBoundaryString(enum BoundaryFlag input_flag);

//----------------------------------------------------------------------------------------
//! \class BoundaryBase
//  \brief Base class for all the BoundaryValues classes

// Mid-2018, BoundaryBase (no virtual methods) is the parent class to 3x derived classes:
// BoundaryValues, GravityBoundaryValues, MGBoundaryValues. However, the key class members
// are all static; is there a reason for the other members to not be static?

// TODO(felker): describe its contents, and why it has so much in common with Multigrid:

// After redesign, it should only be the parent class to 2x classes: MGBoundaryValues
// and BoundaryValues (or new analog). It should probably be merged with BoundaryValues,
// and MGBoundaryValues should be made into a new BoundaryVariables derived class

// Notation question: boundary value vs. boundary condition vs. boundary function
// BE CONSISTENT. What should the classes be named?

// KGF: potential names = BoundaryLogic, BoundaryNeighbors, BoundaryBase
// Used in mesh.hpp, meshblock_tree.hpp (Friend class), bvals_mg/grav*
class BoundaryBase {
 public:
  BoundaryBase(Mesh *pm, LogicalLocation iloc, RegionSize isize,
               enum BoundaryFlag *input_bcs);
  virtual ~BoundaryBase();

  // Currently, these 2x static data members are shared among 3x possible derived class
  // instances: BoundaryValues, MGBoundaryValues, GravityBoundaryValues

  // 1x pair (neighbor index, buffer ID) per entire SET of separate variable buffers
  // (Hydro, Field, Passive Scalar, Gravity, etc.). Greedy allocation for worst-case
  // of refined 3D; only 26 entries needed/initialized if unrefined 3D, e.g.
  static NeighborIndexes ni[56];
  static int bufid[56];

  NeighborBlock neighbor[56];
  int nneighbor;
  int nblevel[3][3][3];
  LogicalLocation loc;
  enum BoundaryFlag block_bcs[6];
  PolarNeighborBlock *polar_neighbor_north, *polar_neighbor_south;

  static int CreateBvalsMPITag(int lid, int phys, int bufid);
  static int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
  static int BufferID(int dim, bool multilevel);
  static int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);

  void SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist);

 protected:
  // 1D refined or unrefined=2
  // 2D refined=12, unrefined=8
  // 3D refined=56, unrefined=26. Refinement adds: 3*6 faces + 1*12 edges = +30 neighbors
  static int maxneighbor_;

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
//  \brief centralized class for interacting with individual variable boundary data

// KGF: BoundaryShared, BoundaryCommon, BoundaryCoupled, (keep?) BoundaryValues
// BoundaryInterface (as in "centralized interface for interacting with
//        BoundaryVariables", but this conflicts with "OO interfaces" and

class BoundaryValues : public BoundaryBase, //public BoundaryPhysics,
                       public BoundaryCommunication {
 public:
  BoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs, ParameterInput *pin);
  ~BoundaryValues();

  // must repeat function declarations to override the pure virtual methods from
  // BoundaryCommunication parent class:

  // called in BoundaryValues() constructor/destructor:
  // void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) final;
  // void DestroyBoundaryData(BoundaryData &bd) final;

  // called before time-stepper:
  void Initialize(void) final; // setup MPI requests
  void StartReceivingForInit(bool cons_and_field) final;
  void ClearBoundaryForInit(bool cons_and_field) final;
  // called during time-stepper:
  void StartReceivingAll(const Real time) final;
  void ClearBoundaryAll(void) final;

  // functions unique to BoundaryValues. these do not exist in individual BoundaryVariable

  // these typically define a coupled interaction of these boundary variables
  // KGF: will need to access AthenaArray<Real> &bcdst = pfield->bcc
  void ApplyPhysicalBoundaries(const Real time, const Real dt);
  void ProlongateBoundaries(const Real time, const Real dt);

  // The following 2x methods are unique to the BoundaryValues class, and serve only to
  // check the user's configuration. Called in Mesh::Initialize() after processing
  // ParameterInput(), before pbval->Initialize() is called in this class.
  void CheckBoundary(void);
  void CheckPolarBoundaries(void);

  // doubly linked list of references to BoundaryVariable instances
  std::vector<BoundaryVariable *> bvars;  // (num_bvars)

 private:
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this BoundaryValues
  // Set in the BoundaryValues() constructor based on block_bcs = input_bcs:
  // (these could probably all be moved to the FaceCenteredBoundaryVariable for the
  // flux_correction_fc.cpp functions. None are used in cc/)
  int num_north_polar_blocks_, num_south_polar_blocks_;

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_;
  // KGF: rename "firsttime_". The variable switch is used in only 2x functions:
  // ReceiveEMFCorrection() and StartReceivingAll()

  // For spherical polar coordinates edge-case: if one MeshBlock wraps entirely around
  // (azimuthally) the pole, shift the k-axis by nx3/2 for cell- and face-centered
  // variables, & emf. Used in bvals_cc.cpp, bvals_fc.cpp. Calculated in BoundaryValues()
  AthenaArray<Real> azimuthal_shift_;

  // Store signed, but positive, integer corresponding to the next unused value to be used
  // as unique tag ID for a BoundaryVariable object's MPI communication (formerly "enum
  // Athena_MPI_Tag"). 5 bits of unsigned integer representation are currently reserved
  // for this "phys" part of the bitfield tag, making 0, ..., 31 legal values
  int bvars_next_tag_;

  BValFunc BoundaryFunction_[6];

  // Shearingbox (shared with Field and Hydro)
  // ShearingBoundaryBlock shbb_;  // shearing block properties: lists etc.
  // Real x1size_,x2size_,x3size_; // mesh_size.x1max-mesh_size.x1min etc. [Lx,Ly,Lz]
  // Real Omega_0_, qshear_;       // orbital freq and shear rate
  // int ShBoxCoord_;              // shearcoordinate type: 1 = xy (default), 2 = xz
  // int joverlap_;                // # of cells the shear runs over one block
  // Real ssize_;                  // # of ghost cells in x-z plane
  // Real eps_;                    // fraction part of the shear

  // int  send_inner_gid_[4], recv_inner_gid_[4]; // gid of meshblocks for communication
  // int  send_inner_lid_[4], recv_inner_lid_[4]; // lid of meshblocks for communication
  // int send_inner_rank_[4],recv_inner_rank_[4]; // rank of meshblocks for communication
  // int  send_outer_gid_[4], recv_outer_gid_[4]; // gid of meshblocks for communication
  // int  send_outer_lid_[4], recv_outer_lid_[4]; // lid of meshblocks for communication
  // int send_outer_rank_[4],recv_outer_rank_[4]; // rank of meshblocks for communication

  // temporary--- Added by @tomidakn on 2015-11-27 in f0f989f85f
  friend class Mesh;
  friend class BoundaryVariable;
  friend class CellCenteredBoundaryVariable;
  friend class FaceCenteredBoundaryVariable;
};
#endif // BVALS_BVALS_HPP_
