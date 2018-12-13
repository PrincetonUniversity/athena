#ifndef BVALS_BVALS_HPP_
#define BVALS_BVALS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals.hpp
//  \brief defines BoundaryValues class used for setting BCs on all data types

// C++ headers
#include <string>   // string

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class Mesh;
class MeshBlock;
class MeshBlockTree;
class Hydro;
class Field;
class ParameterInput;
class Coordinates;
struct RegionSize;
struct FaceField;

// identifiers for all 6 faces of a MeshBlock
enum BoundaryFace {FACE_UNDEF=-1, INNER_X1=0, OUTER_X1=1, INNER_X2=2, OUTER_X2=3,
                   INNER_X3=4, OUTER_X3=5};

// identifiers for boundary conditions
enum BoundaryFlag {BLOCK_BNDRY=-1, BNDRY_UNDEF=0, REFLECTING_BNDRY=1, OUTFLOW_BNDRY=2,
                   USER_BNDRY=3, PERIODIC_BNDRY=4, POLAR_BNDRY=5, POLAR_BNDRY_WEDGE=6,
                   SHEAR_PERIODIC_BNDRY=7};

// identifiers for types of neighbor blocks
enum NeighborType {NEIGHBOR_NONE=0, NEIGHBOR_FACE=1, NEIGHBOR_EDGE=2, NEIGHBOR_CORNER=3};

// identifiers for status of MPI boundary communications
enum BoundaryStatus {BNDRY_WAITING, BNDRY_ARRIVED, BNDRY_COMPLETED};

// flags to mark which variables are reversed across polar boundary
static bool flip_across_pole_hydro[] = {false, false, true, true, false};
static bool flip_across_pole_field[] = {false, true, true};

//----------------------------------------------------------------------------------------
//! \struct NeighborBlock
//  \brief neighbor rank, level, and ids

typedef struct NeighborBlock {
  int rank, level;
  int gid, lid;
  int ox1, ox2, ox3;
  int fi1, fi2;
  int bufid, eid, targetid;
  enum NeighborType type;
  enum BoundaryFace fid;
  bool polar; // flag indicating boundary is across a pole
  bool shear; // flag indicating boundary is attaching shearing periodic boundaries.
  NeighborBlock() : rank(-1), level(-1), gid(-1), lid(-1), ox1(-1), ox2(-1), ox3(-1),
                    fi1(-1), fi2(-1), bufid(-1), eid(-1), targetid(-1),
                    type(NEIGHBOR_NONE), fid(FACE_UNDEF), polar(false), shear(false) {}
  void SetNeighbor(int irank, int ilevel, int igid, int ilid, int iox1, int iox2,
                   int iox3, enum NeighborType itype, int ibid, int itargetid,
                   bool ipolar, bool ishear, int ifi1, int ifi2);
} NeighborBlock;

//----------------------------------------------------------------------------------------
//! \struct PolarNeighborBlock
//  \brief Struct for describing neighbors around pole at same radius and polar angle

typedef struct PolarNeighborBlock {
  int rank;    // MPI rank of neighbor
  int lid;     // local ID of neighbor
  int gid;     // global ID of neighbor
  bool north;  // flag that is true for North pole and false for South pole
} PolarNeighborBlock;

//! \struct NeighborType
//  \brief data to describe MeshBlock neighbors
typedef struct NeighborIndexes {
  int ox1, ox2, ox3; // 3-vector of integer offsets of indices, {-1, 0, +1}
  int fi1, fi2; // 2-vector for identifying refined neighbors, {0, 1}
  enum NeighborType type;
  NeighborIndexes() {
    ox1=0; ox2=0; ox3=0; fi1=0; fi2=0;
    type=NEIGHBOR_NONE;
  }
} NeighborIndexes;

//! \struct BoundaryData
//  \brief structure storing boundary information
// TODO(felker): rename/be more specific--- what kind of data/info?
typedef struct BoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} BoundaryData;

// function to return boundary flag given user's athinput string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);

// KGF: shearing box
// Struct for describing blocks which touched the shearing-periodic boundaries
// typedef struct ShearingBoundaryBlock {
//   int *igidlist, *ilidlist, *irnklist, *ilevlist;
//   int *ogidlist, *olidlist, *ornklist, *olevlist;
//   bool inner, outer; // inner=true if inner blocks
// } ShearingBoundaryBlock;
// end KGF

//----------------------------------------------------------------------------------------
//! \class BoundaryBase
//  \brief Base class for all the BoundaryValues classes

// Mid-2018, BoundaryBase (no virtual methods) is the parent class to 3x derived classes:
// BoundaryValues, GravityBoundaryValues, MGBoundaryValues

// After redesign, it should only be the parent class to 2x classes: MGBoundaryValues
// and BoundaryValues (or new analog)

// TODO(felker): describe its contents, and why it has so much in common with Multigrid:

// Notation question: boundary value vs. boundary condition vs. boundary function
// BE CONSISTENT. What should the classes be named?

// TOOD(felker): convert all class access specifiers (public:, private:, protected:) to
// 1-space indent, following Google C++ style

// KGF: BoundaryLogic, BoundaryNeighbors, BoundaryBase
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

  static unsigned int CreateBvalsMPITag(int lid, int phys, int bufid);
  static unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
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
//! \class BoundaryMemory
//  \brief interface = class containing only pure virtual functions
//                     contains methods for managing buffers, MPI requests
//    BoundaryCommunication

class BoundaryMemory {
 public:
  BoundaryMemory() {}
  virtual ~BoundaryMemory() {}

  // functions called exclusively in the constructor/destructor of same class instance
  virtual void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) = 0;
  virtual void DestroyBoundaryData(BoundaryData &bd) = 0;

  // functions called only at the start of simulation in Mesh::Initialize(res_flag, pin)
  // TODO(felker): rename this function to disambiguate from mesh.cpp, and specify MPI
  virtual void Initialize(void) = 0; // setup MPI requests
  virtual void StartReceivingForInit(bool cons_and_field) = 0;
  virtual void ClearBoundaryForInit(bool cons_and_field) = 0;

  // functions called only in task_list/ during timestepping
  // time: pmesh->time+dtstep, where dtstep is the delta t for current step
  virtual void StartReceivingAll(const Real time) = 0;
  virtual void ClearBoundaryAll(void) = 0;
};


//----------------------------------------------------------------------------------------
//! \class BoundaryValues
//  \brief centralized class for interacting with individual variable boundary data

// KGF: BoundaryShared, BoundaryCommon, BoundaryCoupled,
// BoundaryInterface (as in "centralized interface for interacting with
//                    BoundaryVariables", but this
class BoundaryValues : public BoundaryBase, public BoundaryMemory {
 public:
  BoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs, ParameterInput *pin);
  ~BoundaryValues();

  // override the pure virtual methods from BoundaryMemory parent class:
  // BoundaryValues() constructor
  void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) final;
  void DestroyBoundaryData(BoundaryData &bd) final;
  // before time-stepper:
  void Initialize(void) final; // setup MPI requests
  void StartReceivingForInit(bool cons_and_field) final;
  void ClearBoundaryForInit(bool cons_and_field) final;
  // during time-stepper:
  void StartReceivingAll(const Real time) final;
  void ClearBoundaryAll(void) final;

  // new functions that do not exist in individual variable boundary value classes
  // and define a coupled interaction of these boundary values
  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
       FaceField &bfdst, AthenaArray<Real> &bcdst, const Real time, const Real dt);
  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
       FaceField &bfdst, AthenaArray<Real> &bcdst, const Real time, const Real dt);
  // called in Mesh::Initialize() after processing ParameterInput(), before
  // pbval->Initialize() is called in this class
  void CheckBoundary(void);

 private:
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this BVals
  // Set in the BoundaryValues() constructor based on block_bcs = input_bcs:
  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

  // For spherical polar coordinates edge-case: if one MeshBlock wraps entirely around
  // (azimuthally) the pole, shift the k-axis by nx3/2 for cell- and face-centered
  // variables, & emf. Used in bvals_cc.cpp, bvals_fc.cpp. Calculated in BoundaryValues()
  AthenaArray<Real> azimuthal_shift_;

  BValFunc_t BoundaryFunction_[6];

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

  // temporary--- Added by @tomidakn on 2015-11-27
  friend class Mesh;
};

//----------------------------------------------------------------------------------------
//! \class CellCenteredBoundaryFunctions
//  \brief abstract base class for all the derived classes for applying BC to AthenaArray

class CellCenteredBoundaryFunctions : public BoundaryBase {
  // Allow functions to access most any variable to allow for fully general BC dependence
  friend class Hydro;
  friend class Field;
public:
  CellCenteredBoundaryFunctions();
  virtual ~CellCenteredBoundaryFunctions();

  // standard cell-centered buffer management:
  // required: unrefined
  int LoadCellCenteredBoundaryBufferSameLevel(AthenaArray<Real> &src,
                      int ns, int ne, Real *buf, const NeighborBlock& nb);
  // optional: SMR/AMR
  int LoadCellCenteredBoundaryBufferToCoarser(AthenaArray<Real> &src,
      int ns, int ne, Real *buf, AthenaArray<Real> &cbuf, const NeighborBlock& nb);
  int LoadCellCenteredBoundaryBufferToFiner(AthenaArray<Real> &src,
                      int ns, int ne, Real *buf, const NeighborBlock& nb);

  // required: universal
  void SendCellCenteredBoundaryBuffers(AthenaArray<Real> &src,
                                       enum CCBoundaryType type);

  // required: universal
  bool ReceiveCellCenteredBoundaryBuffers(AthenaArray<Real> &dst,
                                          enum CCBoundaryType type);

  // optional: initialization in mesh.cpp
  void ReceiveCellCenteredBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                                           enum CCBoundaryType type);

  // required: unrefined
  void SetCellCenteredBoundarySameLevel(AthenaArray<Real> &dst, int ns, int ne,
                                        Real *buf, const NeighborBlock& nb, bool *flip);

  // optional: SMR/AMR
  void SetCellCenteredBoundaryFromCoarser(int ns, int ne, Real *buf,
                                          AthenaArray<Real> &cbuf,
                                          const NeighborBlock& nb, bool *flip);
  void SetCellCenteredBoundaryFromFiner(AthenaArray<Real> &dst, int ns, int ne,
                                        Real *buf, const NeighborBlock& nb, bool *flip);

  // optional? compare to PolarSingleField(), PolarSingleEMF()
  // what about PolarAxisFieldAverage()
  void PolarSingleCellCentered(AthenaArray<Real> &dst, int ns, int ne);

  // Cell-centered flux correction functions are much simpler than Field counterpart
  // In addition to 2x simple Send/Recv EMFCorrection() functions, there are:
  // - 6x Load/Set EMF (not correction). No Load to finer, to Set to coarser, but
  //   Set/LoadEMFBoundaryPolarBuffer()
  // - AverageEMFBoundary(), ClearCoarseEMFBoundary(), PolarSingleEMF()
  void SendFluxCorrection(enum FluxCorrectionType type);
  bool ReceiveFluxCorrection(enum FluxCorrectionType type);
  // TODO(felker): FLUX_HYDRO=0 is the only defined FluxCorrectionType enum in athena.hpp


  // Shearingbox Hydro
  // void LoadHydroShearing(AthenaArray<Real> &src, Real *buf, int nb);
  // void SendHydroShearingboxBoundaryBuffersForInit(AthenaArray<Real> &src, bool cons);
  // void SendHydroShearingboxBoundaryBuffers(AthenaArray<Real> &src, bool cons);

  // void SetHydroShearingboxBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
  //                                           const int nb);
  // bool ReceiveHydroShearingboxBoundaryBuffers(AthenaArray<Real> &dst);
  // void FindShearBlock(const Real time);
  // void RemapFlux(const int n, const int k, const int jinner, const int jouter,
  //                const int i, const Real eps, const AthenaArray<Real> &U,
  //                AthenaArray<Real> &Flux);

  //-------------------- prototypes for all BC functions ---------------------------------
  virtual void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);

  virtual void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              AthenaArray<Real> &arr_cc, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);

  virtual void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 AthenaArray<Real> &arr_cc, int il, int iu, int jl,
                                 int ju, int kl, int ku, int nl, int nu);
  virtual void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 AthenaArray<Real> &arr_cc, int il, int iu, int jl,
                                 int ju, int kl, int ku, int nl, int nu);
 private:
  // standard cell-centered and flux BV private variables
  BoundaryData bd_hydro_, bd_flcor_;

  // Shearingbox Hydro
//   enum BoundaryStatus shbox_inner_hydro_flag_[4], shbox_outer_hydro_flag_[4];
//   // working arrays of remapped quantities
//   AthenaArray<Real>  shboxvar_inner_hydro_, shboxvar_outer_hydro_;
//   // Hydro flux from conservative remapping
//   AthenaArray<Real>  flx_inner_hydro_, flx_outer_hydro_;
//   int  send_innersize_hydro_[4], recv_innersize_hydro_[4]; // buffer sizes
//   Real *send_innerbuf_hydro_[4], *recv_innerbuf_hydro_[4]; // send and recv buffers
//   int  send_outersize_hydro_[4], recv_outersize_hydro_[4]; // buffer sizes
//   Real *send_outerbuf_hydro_[4], *recv_outerbuf_hydro_[4]; // send and recv buffers
// #ifdef MPI_PARALLEL
//   // MPI request for send and recv msgs
//   MPI_Request rq_innersend_hydro_[4], rq_innerrecv_hydro_[4];
//   MPI_Request rq_outersend_hydro_[4], rq_outerrecv_hydro_[4];
// #endif

  //protected:

};

//----------------------------------------------------------------------------------------
//! \class FieldBoundaryFunctions
//  \brief abstract base class for all the derived classes for applying BC to AthenaArray

class FieldBoundaryFunctions {
  // Allow functions to access most any variable to allow for fully general BC dependence
  friend class Hydro;
  friend class Field;

 public:
  FieldBoundaryFunctions();
  virtual ~FieldBoundaryFunctions();
  // standard Field buffer management (unrefined and SMR/AMR)
  int LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                     const NeighborBlock& nb);

  void SendFieldBoundaryBuffers(FaceField &src);

  bool ReceiveFieldBoundaryBuffers(FaceField &dst);
  void ReceiveFieldBoundaryBuffersWithWait(FaceField &dst);

  void SetFieldBoundarySameLevel(FaceField &dst, Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(FaceField &dst, Real *buf, const NeighborBlock& nb);

  // above, only PolarSingleCellCentered(), no PolarAxis*Average()
  void PolarSingleField(FaceField &dst);
  void PolarAxisFieldAverage(FaceField &dst);

  // EMF buffer management (unrefined and SMR/AMR)
  // Hydro FluxCorrection only has counterparts for these
  void SendEMFCorrection(void);
  bool ReceiveEMFCorrection(void);

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb);

  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north);

  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  void PolarSingleEMF(void);

  // Shearingbox Field
  // void LoadFieldShearing(FaceField &src, Real *buf, int nb);
  // void SendFieldShearingboxBoundaryBuffersForInit(FaceField &src, bool cons);
  // void SendFieldShearingboxBoundaryBuffers(FaceField &src, bool cons);
  // void SetFieldShearingboxBoundarySameLevel(FaceField &dst, Real *buf, const int nb);
  // bool ReceiveFieldShearingboxBoundaryBuffers(FaceField &dst);
  // void RemapFluxField(const int k, const int jinner, const int jouter, const int i,
  //                     const Real eps, const AthenaArray<Real> &U,
  //                     AthenaArray<Real> &Flux);
  // // Shearingbox EMF
  // void LoadEMFShearing(EdgeField &src, Real *buf, const int nb);
  // void SendEMFShearingboxBoundaryCorrectionForInit(void);
  // void SendEMFShearingboxBoundaryCorrection(void);
  // void SetEMFShearingboxBoundarySameLevel(EdgeField &dst, Real *buf, const int nb);
  // bool ReceiveEMFShearingboxBoundaryCorrection(void);
  // void RemapEMFShearingboxBoundary(void);
  // void ClearEMFShearing(EdgeField &work);
  // void RemapFluxEMF(const int k, const int jinner, const int jouter, const Real eps,
  //                   const AthenaArray<Real> &U, AthenaArray<Real> &Flux);

  //-------------------- prototypes for all BC functions ---------------------------------
  virtual void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);

  virtual void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);
  virtual void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nl, int nu);

  virtual void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 FaceField &b, int il, int iu, int jl,
                                 int ju, int kl, int ku, int nl, int nu);
  virtual void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 FaceField &b, int il, int iu, int jl,
                                 int ju, int kl, int ku, int nl, int nu);
 private:
  // standard Field and emf BV private variables
  BoundaryData bd_field_, bd_emfcor_;
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;

#ifdef MPI_PARALLEL
  MPI_Request *req_emf_north_send_, *req_emf_north_recv_;
  MPI_Request *req_emf_south_send_, *req_emf_south_recv_;
#endif

  // Shearingbox Field
//   enum BoundaryStatus shbox_inner_field_flag_[4], shbox_outer_field_flag_[4];
//   FaceField shboxvar_inner_field_, shboxvar_outer_field_;
//   FaceField flx_inner_field_, flx_outer_field_;
//   int  send_innersize_field_[4], recv_innersize_field_[4];
//   Real *send_innerbuf_field_[4], *recv_innerbuf_field_[4];
//   int  send_outersize_field_[4], recv_outersize_field_[4];
//   Real *send_outerbuf_field_[4], *recv_outerbuf_field_[4];
// #ifdef MPI_PARALLEL
//   MPI_Request rq_innersend_field_[4], rq_innerrecv_field_[4];
//   MPI_Request rq_outersend_field_[4], rq_outerrecv_field_[4];
// #endif
//   // Shearing box EMF correction
//   enum BoundaryStatus shbox_inner_emf_flag_[5], shbox_outer_emf_flag_[5];
//   EdgeField shboxvar_inner_emf_, shboxvar_outer_emf_;
//   EdgeField shboxmap_inner_emf_, shboxmap_outer_emf_;
//   EdgeField flx_inner_emf_, flx_outer_emf_;
//   int  send_innersize_emf_[4], recv_innersize_emf_[4];
//   Real *send_innerbuf_emf_[4], *recv_innerbuf_emf_[4];
//   int  send_outersize_emf_[4], recv_outersize_emf_[4];
//   Real *send_outerbuf_emf_[4], *recv_outerbuf_emf_[4];
// #ifdef MPI_PARALLEL
//   MPI_Request rq_innersend_emf_[4],  rq_innerrecv_emf_[4];
//   MPI_Request rq_outersend_emf_[4],  rq_outerrecv_emf_[4];
// #endif

  // protected:
};


#endif // BVALS_BVALS_HPP_
