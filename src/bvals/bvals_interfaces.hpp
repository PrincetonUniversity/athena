#ifndef BVALS_BVALS_INTERFACES_HPP_
#define BVALS_BVALS_INTERFACES_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_interfaces.hpp
//! \brief defines enums, structs, and abstract classes
//!
//! \todo (felker):
//! - deduplicate forward declarations
//! - consider moving enums and structs in a new file? bvals_structs.hpp?

// C++ headers
#include <string>   // string
#include <vector>   // vector

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
class BoundaryValues;
struct RegionSize;
struct FaceField;

//! \todo (felker):
//! - nest these enum definitions inside bvals/ classes, when possible.

//! \deprecated (felker):
//! - maintain old-style (ALL_CAPS) enumerators as unscoped,unnamed types
//! - Keep for compatibility with user-provided pgen/ files.
//! - Use only new types internally.

// GCC 6 added Enumerator Attr (v6.1 released on 2016-04-27)
//! \todo (felker):
//! - replace with C++14 [[deprecated]] attributes if we ever bump --std=c++14
#if (defined(__GNUC__) &&__GNUC__ >= 6) || (defined(__clang__) && __clang_major__ >= 3)
enum {FACE_UNDEF __attribute__((deprecated)) = -1,
      INNER_X1 __attribute__((deprecated)),
      OUTER_X1 __attribute__((deprecated)),
      INNER_X2 __attribute__((deprecated)),
      OUTER_X2 __attribute__((deprecated)),
      INNER_X3 __attribute__((deprecated)),
      OUTER_X3 __attribute__((deprecated))};
enum {BLOCK_BNDRY __attribute__((deprecated)) = -1,
      BNDRY_UNDEF __attribute__((deprecated)),
      REFLECTING_BNDRY __attribute__((deprecated)),
      OUTFLOW_BNDRY __attribute__((deprecated)),
      USER_BNDRY __attribute__((deprecated)),
      PERIODIC_BNDRY __attribute__((deprecated)),
      POLAR_BNDRY __attribute__((deprecated)),
      POLAR_BNDRY_WEDGE __attribute__((deprecated)),
      SHEAR_PERIODIC_BNDRY __attribute__((deprecated))};
#else
enum {FACE_UNDEF = -1, INNER_X1, OUTER_X1, INNER_X2, OUTER_X2, INNER_X3, OUTER_X3};
enum {BLOCK_BNDRY = -1, BNDRY_UNDEF, REFLECTING_BNDRY, OUTFLOW_BNDRY, USER_BNDRY,
      PERIODIC_BNDRY, POLAR_BNDRY, POLAR_BNDRY_WEDGE, SHEAR_PERIODIC_BNDRY};
#endif

//! identifiers for all 6 faces of a MeshBlock
//! \todo (felker):
//! - BoundaryFace must be unscoped enum, for now. Its enumerators are used as
//!   int to index raw arrays (not AthenaArrays)
//!   --> enumerator vals are explicitly specified
enum BoundaryFace {undef=-1, inner_x1=0, outer_x1=1, inner_x2=2, outer_x2=3,
                   inner_x3=4, outer_x3=5};

//! identifiers for boundary conditions
enum class BoundaryFlag {block=-1, undef, reflect, outflow, user, periodic,
                         polar, polar_wedge, shear_periodic};

//! identifiers for types of neighbor blocks (connectivity with current MeshBlock)
enum class NeighborConnect {none, face, edge, corner}; // degenerate/shared part of block

//! identifiers for status of MPI boundary communications
enum class BoundaryStatus {waiting, arrived, completed};

//! flags to mark which variables are reversed across polar boundary
constexpr const bool flip_across_pole_hydro[] = {false, false, true, true, false};
//! flags to mark which variables are reversed across polar boundary
constexpr const bool flip_across_pole_field[] = {false, true, true};

//----------------------------------------------------------------------------------------
//! \struct SimpleNeighborBlock
//! \brief Struct storing only the basic info about a MeshBlocks neighbors.
//!
//! Typically used for convenience to store redundant info
//! from subset of the more complete NeighborBlock objects,
//! e.g. for describing neighbors around pole at same radius and polar angle

struct SimpleNeighborBlock { // aggregate and POD
  int rank;    //!< MPI rank of neighbor
  int level;   //!< refinement (logical, not physical) level of neighbor
  int lid;     //!< local ID of neighbor
  int gid;     //!< global ID of neighbor
};

//----------------------------------------------------------------------------------------
//! \struct NeighborIndexes
//! \brief data to describe MeshBlock neighbors
//! \note
//! User-provided ctor is unnecessary and prevents the type from being POD and aggregate.
//! This struct's implicitly-defined or defaulted default ctor is trivial, implying that
//! NeighborIndexes is a trivial type. Combined with standard layout --> POD. Advantages:
//!  - No user-provided ctor: value initialization first performs zero initialization
//!    (then default initialization if ctor is non-trivial)
//!  - Aggregate type: supports aggregate initialization {}
//!  - POD type: safely copy objects via memcpy, no memory padding in the beginning of
//!    object, C portability, supports static initialization

struct NeighborIndexes { // aggregate and POD
  int ox1, ox2, ox3; // 3-vec of offsets in {-1,0,+1} relative to this block's (i,j,k)
  int fi1, fi2;      // 2-vec for identifying refined neighbors (up to 4x face neighbors
                     // in 3D), entries in {0, 1}={smaller, larger} LogicalLcation::lxi
  NeighborConnect type;
};

//----------------------------------------------------------------------------------------
//! \struct NeighborBlock
//! \brief

struct NeighborBlock { // aggregate and POD type. Inheritance breaks standard-layout-> POD
                       // : SimpleNeighborBlock, NeighborIndexes {
  // composition:
  SimpleNeighborBlock snb;
  NeighborIndexes ni;

  int bufid, eid, targetid;
  BoundaryFace fid;
  bool polar; //!> flag indicating boundary is across a pole
  bool shear; //!> flag indicating boundary is attaching shearing periodic boundaries.
              //!> (used only in flux_correction_fc.cpp, for now)

  void SetNeighbor(int irank, int ilevel, int igid, int ilid, int iox1, int iox2,
                   int iox3, NeighborConnect itype, int ibid, int itargetid,
                   bool ipolar, bool ishear, int ifi1=0, int ifi2=0);
};

//----------------------------------------------------------------------------------------
//! \struct BoundaryData
//! \brief structure storing boundary information
//!
//! \todo (felker):
//! - consider renaming/be more specific--- what kind of data/info?
//!   one for each type of "BoundaryQuantity" corresponding to BoundaryVariable

template <int n = 56>
struct BoundaryData { // aggregate and POD (even when MPI_PARALLEL is defined)
  static constexpr int kMaxNeighbor = n;
  // KGF: "nbmax" only used in bvals_var.cpp, Init/DestroyBoundaryData()
  int nbmax;  //!> actual maximum number of neighboring MeshBlocks
  // currently, sflag[] is only used by Multgrid (send buffers are reused each stage in
  // red-black comm. pattern; need to check if they are available) and shearing box
  BoundaryStatus flag[kMaxNeighbor], sflag[kMaxNeighbor];
  Real *send[kMaxNeighbor], *recv[kMaxNeighbor];
#ifdef MPI_PARALLEL
  MPI_Request req_send[kMaxNeighbor], req_recv[kMaxNeighbor];
#endif
};

using ShearingBoundaryData = BoundaryData<4>;
using ShearingFluxBoundaryData = BoundaryData<3>;

//----------------------------------------------------------------------------------------
//! \struct ShearNeighborData
//! \brief structure storing shearing boundary information

template <int n = 4>
struct ShearNeighborData {
  static constexpr int kMaxNeighbor = n;
  SimpleNeighborBlock send_neighbor[kMaxNeighbor], recv_neighbor[kMaxNeighbor];
  int send_count[kMaxNeighbor], recv_count[kMaxNeighbor];
  int jmin_send[kMaxNeighbor], jmax_send[kMaxNeighbor];
  int jmin_recv[kMaxNeighbor], jmax_recv[kMaxNeighbor];
};

// Struct for describing blocks which touch the shearing-periodic boundaries
// struct ShearingBoundaryBlock {
//   int *igidlist, *ilidlist, *irnklist, *ilevlist;
//   int *ogidlist, *olidlist, *ornklist, *olevlist;
//   bool inner, outer;
// };

//----------------------------------------------------------------------------------------
// Interfaces = abstract classes containing ONLY pure virtual functions
//              Merely lists functions and their argument lists that must be implemented
//              in derived classes to form a somewhat-strict contract for functionality
//              (can always implement as a do-nothing/no-op function if a derived class
//              instance is the exception to the rule and does not use a particular
//              interface function)
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//! \class BoundaryCommunication
//! \brief contains methods for managing BoundaryStatus flags and MPI requests

class BoundaryCommunication {
 public:
  BoundaryCommunication() {}
  virtual ~BoundaryCommunication() {}
  //! create unique tags for each MeshBlock/buffer/quantity and initalize MPI requests:
  virtual void SetupPersistentMPI() = 0;
  //! call MPI_Start() on req_recv[]
  virtual void StartReceiving(BoundaryCommSubset phase) = 0;
  //! call MPI_Wait() on req_send[] and set flag[] to BoundaryStatus::waiting
  virtual void ClearBoundary(BoundaryCommSubset phase) = 0;

  virtual void StartReceivingShear(BoundaryCommSubset phase) = 0;
};

//----------------------------------------------------------------------------------------
//! \class BoundaryBuffer
//! \brief contains methods for managing MPI send/recvs and associated loads/stores from
//! communication buffers.
//!
//! \todo (felker):
//! - Merge with above BoundaryCommunication interface?

class BoundaryBuffer {
 public:
  BoundaryBuffer() {}
  virtual ~BoundaryBuffer() {}

  // universal buffer management methods for Cartesian grids (unrefined and SMR/AMR)
  virtual void SendBoundaryBuffers() = 0;
  virtual bool ReceiveBoundaryBuffers() = 0;
  // this next fn is used only during problem initialization in mesh.cpp:
  virtual void ReceiveAndSetBoundariesWithWait() = 0;
  virtual void SetBoundaries() = 0;

  virtual void SendFluxCorrection() = 0;
  virtual bool ReceiveFluxCorrection() = 0;

 protected:
  // universal buffer management methods for Cartesian grids (unrefined and SMR/AMR):
  virtual int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) = 0;
  virtual void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) = 0;

  // SMR/AMR-exclusive buffer management methods:
  virtual int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) = 0;
  virtual int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) = 0;
  virtual void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) = 0;
  virtual void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) = 0;

  // optional extensions: spherical-polar-like coordinates, shearing box, etc.:
  virtual void PolarBoundarySingleAzimuthalBlock() = 0;
};

//----------------------------------------------------------------------------------------
//! \class BoundaryPhysics
//! \brief defines methods for handling non-periodic domain limits, including:
//!    far-field/freestream (outflow, reflect), coordinate (spherical polar), etc.

class BoundaryPhysics {
 public:
  BoundaryPhysics() {}
  virtual ~BoundaryPhysics() {}

  //--------- prototypes for all required physical/coordinate BC functions ---------------
  virtual void ReflectInnerX1(Real time, Real dt,
                              int il, int jl, int ju, int kl, int ku, int ngh) = 0;
  virtual void ReflectOuterX1(Real time, Real dt,
                              int iu, int jl, int ju, int kl, int ku, int ngh) = 0;
  virtual void ReflectInnerX2(Real time, Real dt,
                              int il, int iu, int jl, int kl, int ku, int ngh) = 0;
  virtual void ReflectOuterX2(Real time, Real dt,
                              int il, int iu, int ju, int kl, int ku, int ngh) = 0;
  virtual void ReflectInnerX3(Real time, Real dt,
                              int il, int iu, int jl, int ju, int kl, int ngh) = 0;
  virtual void ReflectOuterX3(Real time, Real dt,
                              int il, int iu, int jl, int ju, int ku, int ngh) = 0;

  virtual void OutflowInnerX1(Real time, Real dt,
                              int il, int jl, int ju, int kl, int ku, int ngh) = 0;
  virtual void OutflowOuterX1(Real time, Real dt,
                              int iu, int jl, int ju, int kl, int ku, int ngh) = 0;
  virtual void OutflowInnerX2(Real time, Real dt,
                              int il, int iu, int jl, int kl, int ku, int ngh) = 0;
  virtual void OutflowOuterX2(Real time, Real dt,
                              int il, int iu, int ju, int kl, int ku, int ngh) = 0;
  virtual void OutflowInnerX3(Real time, Real dt,
                              int il, int iu, int jl, int ju, int kl, int ngh) = 0;
  virtual void OutflowOuterX3(Real time, Real dt,
                              int il, int iu, int jl, int ju, int ku, int ngh) = 0;

  virtual void PolarWedgeInnerX2(Real time, Real dt,
                                 int il, int iu, int jl, int kl, int ku, int ngh) = 0;
  virtual void PolarWedgeOuterX2(Real time, Real dt,
                                 int il, int iu, int ju, int kl, int ku, int ngh) = 0;
};

//----------------------------------------------------------------------------------------
// Abstract classes containing mix of pure virtual, virtual, and concrete functoins
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
//! \class BoundaryVariable (abstract)
//! \brief

class BoundaryVariable : public BoundaryCommunication, public BoundaryBuffer,
                         public BoundaryPhysics {
 public:
  explicit BoundaryVariable(MeshBlock *pmb);
  virtual ~BoundaryVariable() = default;

  // (usuallly the std::size_t unsigned integer type)
  std::vector<BoundaryVariable *>::size_type bvar_index;

  virtual int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) = 0;
  virtual int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) = 0;

  //!@{
  //! BoundaryBuffer public functions with shared implementations
  void SendBoundaryBuffers() override;
  bool ReceiveBoundaryBuffers() override;
  void ReceiveAndSetBoundariesWithWait() override;
  void SetBoundaries() override;
  //!@}

 protected:
  // deferred initialization of BoundaryData objects in derived class constructors
  BoundaryData<> bd_var_, bd_var_flcor_;
  // derived class dtors are also responsible for calling DestroyBoundaryData(bd_var_)

  MeshBlock *pmy_block_;   // ptr to MeshBlock containing this BoundaryVariable
  Mesh *pmy_mesh_;
  BoundaryValues *pbval_;  // ptr to BoundaryValues that aggregates these
                           // BoundaryVariable objects

  void CopyVariableBufferSameProcess(NeighborBlock& nb, int ssize);
  void CopyFluxCorrectionBufferSameProcess(NeighborBlock& nb, int ssize);

  void InitBoundaryData(BoundaryData<> &bd, BoundaryQuantity type);
  void DestroyBoundaryData(BoundaryData<> &bd);

  ShearingBoundaryData shear_bd_var_[2];
  ShearingFluxBoundaryData shear_bd_flux_[2];
  // TODO(felker): combine 4x Copy*SameProcess() functions
  void CopyShearBufferSameProcess(SimpleNeighborBlock& snb, int ssize, int bufid,
                                  bool upper);
  void CopyShearFluxSameProcess(SimpleNeighborBlock& snb, int ssize, int bufid,
                               bool upper);
  void SetCompletedFlagSameProcess(NeighborBlock& nb);
  // private:
};

#endif // BVALS_BVALS_INTERFACES_HPP_
