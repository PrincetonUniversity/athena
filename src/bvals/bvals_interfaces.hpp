#ifndef BVALS_BVALS_INTERFACES_HPP_
#define BVALS_BVALS_INTERFACES_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_interfaces.hpp
//  \brief defines enums, structs, and abstract bvals/ classes

// TODO(felker): deduplicate headers and forward declarations
// TODO(felker): consider moving enums and structs in a new file? bvals_structs.hpp?

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
// one for each type of "enum BoundaryType" corresponding to BoundaryVariable
typedef struct BoundaryData {
  int nbmax;  // actual maximum number of neighboring MeshBlocks (at most 56)
  enum BoundaryStatus flag[56];
  // "static int bufid[56]" corresponds to these (for all BoundaryData instances per
  // MeshBlock, e.g. one per BoundaryType / variable
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} BoundaryData;

// KGF: shearing box
// Struct for describing blocks which touched the shearing-periodic boundaries
// typedef struct ShearingBoundaryBlock {
//   int *igidlist, *ilidlist, *irnklist, *ilevlist;
//   int *ogidlist, *olidlist, *ornklist, *olevlist;
//   bool inner, outer; // inner=true if inner blocks
// } ShearingBoundaryBlock;
// end KGF

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
//  \brief contains methods for managing BoundaryStatus flags and MPI requests
//    alternative name ideas: BoundaryMemory, BoundaryMPI,
class BoundaryCommunication {
 public:
  // allocate
  BoundaryCommunication() {}
  virtual ~BoundaryCommunication() {}

  // functions called exclusively in the constructor/destructor of same class instance
  virtual void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) = 0;
  virtual void DestroyBoundaryData(BoundaryData &bd) = 0;

  // functions called only at the start of simulation in Mesh::Initialize(res_flag, pin)
  // TODO(felker): rename this function to disambiguate from mesh.cpp, and specify MPI
  virtual void Initialize(void) = 0; // setup MPI requests
  virtual void StartReceivingForInit(bool cons_and_field) = 0; // Call MPI_Start()
  virtual void ClearBoundaryForInit(bool cons_and_field) = 0;

  // functions called only in task_list/ during timestepping
  // time: pmesh->time+dtstep, where dtstep is the delta t for current step
  virtual void StartReceivingAll(const Real time) = 0;
  virtual void ClearBoundaryAll(void) = 0;
};

//----------------------------------------------------------------------------------------
//! \class BoundaryBuffer
//  \brief contains methods for managing buffers, MPI requests


//----------------------------------------------------------------------------------------
//! \class BoundaryPhysics
//  \brief defines methods for handling non-periodic domain limits, including:
//     far-field/freestream (outflow, reflect), coordinate (spherical polar), etc.
class BoundaryPhysics {
 public:
  BoundaryPhysics() {}
  virtual ~BoundaryPhysics() {}

  //-------------------- prototypes for all BC functions ---------------------------------
  // TODO(felker): somehow, convert "FaceField &b" argument to some more generic type:
  // "BoundaryArray &dst" that can be either AthenaArray or FaceField = struct 3x
  // AthenaArray. Not quite templating, since each derived class will only support a
  // single type, but all of the other function arguments should be identical.

  // A similar mechanism will likely be used for the generic "BoundaryVariable" abstract
  // base class. E.g. "BoundaryArray *dst = phydro->w" in the class constructor
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
};

//----------------------------------------------------------------------------------------
// Abstract classes containing mix of pure virtual, virtual, and concrete functoins
//----------------------------------------------------------------------------------------



#endif // BVALS_BVALS_INTERFACES_HPP_