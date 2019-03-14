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

// TODO(felker): nest these enum definitions inside bvals/ classes, when possible.

// identifiers for all 6 faces of a MeshBlock
enum BoundaryFace {undef=-1, inner_x1=0, outer_x1=1, inner_x2=2, outer_x2=3,
                   inner_x3=4, outer_x3=5};
// TODO(felker): BoundaryFace must be unscoped enum, for now. Its enumerators are used as
// int to index regular arrays (not AthenaArrays). Hence enumerator values are specified.

// identifiers for boundary conditions
enum class BoundaryFlag {block=-1, undef, reflect, outflow, user, periodic,
                         polar, polar_wedge, shear_periodic};

// identifiers for types of neighbor blocks (connectivity with current MeshBlock)
enum class NeighborConnect {none, face, edge, corner}; // degenerate/shared part of block

// identifiers for status of MPI boundary communications
enum class BoundaryStatus {waiting, arrived, completed};

// flags to mark which variables are reversed across polar boundary
constexpr const bool flip_across_pole_hydro[] = {false, false, true, true, false};
constexpr const bool flip_across_pole_field[] = {false, true, true};

//----------------------------------------------------------------------------------------
//! \struct NeighborBlock
//  \brief neighbor rank, level, and ids

struct NeighborBlock { // not aggregate nor POD type
  int rank, level;
  int gid, lid;
  int ox1, ox2, ox3;
  int fi1, fi2;
  int bufid, eid, targetid;
  NeighborConnect type;
  BoundaryFace fid;
  bool polar; // flag indicating boundary is across a pole
  bool shear; // flag indicating boundary is attaching shearing periodic boundaries.
  NeighborBlock() : rank(-1), level(-1), gid(-1), lid(-1), ox1(-1), ox2(-1), ox3(-1),
                    fi1(-1), fi2(-1), bufid(-1), eid(-1), targetid(-1),
                    type(NeighborConnect::none), fid(BoundaryFace::undef), polar(false),
                    shear(false) {}
  void SetNeighbor(int irank, int ilevel, int igid, int ilid, int iox1, int iox2,
                   int iox3, NeighborConnect itype, int ibid, int itargetid,
                   bool ipolar, bool ishear, int ifi1, int ifi2);
};

//----------------------------------------------------------------------------------------
//! \struct PolarNeighborBlock
//  \brief Struct for describing neighbors around pole at same radius and polar angle

struct PolarNeighborBlock { // aggregate and POD
  int rank;    // MPI rank of neighbor
  int lid;     // local ID of neighbor
  int gid;     // global ID of neighbor
  bool north;  // flag that is true for North pole and false for South pole
};

//! \struct NeighborConnect
//  \brief data to describe MeshBlock neighbors
struct NeighborIndexes { // aggregate and POD
  int ox1, ox2, ox3; // 3-vector of integer offsets of indices, {-1, 0, +1}
  int fi1, fi2; // 2-vector for identifying refined neighbors, {0, 1}
  NeighborConnect type;
  // User-provided ctor is unnecessary and prevents the type from being POD and aggregate:
  // NeighborIndexes() {
  //   ox1=0; ox2=0; ox3=0; fi1=0; fi2=0;
  //   type=NeighborConnect::none;
  // }

  // This struct's implicitly-defined or defaulted default ctor is trivial, implying that
  // NeighborIndexes is a trivial type. Combined with standard layout --> POD. Advantages:

  // No user-provided ctor: value initialization first performs zero initialization (then
  // default initialization if ctor is non-trivial)

  // Aggregate type: supports aggregate initialization {}

  // POD type: safely copy objects via memcpy, no memory padding in the beginning of
  // object, C portability, supports static initialization
};

//! \struct BoundaryData
//  \brief structure storing boundary information
// TODO(felker): rename/be more specific--- what kind of data/info?
// one for each type of "BoundaryQuantity" corresponding to BoundaryVariable
struct BoundaryData { // aggregate and POD (even when MPI_PARALLEL is defined)
  int nbmax;  // actual maximum number of neighboring MeshBlocks (at most 56)
  BoundaryStatus flag[56];
  // "static int bufid[56]" corresponds to these (for all BoundaryData instances per
  // MeshBlock, e.g. one per BoundaryQuantity / variable
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
};

// KGF: shearing box
// Struct for describing blocks which touched the shearing-periodic boundaries
// struct ShearingBoundaryBlock {
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
  BoundaryCommunication() {}
  virtual ~BoundaryCommunication() {}

  // functions called only at the start of simulation in Mesh::Initialize(res_flag, pin)
  // TODO(felker): rename this function to disambiguate from mesh.cpp, and specify MPI
  virtual void Initialize() = 0; // setup persistent MPI requests
  virtual void StartReceivingForInit(bool cons_and_field) = 0; // Call MPI_Start()
  virtual void ClearBoundaryForInit(bool cons_and_field) = 0;

  // functions called only in task_list/ during timestepping
  // time: pmesh->time+dtstep, where dtstep is the delta t for current step
  virtual void StartReceivingAll(const Real time) = 0;
  virtual void ClearBoundaryAll() = 0;
};

//----------------------------------------------------------------------------------------
//! \class BoundaryBuffer
//  \brief contains methods for managing MPI send/recvs and associated loads/stores from
//  communication buffers
class BoundaryBuffer {
 public:
  BoundaryBuffer() {}
  virtual ~BoundaryBuffer() {}

  // universal buffer management methods for Cartesian grids (unrefined and SMR/AMR)
  virtual void SendBoundaryBuffers() = 0; // client-facing
  virtual bool ReceiveBoundaryBuffers() = 0; // client-facing
  // this next fn is used only during problem initialization in mesh.cpp:
  virtual void ReceiveAndSetBoundariesWithWait() = 0; // client-facing
  virtual void SetBoundaries() = 0; // client-facing

  // TODO(felker): handle the 6x unique Field-related flux correction functions
  virtual void SendFluxCorrection() = 0;
  virtual bool ReceiveFluxCorrection() = 0;

 private:
  // universal buffer management methods for Cartesian grids (unrefined and SMR/AMR)
  virtual int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) = 0;
  virtual void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) = 0;

  // SMR/AMR-exclusive buffer management methods
  virtual int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) = 0;
  virtual int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) = 0;
  virtual void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) = 0;
  virtual void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) = 0;

  // optional extensions: spherical-polar-like coordinates, shearing box, etc.
  virtual void PolarBoundarySingleAzimuthalBlock() = 0;

  // compare to PolarBoundarySingleAzimuthalBlockField(),
  //                      PolarBoundarySingleAzimuthalBlockEMF()
  // what about PolarBoundaryAverageField()?
};

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
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;

  virtual void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;
  virtual void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              int il, int iu, int jl, int ju,
                              int kl, int ku, int ngh) = 0;

  virtual void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 int il, int iu, int jl,
                                 int ju, int kl, int ku, int ngh) = 0;
  virtual void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 int il, int iu, int jl,
                                 int ju, int kl, int ku, int ngh) = 0;
  //virtual void ApplyUserDefinedBoundary
};

//----------------------------------------------------------------------------------------
// Abstract classes containing mix of pure virtual, virtual, and concrete functoins
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
//! \class BoundaryVariable (abstract)
//  \brief nodes in a linked list of BoundaryVariable derived class instances

// this class will mostly contain pure virtual functions
// could be eliminated, in favor of having the CellCenteredBoundaryVariable and
// FaceCenteredBoundaryVariable derived classes directly inheriting from the 3x interface
// classes. But then, won't be able to treat instances of those two classes as equivalent
// BoundaryVariable instances in a linked list (subtyping polymorphism)...

// BoundaryVariable will only provide default implementations of some member functions in
// the BoundaryCommunication interface.
class BoundaryVariable : public BoundaryCommunication, public BoundaryBuffer,
                         public BoundaryPhysics {
 public:
  explicit BoundaryVariable(MeshBlock *pmb);  // , BoundaryQuantity type);
  virtual ~BoundaryVariable() = default; // used to call DestroyBoundaryData(bd_var_)

  // KGF: this is usuallly the std::size_t unsigned integer type
  std::vector<BoundaryVariable *>::size_type bvar_index;

  virtual int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) = 0;
  virtual int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) = 0;

 protected:
  // KGF option 1:
  // deferred initialization of full BoundaryData objects in derived class destructor
  BoundaryData bd_var_, bd_var_flcor_;

  // KGF option 2:
  // BoundaryData *pbd_var_;
  // BoundaryData *pbd_var_flcor_;  // (make optional? not used for unrefined Hydro)
  // BoundaryQuantity btype_;

  BoundaryValues *pbval_;  // ptr to BoundaryValues containing this linked list
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this BoundaryVariable
  // KGF: clean up mixed/duplicated locations of pointers to mesh/ classes
  // KGF: pmy_mesh_=protected member of BoundaryBase, pmy_block_=private in BoundaryValues
  Mesh *pmy_mesh_;  // KGF: replace pbval->pmy_mesh_ usages in cc/ and fc/ function defs

  void CopyVariableBufferSameProcess(NeighborBlock& nb, int ssize);
  void CopyFluxCorrectionBufferSameProcess(NeighborBlock& nb, int ssize);

  void InitBoundaryData(BoundaryData &bd, BoundaryQuantity type);
  void DestroyBoundaryData(BoundaryData &bd);

  // void InitBoundaryData(BoundaryData &bd, BoundaryQuantity type);
  // void DestroyBoundaryData(BoundaryData &bd);
  // (originally were virtual functions and a part of the BoundaryCommunication interface)
  // The above 2x functions should probably be defined outside this abstract class. Could
  // split them both up and define in the concrete derived classes, like the rest of the
  // BoundaryCommunication interface. However, the BoundaryValues class (which realizes
  // the BoundaryCommunication interface) does not need these 2x functions.

  // Or, make them the constructor/destructor of the
  // BoundaryData struct? But InitBoundaryData must access many data members of BValues

 private:
};

#endif // BVALS_BVALS_INTERFACES_HPP_
