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
  int rank, level, gid, lid, ox1, ox2, ox3, fi1, fi2, bufid, eid, targetid;
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
  int ox1, ox2, ox3, fi1, fi2;
  enum NeighborType type;
  NeighborIndexes() {
    ox1=0; ox2=0; ox3=0; fi1=0; fi2=0;
    type=NEIGHBOR_NONE;
  }
} NeighborIndexes;

//! \struct BoundaryData
//  \brief structure storing boundary information
typedef struct BoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} BoundaryData;


//---------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);

void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);


// function to return boundary flag given input string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);

// Struct for describing blocks which touched the shearing-periodic boundaries
typedef struct ShearingBoundaryBlock {
  int *igidlist, *ilidlist, *irnklist, *ilevlist;
  int *ogidlist, *olidlist, *ornklist, *olevlist;
  bool inner, outer; // inner=true if inner blocks
} ShearingBoundaryBlock;

//----------------------------------------------------------------------------------------
//! \class BoundaryBase
//  \brief Base class for all the BoundaryValues classes

class BoundaryBase {
public:
  BoundaryBase(Mesh *pm, LogicalLocation iloc, RegionSize isize,
               enum BoundaryFlag *input_bcs);
  virtual ~BoundaryBase();

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
  static int maxneighbor_;
  Mesh *pmy_mesh_;
  RegionSize block_size_;
  AthenaArray<Real> sarea_[2];

private:
  static bool called_;
};

//----------------------------------------------------------------------------------------
//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues : public BoundaryBase {
public:
  BoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs, ParameterInput *pin);
  ~BoundaryValues();

  void InitBoundaryData(BoundaryData &bd, enum BoundaryType type);
  void DestroyBoundaryData(BoundaryData &bd);
  void Initialize(void);
  void CheckBoundary(void);
  void StartReceivingForInit(bool cons_and_field);
  // time: pmesh->time+dtstep, where dtstep is the delta t for current step
  void StartReceivingAll(const Real time);
  void ClearBoundaryForInit(bool cons_and_field);
  void ClearBoundaryAll(void);
  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
       FaceField &bfdst, AthenaArray<Real> &bcdst, const Real time, const Real dt);
  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
       FaceField &bfdst, AthenaArray<Real> &bcdst, const Real time, const Real dt);

  int LoadCellCenteredBoundaryBufferSameLevel(AthenaArray<Real> &src,
                      int ns, int ne, Real *buf, const NeighborBlock& nb);
  int LoadCellCenteredBoundaryBufferToCoarser(AthenaArray<Real> &src,
      int ns, int ne, Real *buf, AthenaArray<Real> &cbuf, const NeighborBlock& nb);
  int LoadCellCenteredBoundaryBufferToFiner(AthenaArray<Real> &src,
                      int ns, int ne, Real *buf, const NeighborBlock& nb);
  void SendCellCenteredBoundaryBuffers(AthenaArray<Real> &src,
                                       enum CCBoundaryType type);
  void SetCellCenteredBoundarySameLevel(AthenaArray<Real> &dst, int ns, int ne,
                                  Real *buf, const NeighborBlock& nb, bool *flip);
  void SetCellCenteredBoundaryFromCoarser(int ns, int ne, Real *buf,
                      AthenaArray<Real> &cbuf, const NeighborBlock& nb, bool *flip);
  void SetCellCenteredBoundaryFromFiner(AthenaArray<Real> &dst, int ns, int ne,
                                  Real *buf, const NeighborBlock& nb, bool *flip);
  bool ReceiveCellCenteredBoundaryBuffers(AthenaArray<Real> &dst,
                                          enum CCBoundaryType type);
  void ReceiveCellCenteredBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                                           enum CCBoundaryType type);
  void PolarSingleCellCentered(AthenaArray<Real> &dst, int ns, int ne);

  int LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendFieldBoundaryBuffers(FaceField &src);
  void SetFieldBoundarySameLevel(FaceField &dst, Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(FaceField &dst, Real *buf, const NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(FaceField &dst);
  void ReceiveFieldBoundaryBuffersWithWait(FaceField &dst);
  void PolarSingleField(FaceField &dst);
  void PolarAxisFieldAverage(FaceField &dst);

  void SendFluxCorrection(enum FluxCorrectionType type);
  bool ReceiveFluxCorrection(enum FluxCorrectionType type);

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb);
  void SendEMFCorrection(void);
  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north);
  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  void PolarSingleEMF(void);
  bool ReceiveEMFCorrection(void);

  // Shearingbox Hydro
  void LoadHydroShearing(AthenaArray<Real> &src, Real *buf, int nb);
  void SendHydroShearingboxBoundaryBuffersForInit(AthenaArray<Real> &src, bool cons);
  void SendHydroShearingboxBoundaryBuffers(AthenaArray<Real> &src, bool cons);

  void SetHydroShearingboxBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                            const int nb);
  bool ReceiveHydroShearingboxBoundaryBuffers(AthenaArray<Real> &dst);
  void FindShearBlock(const Real time);
  void RemapFlux(const int n, const int k, const int jinner, const int jouter,
                 const int i, const Real eps, const AthenaArray<Real> &U,
                 AthenaArray<Real> &Flux);
  // Shearingbox Field
  void LoadFieldShearing(FaceField &src, Real *buf, int nb);
  void SendFieldShearingboxBoundaryBuffersForInit(FaceField &src, bool cons);
  void SendFieldShearingboxBoundaryBuffers(FaceField &src, bool cons);
  void SetFieldShearingboxBoundarySameLevel(FaceField &dst, Real *buf, const int nb);
  bool ReceiveFieldShearingboxBoundaryBuffers(FaceField &dst);
  void RemapFluxField(const int k, const int jinner, const int jouter, const int i,
                      const Real eps, const AthenaArray<Real> &U,
                      AthenaArray<Real> &Flux);
  // Shearingbox EMF
  void LoadEMFShearing(EdgeField &src, Real *buf, const int nb);
  void SendEMFShearingboxBoundaryCorrectionForInit(void);
  void SendEMFShearingboxBoundaryCorrection(void);
  void SetEMFShearingboxBoundarySameLevel(EdgeField &dst, Real *buf, const int nb);
  bool ReceiveEMFShearingboxBoundaryCorrection(void);
  void RemapEMFShearingboxBoundary(void);
  void ClearEMFShearing(EdgeField &work);
  void RemapFluxEMF(const int k, const int jinner, const int jouter, const Real eps,
                    const AthenaArray<Real> &U, AthenaArray<Real> &Flux);

private:
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this BVals
  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_;

  BoundaryData bd_hydro_, bd_field_, bd_flcor_, bd_emfcor_;
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;
  AthenaArray<Real> exc_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

#ifdef MPI_PARALLEL
  MPI_Request *req_emf_north_send_, *req_emf_north_recv_;
  MPI_Request *req_emf_south_send_, *req_emf_south_recv_;
#endif

  BValFunc_t BoundaryFunction_[6];

// Shearingbox
  ShearingBoundaryBlock shbb_;  // shearing block properties: lists etc.
  Real x1size_,x2size_,x3size_; // mesh_size.x1max-mesh_size.x1min etc. [Lx,Ly,Lz]
  Real Omega_0_, qshear_;       // orbital freq and shear rate
  int ShBoxCoord_;              // shearcoordinate type: 1 = xy (default), 2 = xz
  int joverlap_;                // # of cells the shear runs over one block
  Real ssize_;                  // # of ghost cells in x-z plane
  Real eps_;                    // fraction part of the shear
  int  send_inner_gid_[4], recv_inner_gid_[4]; // gid of meshblocks for communication
  int  send_inner_lid_[4], recv_inner_lid_[4]; // lid of meshblocks for communication
  int send_inner_rank_[4],recv_inner_rank_[4]; // rank of meshblocks for communication
  int  send_outer_gid_[4], recv_outer_gid_[4]; // gid of meshblocks for communication
  int  send_outer_lid_[4], recv_outer_lid_[4]; // lid of meshblocks for communication
  int send_outer_rank_[4],recv_outer_rank_[4]; // rank of meshblocks for communication

  // Hydro
  enum BoundaryStatus shbox_inner_hydro_flag_[4], shbox_outer_hydro_flag_[4];
  // working arrays of remapped quantities
  AthenaArray<Real>  shboxvar_inner_hydro_, shboxvar_outer_hydro_;
  // flux from conservative remapping
  AthenaArray<Real>  flx_inner_hydro_, flx_outer_hydro_;
  int  send_innersize_hydro_[4], recv_innersize_hydro_[4]; // buffer sizes
  Real *send_innerbuf_hydro_[4], *recv_innerbuf_hydro_[4]; // send and recv buffers
  int  send_outersize_hydro_[4], recv_outersize_hydro_[4]; // buffer sizes
  Real *send_outerbuf_hydro_[4], *recv_outerbuf_hydro_[4]; // send and recv buffers
#ifdef MPI_PARALLEL
  // MPI request for send and recv msgs
  MPI_Request rq_innersend_hydro_[4], rq_innerrecv_hydro_[4];
  MPI_Request rq_outersend_hydro_[4], rq_outerrecv_hydro_[4];
#endif
  // Field
  enum BoundaryStatus shbox_inner_field_flag_[4], shbox_outer_field_flag_[4];
  FaceField shboxvar_inner_field_, shboxvar_outer_field_;
  FaceField flx_inner_field_, flx_outer_field_;
  int  send_innersize_field_[4], recv_innersize_field_[4];
  Real *send_innerbuf_field_[4], *recv_innerbuf_field_[4];
  int  send_outersize_field_[4], recv_outersize_field_[4];
  Real *send_outerbuf_field_[4], *recv_outerbuf_field_[4];
#ifdef MPI_PARALLEL
  MPI_Request rq_innersend_field_[4], rq_innerrecv_field_[4];
  MPI_Request rq_outersend_field_[4], rq_outerrecv_field_[4];
#endif
  // EMF correction
  enum BoundaryStatus shbox_inner_emf_flag_[5], shbox_outer_emf_flag_[5];
  EdgeField shboxvar_inner_emf_, shboxvar_outer_emf_;
  EdgeField shboxmap_inner_emf_, shboxmap_outer_emf_;
  EdgeField flx_inner_emf_, flx_outer_emf_;
  int  send_innersize_emf_[4], recv_innersize_emf_[4];
  Real *send_innerbuf_emf_[4], *recv_innerbuf_emf_[4];
  int  send_outersize_emf_[4], recv_outersize_emf_[4];
  Real *send_outerbuf_emf_[4], *recv_outerbuf_emf_[4];
#ifdef MPI_PARALLEL
  MPI_Request rq_innersend_emf_[4],  rq_innerrecv_emf_[4];
  MPI_Request rq_outersend_emf_[4],  rq_outerrecv_emf_[4];
#endif

  // temporary
  friend class Mesh;
};

#endif // BVALS_BVALS_HPP_
