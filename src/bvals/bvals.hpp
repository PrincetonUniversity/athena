#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file bvals.hpp
//  \brief defines BoundaryValues class used for setting BCs on all data types
//======================================================================================

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
class MeshBlock;
class Hydro;
class Field;
class ParameterInput;
class Coordinates;
struct FaceField;
struct NeighborBlock;
struct PolarNeighborBlock;

// identifiers for all 6 faces of a MeshBlock
enum BoundaryFace {FACE_UNDEF=-1, INNER_X1=0, OUTER_X1=1, INNER_X2=2, OUTER_X2=3,
  INNER_X3=4, OUTER_X3=5};

// identifiers for boundary conditions
enum BoundaryFlag {BLOCK_BNDRY=-1, BNDRY_UNDEF=0, REFLECTING_BNDRY=1, OUTFLOW_BNDRY=2,
  USER_BNDRY=3, PERIODIC_BNDRY=4, POLAR_BNDRY=5};

// identifiers for types of neighbor blocks
enum NeighborType {NEIGHBOR_NONE, NEIGHBOR_FACE, NEIGHBOR_EDGE, NEIGHBOR_CORNER};

// identifiers for status of MPI boundary communications
enum BoundaryStatus {BNDRY_WAITING, BNDRY_ARRIVED, BNDRY_COMPLETED};

// flags to mark which variables are reversed across polar boundary
static bool flip_across_pole_hydro[] = {false, false, true, true, false};
static bool flip_across_pole_field[] = {false, true, true};

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

//-------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, int is, int ie, int js, int je, int ks, int ke);

// function to return boundary flag given input string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);

//[JMSHI
// Struct for describing blocks which touched the shearing-periodic boundaries
// For now, we caclulate the list info for every boundary blocks, which might
// cost some extra memory.
typedef struct ShearingBoundaryBlock {
  int *igidlist, *ilidlist, *irnklist, *ilevlist;
  int *ogidlist, *olidlist, *ornklist, *olevlist;
  bool inner, outer; //inner=true is inner blocks
  int cnid, cnrk; // corner block gid and rank
} ShearingBoundaryBlock;

//JMSHI]

//--------------------------------------------------------------------------------------
//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  static NeighborIndexes ni[56];
  static int bufid[56];

  void Initialize(void);
  void CheckBoundary(void);
  void StartReceivingForInit(void);
  void StartReceivingAll(void);
  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);
  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
                               FaceField &bfdst, AthenaArray<Real> &bcdst);
  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
                            FaceField &bfdst, AthenaArray<Real> &bcdst);

  int LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb, bool cons);
  int LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendHydroBoundaryBuffers(AthenaArray<Real> &src, bool cons);
  void SetHydroBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                 const NeighborBlock& nb);
  void SetHydroBoundaryFromCoarser(Real *buf, const NeighborBlock& nb, bool cons);
  void SetHydroBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                 const NeighborBlock& nb);
  bool ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst);
  void ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst, bool cons);
  void PolarSingleHydro(AthenaArray<Real> &dst);
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

  void SendFluxCorrection(void);
  bool ReceiveFluxCorrection(void);

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
//[JMSHI
  void LoadHydroShearing(AthenaArray<Real> &src, Real *buf, int nb);
  void SendHydroShearingboxBoundaryBuffersForInit(AthenaArray<Real> &src, bool cons);
  void SendHydroShearingboxBoundaryBuffers(AthenaArray<Real> &src, bool cons);

  void SetHydroShearingboxBoundarySameLevel(AthenaArray<Real> &dst, Real *buf, const int nb);
  bool ReceiveHydroShearingboxBoundaryBuffers(AthenaArray<Real> &dst);
  void ReceiveHydroShearingboxBoundaryBuffersWithWait(AthenaArray<Real> &dst, bool cons);
  void FindShearBlock(void);
  void RemapFlux(const int n, const int k, const int jinner, const int jouter, const int i, const Real eps, const AthenaArray<Real> &U, AthenaArray<Real> &Flux);
//JMSHI]

private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_;

  enum BoundaryStatus hydro_flag_[56], field_flag_[56];
  enum BoundaryStatus flcor_flag_[6][2][2];
  enum BoundaryStatus emfcor_flag_[48];
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real *hydro_send_[56],  *hydro_recv_[56];
  Real *field_send_[56],  *field_recv_[56];
  Real *flcor_send_[6],   *flcor_recv_[6][2][2];
  Real *emfcor_send_[48], *emfcor_recv_[48];
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;
  AthenaArray<Real> sarea_[2];
  AthenaArray<Real> exc_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

#ifdef MPI_PARALLEL
  MPI_Request req_hydro_send_[56],  req_hydro_recv_[56];
  MPI_Request req_field_send_[56],  req_field_recv_[56];
  MPI_Request req_flcor_send_[6],   req_flcor_recv_[6][2][2];
  MPI_Request req_emfcor_send_[48], req_emfcor_recv_[48];
  MPI_Request *req_emf_north_send_, *req_emf_north_recv_;
  MPI_Request *req_emf_south_send_, *req_emf_south_recv_;
#endif

  BValFunc_t BoundaryFunction_[6];

//[JMSHI
  enum BoundaryStatus shbox_inner_hydro_flag_[2], shbox_outer_hydro_flag_[2];
  Real x1size_,x2size_,x3size_; // mesh_size.x1max-mesh_size.x1min etc. [Lx,Ly,Lz]
  Real Omega_0_, qshear_; // orbital freq and shear rate
  int ShBoxCoord_; // shearcoordinate type: 1 = xy (default), 2 = xz
  Real ssize_; // # of ghost cells in x-z plane

  ShearingBoundaryBlock shbb_; // shearing block properties: lists etc.
  AthenaArray<Real> shboxvar_inner_hydro_, shboxvar_outer_hydro_; // working arrays of ghost zones with remapped quantities
  AthenaArray<Real> flx_inner_hydro_, flx_outer_hydro_; // flux from conservative remapping (fractional part of grid cells)

  int joverlap_; // # of cells the shear runs over one block
  Real eps_; // fraction part of the shear

  int send_inner_gid_[2],recv_inner_gid_[2]; // gid of meshblocks for communication
  int send_inner_lid_[2], recv_inner_lid_[2]; // lid of meshblocks for communication
  int send_inner_rank_[2],recv_inner_rank_[2]; // rank of meshblocks for communication
  int send_innersize_hydro_[2], recv_innersize_hydro_[2]; //MPI buffer sizes
  Real *send_innerbuf_hydro_[2], *recv_innerbuf_hydro_[2]; //MPI send and recv buffers
#ifdef MPI_PARALLEL
  MPI_Request rq_innersend_hydro_[2],  rq_innerrecv_hydro_[2];//MPI request for send and recv msgs
#endif
  int send_outer_gid_[2], recv_outer_gid_[2]; // gid of meshblocks for communication
  int send_outer_lid_[2], recv_outer_lid_[2]; // lid of meshblocks for communication
  int send_outer_rank_[2],recv_outer_rank_[2]; // rank of meshblocks for communication
  int send_outersize_hydro_[2], recv_outersize_hydro_[2]; //MPI buffer sizes
  Real *send_outerbuf_hydro_[2], *recv_outerbuf_hydro_[2]; //MPI send and recv buffers
#ifdef MPI_PARALLEL
  MPI_Request rq_outersend_hydro_[2],  rq_outerrecv_hydro_[2];//MPI request for send and recv msgs
#endif

//JMSHI]
  // temporary
  friend class Mesh;
};

unsigned int CreateBvalsMPITag(int lid, int phys, int bufid);
unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
int BufferID(int dim, bool multilevel);
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax);

#endif // BOUNDARY_VALUES_HPP
