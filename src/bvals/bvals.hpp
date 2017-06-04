#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
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
  USER_BNDRY=3, PERIODIC_BNDRY=4, POLAR_BNDRY=5, POLAR_BNDRY_WEDGE=6};

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

//! \struct BoundaryData
//  \brief structure storing boundary information
typedef struct BoundaryData {
  enum BoundaryStatus flag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} BoundaryData;

//---------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


// function to return boundary flag given input string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);


//----------------------------------------------------------------------------------------
//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  static NeighborIndexes ni[56];
  static int bufid[56];

  void InitBoundaryData(BoundaryData &bd, enum BoundaryType type);
  void DestroyBoundaryData(BoundaryData &bd);
  void Initialize(void);
  void CheckBoundary(void);
  void StartReceivingForInit(bool cons_and_field);
  void StartReceivingAll(void);
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

// Gravity
  int LoadGravityBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb);
  void SendGravityBoundaryBuffers(AthenaArray<Real> &src);
  void SetGravityBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                 const NeighborBlock& nb);
  bool ReceiveGravityBoundaryBuffers(AthenaArray<Real> &dst);
  void ReceiveGravityBoundaryBuffersWithWait(AthenaArray<Real> &dst);

// Multigrid
  void StartReceivingMultigrid(int nc, enum MGBoundaryType type);
  void ClearBoundaryMultigrid(enum MGBoundaryType type);
  int LoadMultigridBoundaryBufferSameLevel(AthenaArray<Real> &src,
                   int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb);
  void SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
                                    int nc, enum MGBoundaryType type);
  void SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
                   int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb);
  bool ReceiveCellCenteredBoundaryBuffers(AthenaArray<Real> &dst,
                                          int nc, enum MGBoundaryType type);

private:
  MeshBlock *pmy_block_;  // ptr to MeshBlock containing this BVals

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_;

  BoundaryData bd_hydro_, bd_field_, bd_gravity_, bd_mggrav_;
  enum BoundaryStatus flcor_flag_[6][2][2];
  enum BoundaryStatus emfcor_flag_[48];
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real *flcor_send_[6],   *flcor_recv_[6][2][2];
  Real *emfcor_send_[48], *emfcor_recv_[48];
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;
  AthenaArray<Real> sarea_[2];
  AthenaArray<Real> exc_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

#ifdef MPI_PARALLEL
  MPI_Request req_flcor_send_[6],   req_flcor_recv_[6][2][2];
  MPI_Request req_emfcor_send_[48], req_emfcor_recv_[48];
  MPI_Request *req_emf_north_send_, *req_emf_north_recv_;
  MPI_Request *req_emf_south_send_, *req_emf_south_recv_;
#endif

  BValFunc_t BoundaryFunction_[6];

  // temporary
  friend class Mesh;
};

unsigned int CreateBvalsMPITag(int lid, int phys, int bufid);
unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
int BufferID(int dim, bool multilevel);
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax);

#endif // BOUNDARY_VALUES_HPP
