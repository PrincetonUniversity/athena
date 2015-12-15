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

// identifiers for all 6 boundaries of a MeshBlock
enum BoundarySide {SIDE_UNDEF=-1, INNER_X1=0, OUTER_X1=1, INNER_X2=2, OUTER_X2=3, 
  INNER_X3=4, OUTER_X3=5};

//-------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);

typedef void (*BValHydro_t)(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                            int is, int ie, int js, int je, int ks, int ke);
typedef void (*BValField_t)(MeshBlock *pmb, Coordinates *pco, FaceField &buf,
                            int is, int ie, int js, int je, int ks, int ke);

//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  void Initialize(void);
  void StartReceivingForInit(void);
  void StartReceivingAll(void);

  void EnrollHydroBoundaryFunction (enum BoundarySide edge, BValHydro_t  my_bc);
  void EnrollFieldBoundaryFunction(enum BoundarySide edge, BValField_t my_bc);
  void CheckBoundary(void);

  int LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step);
  void SetHydroBoundarySameLevel(AthenaArray<Real> &dst, Real *buf, const NeighborBlock& nb);
  void SetHydroBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetHydroBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf, const NeighborBlock& nb);
  bool ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst, int step);
  void ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst, int step);
  void RestrictHydro(AthenaArray<Real> &src,
                     int si, int ei, int sj, int ej, int sk, int ek);

  void SendFluxCorrection(int step);
  bool ReceiveFluxCorrection(int step);

  int LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendFieldBoundaryBuffers(FaceField &src, int step);
  void SetFieldBoundarySameLevel(FaceField &dst, Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(FaceField &dst, Real *buf, const NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(FaceField &dst, int step);
  void ReceiveFieldBoundaryBuffersWithWait(FaceField &dst, int step);

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  void SendEMFCorrection(int step);
  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  bool ReceiveEMFCorrection(int step);

  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst, 
                            FaceField &bfdst, AthenaArray<Real> &bcdst);

  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
                               FaceField &bfdst, AthenaArray<Real> &bcdst);

  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);


private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  BValHydro_t HydroBoundary_[6];
  BValField_t FieldBoundary_[6];

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_[NSTEP];

  enum boundary_status hydro_flag_[NSTEP][56], field_flag_[NSTEP][56];
  enum boundary_status flcor_flag_[NSTEP][6][2][2];
  enum boundary_status emfcor_flag_[NSTEP][48];
  Real *hydro_send_[NSTEP][56],  *hydro_recv_[NSTEP][56];
  Real *field_send_[NSTEP][56],  *field_recv_[NSTEP][56];
  Real *flcor_send_[NSTEP][6],   *flcor_recv_[NSTEP][6][2][2];
  Real *emfcor_send_[NSTEP][48], *emfcor_recv_[NSTEP][48];
  AthenaArray<Real> sarea_[2];

#ifdef MPI_PARALLEL
  MPI_Request req_hydro_send_[NSTEP][56],  req_hydro_recv_[NSTEP][56];
  MPI_Request req_field_send_[NSTEP][56],  req_field_recv_[NSTEP][56];
  MPI_Request req_flcor_send_[NSTEP][6],   req_flcor_recv_[NSTEP][6][2][2];
  MPI_Request req_emfcor_send_[NSTEP][48], req_emfcor_recv_[NSTEP][48];
#endif
};

unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
unsigned int CreateMPITag(int lid, int flag, int phys, int bufid);
int BufferID(int dim, bool multilevel, bool face_only);
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax);

typedef struct NeighborIndexes
{
  int ox1, ox2, ox3, fi1, fi2;
  enum neighbor_type type;
} NeighborIndexes;


#endif // BOUNDARY_VALUES_HPP
