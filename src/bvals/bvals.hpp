#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file bvals.hpp
 *  \brief defines BoundaryValues class used for setting BCs on all data types
 *====================================================================================*/

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class MeshBlock;
class Fluid;
class FluidIntegrator;
class ParameterInput;
struct InterfaceField;
struct NeighborBlock;

//-------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void ReflectInnerX1(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, InterfaceField &buf,
                    int is, int ie, int js, int je, int ks, int ke);

typedef void (*BValFluid_t)(MeshBlock *pmb, AthenaArray<Real> &buf,
                            int is, int ie, int js, int je, int ks, int ke);
typedef void (*BValField_t)(MeshBlock *pmb, InterfaceField &buf,
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

  void EnrollFluidBoundaryFunction (enum direction edge, BValFluid_t  my_bc);
  void EnrollFieldBoundaryFunction(enum direction edge, BValField_t my_bc);
  void CheckBoundary(void);

  int LoadFluidBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                       NeighborBlock& nb);
  int LoadFluidBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                       NeighborBlock& nb);
  int LoadFluidBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                     NeighborBlock& nb);
  void SendFluidBoundaryBuffers(AthenaArray<Real> &src, int step);
  void SetFluidBoundarySameLevel(AthenaArray<Real> &dst, Real *buf, NeighborBlock& nb);
  void SetFluidBoundaryFromCoarser(Real *buf, NeighborBlock& nb);
  void SetFluidBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf, NeighborBlock& nb);
  bool ReceiveFluidBoundaryBuffers(AthenaArray<Real> &dst, int step);
  void ReceiveFluidBoundaryBuffersWithWait(AthenaArray<Real> &dst, int step);
  void RestrictFluid(AthenaArray<Real> &src,
                     int si, int ei, int sj, int ej, int sk, int ek);
  void SendFluxCorrection(int step);
  bool ReceiveFluxCorrection(AthenaArray<Real> &dst, int step);
  void ProlongateFluidBoundaries(AthenaArray<Real> &dst);
  int LoadFieldBoundaryBufferSameLevel(InterfaceField &src, Real *buf,
                                       NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(InterfaceField &src, Real *buf,
                                       NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(InterfaceField &src, Real *buf,
                                     NeighborBlock& nb);
  void SendFieldBoundaryBuffers(InterfaceField &src, int step);
  void SetFieldBoundarySameLevel(InterfaceField &dst, Real *buf, NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(InterfaceField &dst, Real *buf, NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(InterfaceField &dst, int step);
  void ReceiveFieldBoundaryBuffersWithWait(InterfaceField &dst, int step);
  void RestrictFieldX1(AthenaArray<Real> &bx1f,
                       int si, int ei, int sj, int ej, int sk, int ek);
  void RestrictFieldX2(AthenaArray<Real> &bx2f,
                       int si, int ei, int sj, int ej, int sk, int ek);
  void RestrictFieldX3(AthenaArray<Real> &bx3f,
                       int si, int ei, int sj, int ej, int sk, int ek);
  void ProlongateFieldBoundaries(AthenaArray<Real> &dst);
  void SendEMFCorrection(int step);
  bool ReceiveEMFCorrection(int step);

  void FluidPhysicalBoundaries(AthenaArray<Real> &dst);
  void FieldPhysicalBoundaries(InterfaceField &dst);

  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);


private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  BValFluid_t FluidBoundary_[6];
  BValField_t FieldBoundary_[6];

  int nface_, nedge_;
  bool edge_flag_[12];

  enum boundary_status fluid_flag_[NSTEP][56], field_flag_[NSTEP][56];
  enum boundary_status flcor_flag_[NSTEP][6][2][2];
  enum boundary_status emfcor_fflag_[NSTEP][6][2][2];
  enum boundary_status emfcor_eflag_[NSTEP][12][2];
  Real *fluid_send_[NSTEP][56],   *fluid_recv_[NSTEP][56];
  Real *field_send_[NSTEP][56],   *field_recv_[NSTEP][56];
  Real *flcor_send_[NSTEP][6],    *flcor_recv_[NSTEP][6][2][2];
  Real *emfcor_fsend_[NSTEP][6],  *emfcor_frecv_[NSTEP][6][2][2];
  Real *emfcor_esend_[NSTEP][12], *emfcor_erecv_[NSTEP][12][2];
  AthenaArray<Real> coarse_cons_, coarse_prim_;
  AthenaArray<Real> fvol_[2][2], sarea_[2];
  AthenaArray<Real> surface_flux_[6];
  InterfaceField coarse_b_;

#ifdef MPI_PARALLEL
  MPI_Request req_fluid_send_[NSTEP][56],   req_fluid_recv_[NSTEP][56];
  MPI_Request req_field_send_[NSTEP][56],   req_field_recv_[NSTEP][56];
  MPI_Request req_flcor_send_[NSTEP][6],    req_flcor_recv_[NSTEP][6][2][2];
  MPI_Request req_emfcor_fsend_[NSTEP][6],  req_emfcor_frecv_[NSTEP][6][2][2];
  MPI_Request req_emfcor_fsend_[NSTEP][12], req_emfcor_frecv_[NSTEP][12][2];
#endif

  friend class FluidIntegrator;
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


#endif
