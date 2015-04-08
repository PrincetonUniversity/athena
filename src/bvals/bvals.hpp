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
//  int LoadFieldBoundaryBufferSameLevel(InterfaceField &src, Real *buf,
//                                       NeighborBlock& nb);
//  int LoadFieldBoundaryBufferToCoarser(InterfaceField &src, Real *buf,
//                                       NeighborBlock& nb);
//  int LoadFieldBoundaryBufferToFiner(InterfaceField &src, Real *buf,
//                                     NeighborBlock& nb);
  void SendFieldBoundaryBuffers(InterfaceField &src, int step);
//  void SetFieldBoundarySameLevel(InterfaceField &dst, Real *buf, NeighborBlock& nb);
//  void SetFieldBoundaryFromCoarser(Real *buf, NeighborBlock& nb);
//  void SetFieldBoundaryFromFiner(InterfaceField &dst, Real *buf, NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(InterfaceField &dst, int step);
  void ReceiveFieldBoundaryBuffersWithWait(InterfaceField &dst, int step);

  void FluidPhysicalBoundaries(AthenaArray<Real> &dst);
  void FieldPhysicalBoundaries(InterfaceField &dst);

  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);


private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  BValFluid_t FluidBoundary_[6];
  BValField_t FieldBoundary_[6];

  enum boundary_status fluid_flag_[NSTEP][56], field_flag_[NSTEP][56];
  Real *fluid_send_[NSTEP][56], *fluid_recv_[NSTEP][56];
  Real *field_send_[NSTEP][56], *field_recv_[NSTEP][56];

#ifdef MPI_PARALLEL
  MPI_Request req_fluid_send_[NSTEP][56], req_fluid_recv_[NSTEP][56];
  MPI_Request req_field_send_[NSTEP][56], req_field_recv_[NSTEP][56];
#endif
};

int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
int CreateMPITag(int lid, int flag, int phys, int ox1, int ox2, int ox3, int fi1, int fi2);
void CalculateTargetBufferID(int dim, bool multilevel, bool face_only);

#endif
