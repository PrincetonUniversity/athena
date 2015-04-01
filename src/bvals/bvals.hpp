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

//-------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, AthenaArray<Real> &buf);
void ReflectInnerX2(MeshBlock *pmb, AthenaArray<Real> &buf);
void ReflectInnerX3(MeshBlock *pmb, AthenaArray<Real> &buf);
void ReflectOuterX1(MeshBlock *pmb, AthenaArray<Real> &buf);
void ReflectOuterX2(MeshBlock *pmb, AthenaArray<Real> &buf);
void ReflectOuterX3(MeshBlock *pmb, AthenaArray<Real> &buf);

void ReflectInnerX1(MeshBlock *pmb, InterfaceField &buf);
void ReflectInnerX2(MeshBlock *pmb, InterfaceField &buf);
void ReflectInnerX3(MeshBlock *pmb, InterfaceField &buf);
void ReflectOuterX1(MeshBlock *pmb, InterfaceField &buf);
void ReflectOuterX2(MeshBlock *pmb, InterfaceField &buf);
void ReflectOuterX3(MeshBlock *pmb, InterfaceField &buf);

void OutflowInnerX1(MeshBlock *pmb, AthenaArray<Real> &buf);
void OutflowInnerX2(MeshBlock *pmb, AthenaArray<Real> &buf);
void OutflowInnerX3(MeshBlock *pmb, AthenaArray<Real> &buf);
void OutflowOuterX1(MeshBlock *pmb, AthenaArray<Real> &buf);
void OutflowOuterX2(MeshBlock *pmb, AthenaArray<Real> &buf);
void OutflowOuterX3(MeshBlock *pmb, AthenaArray<Real> &buf);

void OutflowInnerX1(MeshBlock *pmb, InterfaceField &buf);
void OutflowInnerX2(MeshBlock *pmb, InterfaceField &buf);
void OutflowInnerX3(MeshBlock *pmb, InterfaceField &buf);
void OutflowOuterX1(MeshBlock *pmb, InterfaceField &buf);
void OutflowOuterX2(MeshBlock *pmb, InterfaceField &buf);
void OutflowOuterX3(MeshBlock *pmb, InterfaceField &buf);

typedef void (*BValFluid_t)(MeshBlock *pmb, AthenaArray<Real> &buf);
typedef void (*BValField_t)(MeshBlock *pmb, InterfaceField &buf);

void InitBoundaryBuffer(int nx1, int nx2, int nx3);

//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  void Initialize(void);
  void StartReceivingForInit(void);
  void StartReceivingAll(void);

  void LoadAndSendFluidBoundaryBuffer(enum direction dir,
                                      AthenaArray<Real> &src, int flag);
  bool ReceiveAndSetFluidBoundary(enum direction dir, AthenaArray<Real> &dst, int flag);
  bool ReceiveAndSetFluidBoundaryWithWait(enum direction dir,
                                          AthenaArray<Real> &dst, int flag);
  void LoadAndSendFieldBoundaryBuffer(enum direction dir,
                                      InterfaceField &src, int flag);
  bool ReceiveAndSetFieldBoundary(enum direction dir, InterfaceField &dst, int flag);
  bool ReceiveAndSetFieldBoundaryWithWait(enum direction dir,
                                          InterfaceField &dst, int flag);

  void EnrollFluidBoundaryFunction (enum direction edge, BValFluid_t  my_bc);
  void EnrollFieldBoundaryFunction(enum direction edge, BValField_t my_bc);

  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);

  void CheckBoundary(void);

private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  BValFluid_t FluidBoundary_[6];
  BValField_t FieldBoundary_[6];

  char fluid_flag_[NSTEP][6][2][2];
  char field_flag_[NSTEP][6][2][2];
  Real *fluid_send_[NSTEP][6], *fluid_recv_[NSTEP][6];
  Real *field_send_[NSTEP][6], *field_recv_[NSTEP][6];

#ifdef MPI_PARALLEL
  MPI_Request req_fluid_send_[NSTEP][6][2][2], req_fluid_recv_[NSTEP][6][2][2];
  MPI_Request req_field_send_[NSTEP][6][2][2], req_field_recv_[NSTEP][6][2][2];
#endif
};

#ifdef MPI_PARALLEL
inline int CreateMPITag(int lid, int flag, int dir, int phys, int fb1 = 0, int fb2 = 0)
{
// tag = local id of destination (20) + flag (2) + dir (3) + face (2) + physics (4)
  return (lid<<11) | (flag<<9) | (dir<<6) | (fb1<<5) | (fb2<<4) | phys;
}
#endif

#endif
