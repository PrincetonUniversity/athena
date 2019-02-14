#ifndef BVALS_BVALS_GRAV_HPP_
#define BVALS_BVALS_GRAV_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_grav.hpp
//  \brief defines GravityBoundaryValues class

// C headers

// C++ headers
#include <string>   // string

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../bvals.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class FFTBlock;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;

//! \struct GravityBoundaryData
//  \brief structure storing multigrid boundary information
struct GravityBoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56], sflag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
};

//----------------------------------------------------------------------------------------
//! \class GravityBoundaryValues
//  \brief BVals data and functions

class GravityBoundaryValues : public BoundaryBase {
 public:
  // compare with bvals.hpp:: BoundaryValues derived class
  // BoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs, ParameterInput *pin);
  GravityBoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs);
  ~GravityBoundaryValues();

  // Compare to bvals.hpp:
  // ----------------

  // new interface class BoundaryMemory functions:
  // --------
  // void InitBoundaryData(BoundaryData &bd, enum BoundaryType type);
  // void DestroyBoundaryData(BoundaryData &bd);
  void InitBoundaryData(GravityBoundaryData &bd);
  void DestroyBoundaryData(GravityBoundaryData &bd);

  // missing counterparts to:
  // void Initialize(void)
  // void StartReceivingForInit(bool cons_and_field)
  // void ClearBoundaryForInit(bool cons_and_field)

  // void StartReceivingAll(const Real time) final;
  void StartReceivingGravity(void);
  void ClearBoundaryGravity(void);

  // shared class BoundaryValues/BoundaryShared functions:
  // --------
  void ApplyPhysicalBoundaries(void);
  // missing counterparts to: ProlongateBoundaries, CheckBoundary

  // individual VariableBoundary classes: CellCenteredVariableBoundary,
  // FaceCenteredVariableBoundary
  // --------

  // minimum requirement of all 4x unrefined functions: load, send, set, recv
  int LoadGravityBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                         const NeighborBlock& nb);
  bool SendGravityBoundaryBuffers(AthenaArray<Real> &src);
  void SetGravityBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                   const NeighborBlock& nb);
  bool ReceiveGravityBoundaryBuffers(AthenaArray<Real> &dst);
  // No 2x Load*() + 2x Set*() SMR/AMR functions

  // No "Receive*BoundaryBuffersWithWait()" counterpart. perhaps optional, and must
  // manually be added to Mesh::Initialize() at this time. Should a function be required,
  // even if it does nothing?

  // -- gravity must do something in mesh initialization...

  // No 14x "physical Boundary Functions": 6x outflow + 6x reflect + 2x polar. Only
  // periodic BCs allowed for gravity

  // optional interfaces: RefinedBoundary, PolarBoundary, ...

 private:
  MeshBlock *pmy_block_;
  GravityBoundaryFunc GravityBoundaryFunction_[6];
  GravityBoundaryData bd_gravity_;
};

#endif // BVALS_BVALS_GRAV_HPP_
