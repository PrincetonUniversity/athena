#ifndef BVALS_SIXRAY_BVALS_SIXRAY_HPP_
#define BVALS_SIXRAY_BVALS_SIXRAY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone!@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_sixray.hpp
//! \brief handle boundaries for six-ray column density

// C headers

// C++ headers

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../bvals.hpp"
#include "../bvals_interfaces.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \class SixRayBoundaryVariable
//! \brief Class for six-ray boundary

class SixRayBoundaryVariable : public BoundaryVariable {
 public:
  // only works for uniform mesh for now
  SixRayBoundaryVariable(MeshBlock *pmb, AthenaArray<Real> *var);
  ~SixRayBoundaryVariable();

  AthenaArray<Real> *var;

  //! maximum number of reserved unique "physics ID" component of MPI tag bitfield
  static constexpr int max_phys_id = 1;

  //!@{
  //! BoundaryVariable:
  int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) override;
  int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni,
                                      int cng) override {return 0;};
  //!@}

  //!@{
  //! BoundaryCommunication:
  void SetupPersistentMPI() override;
  void StartReceiving(BoundaryCommSubset phase) override;
  void ClearBoundary(BoundaryCommSubset phase) override;
  void StartReceivingShear(BoundaryCommSubset phase) override {return;};
  //!@}

  //!@{
  //! BoundaryBuffer:
  void SendFluxCorrection() override {return;};
  bool ReceiveFluxCorrection() override {return true;};
  //!@}

  //!@{
  //! BoundaryPhysics:
  void ReflectInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override {return;};
  void ReflectOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override {return;};
  void ReflectInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override {return;};
  void ReflectOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override {return;};
  void ReflectInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override {return;};
  void ReflectOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override {return;};

  void OutflowInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override {return;};
  void OutflowOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override {return;};
  void OutflowInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override {return;};
  void OutflowOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override {return;};
  void OutflowInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override {return;};
  void OutflowOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override {return;};

  void VacuumInnerX1(Real time, Real dt,
                     int il, int jl, int ju, int kl, int ku, int ngh) override {return;};
  void VacuumOuterX1(Real time, Real dt,
                     int iu, int jl, int ju, int kl, int ku, int ngh) override {return;};
  void VacuumInnerX2(Real time, Real dt,
                     int il, int iu, int jl, int kl, int ku, int ngh) override {return;};
  void VacuumOuterX2(Real time, Real dt,
                     int il, int iu, int ju, int kl, int ku, int ngh) override {return;};
  void VacuumInnerX3(Real time, Real dt,
                     int il, int iu, int jl, int ju, int kl, int ngh) override {return;};
  void VacuumOuterX3(Real time, Real dt,
                     int il, int iu, int jl, int ju, int ku, int ngh) override {return;};

  void PolarWedgeInnerX2(Real time, Real dt, int il, int iu, int jl,
                         int kl, int ku, int ngh) override {return;};
  void PolarWedgeOuterX2(Real time, Real dt, int il, int iu, int ju,
                         int kl, int ku, int ngh) override {return;};
  //!@}

  // send to specific direction
  void SendSixRayBoundaryBuffers(const BoundaryFace direction);
  // receive from specific direction
  bool ReceiveAndSetSixRayBoundaryBuffers(const BoundaryFace direction);
  // get opposite direction of face boundary
  BoundaryFace GetOppositeBoundaryFace(const BoundaryFace direction);
  // get face neighbor
  NeighborBlock *GetFaceNeighbor(const BoundaryFace direction);

 private:
  int mu_, ml_;  // indexing of the first dimension (chemical species)
#ifdef MPI_PARALLEL
  int sixray_phys_id_;
#endif

  //!@{
  //! BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;

  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;

  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;
  void PolarBoundarySingleAzimuthalBlock() override {return;};
  //!@}
};

#endif // BVALS_SIXRAY_BVALS_SIXRAY_HPP_
