#ifndef BVALS_CC_BVALS_CC_HPP_
#define BVALS_CC_BVALS_CC_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.hpp
//  \brief handle boundaries for any AthenaArray type variable that represents a physical
//         quantity indexed along / located around cell-centers

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
//! \class CellCenteredBoundaryVariable
//  \brief

class CellCenteredBoundaryVariable : public BoundaryVariable {
 public:
  CellCenteredBoundaryVariable(MeshBlock *pmb,
                               AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
                               AthenaArray<Real> *var_flux);
  ~CellCenteredBoundaryVariable();

  AthenaArray<Real> *var_cc;
  AthenaArray<Real> coarse_buf;

  // currently, no need to ever switch flux[]. Keep as non-pointer member, shallow-copied
  AthenaArray<Real> x1flux, x2flux, x3flux;

  // maximum number of reserved unique "physics ID" component of MPI tag bitfield
  // (CellCenteredBoundaryVariable only actually uses 1x if multilevel==false, no shear)
  // must correspond to the # of "int *phys_id_" private members, below. Convert to array?
  static constexpr int max_phys_id = 3;

  // BoundaryVariable:
  int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) override;
  int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) override;

  // BoundaryCommunication:
  void SetupPersistentMPI() override;
  void StartReceiving(BoundaryCommSubset phase) override;
  void ClearBoundary(BoundaryCommSubset phase) override;
  void StartReceivingShear(BoundaryCommSubset phase) override;
  void ComputeShear(const Real time) override;

  // BoundaryBuffer:
  void SendFluxCorrection() override;
  bool ReceiveFluxCorrection() override;

  // Shearing box
  void SendShearingBoxBoundaryBuffersForInit();
  void SendShearingBoxBoundaryBuffers();
  bool ReceiveShearingBoxBoundaryBuffers();

  // BoundaryPhysics:
  void ReflectInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void ReflectOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void ReflectInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void ReflectOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void ReflectInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void ReflectOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;

  void OutflowInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void OutflowOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void OutflowInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void OutflowOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void OutflowInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void OutflowOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;

  void PolarWedgeInnerX2(Real time, Real dt,
                         int il, int iu, int jl, int kl, int ku, int ngh) override;
  void PolarWedgeOuterX2(Real time, Real dt,
                         int il, int iu, int ju, int kl, int ku, int ngh) override;

 protected:
  int nl_, nu_;
  const bool *flip_across_pole_;

 private:
  // BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;

  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;

  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;

  void PolarBoundarySingleAzimuthalBlock() override;

#ifdef MPI_PARALLEL
  int cc_phys_id_, cc_flx_phys_id_;
#endif

  // Shearing box
  // working arrays of remapped quantities
  AthenaArray<Real>  shear_cc_[2];
  // flux from conservative remapping
  AthenaArray<Real>  shear_flx_cc_[2];
  // KGF: these should probably be combined into a struct or array with send/recv switch
  int shear_send_count_cc_[2][4], shear_recv_count_cc_[2][4]; // buffer sizes

#ifdef MPI_PARALLEL
  int shear_cc_phys_id_;
#endif

  void LoadShearing(AthenaArray<Real> &src, Real *buf, int nb);
  void SetShearingBoxBoundarySameLevel(Real *buf, const int nb);
  // KGF: AthenaArray<Real>: shboxvar_inner/outer_hydro_, flx_inner/outer_hydro_
  void RemapFlux(const int n, const int k, const int jinner, const int jouter,
                 const int i, const Real eps, const AthenaArray<Real> &var,
                 AthenaArray<Real> &flux);
};

#endif // BVALS_CC_BVALS_CC_HPP_
