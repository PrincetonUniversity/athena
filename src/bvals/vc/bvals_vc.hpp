#ifndef BVALS_VC_BVALS_VC_HPP_
#define BVALS_VC_BVALS_VC_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_vc.hpp
//  \brief handle boundaries for any AthenaArray type variable that represents a physical
//         quantity indexed along / located around vertices

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
//! \class VertexCenteredBoundaryVariable
//  \brief

class VertexCenteredBoundaryVariable : public BoundaryVariable {
 public:
  VertexCenteredBoundaryVariable(MeshBlock *pmb,
                               AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
                                 AthenaArray<Real> *var_flux);
  ~VertexCenteredBoundaryVariable();

  // may want to rebind var_vc to u,u1,u2,w,w1, etc. registers for time integrator logic.
  // Also, derived class HydroBoundaryVariable needs to keep switching var and coarse_var
  // arrays between primitive and conserved variables ---> ptr members, not references
  AthenaArray<Real> *var_vc;
  AthenaArray<Real> *coarse_buf;  // may pass nullptr if mesh refinement is unsupported

  // currently, no need to ever switch flux[] ---> keep as reference members (not ptrs)
  // flux[3] w/ 3x empty AthenaArrays may be passed if mesh refinement is unsupported, but
  // nullptr is not allowed
  AthenaArray<Real> &x1flux, &x2flux, &x3flux;

  // maximum number of reserved unique "physics ID" component of MPI tag bitfield
  // (VertexCenteredBoundaryVariable only actually uses 1x if multilevel==false, no shear)
  // must correspond to the # of "int *phys_id_" private members, below. Convert to array?
  static constexpr int max_phys_id = 3;

  // BoundaryVariable:
  int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) override;
  // VC
  int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) {return 0;};

  // BoundaryCommunication:
  void SetupPersistentMPI() override;
  void StartReceiving(BoundaryCommSubset phase) override;
  void ClearBoundary(BoundaryCommSubset phase) override;
  // VC
  void StartReceivingShear(BoundaryCommSubset phase) {return;};
  void ComputeShear(const Real time) {return;};

  // VC
  // BoundaryBuffer:
  void ReceiveAndSetBoundariesWithWait() override;
  void SetBoundaries() override;
  void SendFluxCorrection() {return;};
  bool ReceiveFluxCorrection() {return false;};

  void ZeroVertexGhosts();
  void FinalizeVertexConsistency();

  // Internally increment Outflow / etc conditions to avoid interface changes
  // elsewhere
  inline int IncrementIfNonzero(int idx) {
    return (idx > 0) ? idx + 1 : idx;
  }

  // // Shearing box
  // void SendShearingBoxBoundaryBuffers();
  // bool ReceiveShearingBoxBoundaryBuffers();

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

  void ExtrapolateOutflowInnerX1(Real time, Real dt,
                                 int il, int jl, int ju, int kl, int ku,
                                 int ngh) override;
  void ExtrapolateOutflowOuterX1(Real time, Real dt,
                                 int iu, int jl, int ju, int kl, int ku,
                                 int ngh) override;
  void ExtrapolateOutflowInnerX2(Real time, Real dt,
                                 int il, int iu, int jl, int kl, int ku,
                                 int ngh) override;
  void ExtrapolateOutflowOuterX2(Real time, Real dt,
                                 int il, int iu, int ju, int kl, int ku,
                                 int ngh) override;
  void ExtrapolateOutflowInnerX3(Real time, Real dt,
                                 int il, int iu, int jl, int ju, int kl,
                                 int ngh) override;
  void ExtrapolateOutflowOuterX3(Real time, Real dt,
                                 int il, int iu, int jl, int ju, int ku,
                                 int ngh) override;

  void PolarWedgeInnerX2(Real time, Real dt,
                         int il, int iu, int jl, int kl, int ku, int ngh) override;
  void PolarWedgeOuterX2(Real time, Real dt,
                         int il, int iu, int ju, int kl, int ku, int ngh) override;

protected:
  int nl_, nu_;
  const bool *flip_across_pole_;

  // shearing box:
  // working arrays of remapped quantities
  // AthenaArray<Real>  shear_cc_[2];

private:
  // BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;

  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;

  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;

  void PolarBoundarySingleAzimuthalBlock() override;

  void ErrorUnknownMultiplicity();
  void ErrorIfPolarNotImplemented(const NeighborBlock& nb);
  void ErrorIfShearingBoxNotImplemented();

  void _FinalizeVertexConsistency3(int ox1, int ox2, int ox3);
  void _FinalizeVert3();
  void _FinalizeVert3a();

  void _FinalizeVert3noref();
  void _FinalizeVert2();
  void _FinalizeVert1();

#ifdef MPI_PARALLEL
  int vc_phys_id_; //, cc_flx_phys_id_;
#endif
// VC
//   // shearing box:
//   // flux from conservative remapping
//   AthenaArray<Real>  shear_flx_cc_[2];
//   // KGF: these should probably be combined into a struct or array with send/recv switch
//   int shear_send_count_cc_[2][4], shear_recv_count_cc_[2][4]; // buffer sizes

// #ifdef MPI_PARALLEL
//   int shear_cc_phys_id_;
// #endif

//   void LoadShearing(AthenaArray<Real> &src, Real *buf, int nb);
//   virtual void ShearQuantities(AthenaArray<Real> &shear_cc_, bool upper) {}
//   void SetShearingBoxBoundarySameLevel(Real *buf, const int nb);
//   // KGF: AthenaArray<Real>: shboxvar_inner/outer_hydro_, flx_inner/outer_hydro_
//   void RemapFlux(const int n, const int k, const int jinner, const int jouter,
//                  const int i, const Real eps, const AthenaArray<Real> &var,
//                  AthenaArray<Real> &flux);
};

#endif // BVALS_VC_BVALS_VC_HPP_
