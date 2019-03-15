#ifndef BVALS_FC_BVALS_FC_HPP_
#define BVALS_FC_BVALS_FC_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.hpp
//  \brief handle boundaries for any FaceField type variable that represents a physical
//         quantity indexed along / located around face-centers of cells

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
//! \class FaceCenteredBoundaryVariable
//  \brief

class FaceCenteredBoundaryVariable : public BoundaryVariable {
 public:
  FaceCenteredBoundaryVariable(MeshBlock *pmb, FaceField &var, EdgeField &var_flux);
  ~FaceCenteredBoundaryVariable();

  FaceField var_fc;
  // KGF: TODO mimic new CellCenteredBoundaryVariable member &coarse_buf (originally
  // passed as a function parameter as "cbuf"); Store reference to MeshRefinement
  // "pmr->coarse_b_"
  //FaceField &coarse_buf;

  // KGF: rename/change to single "EdgeField var_fc_flux"?
  AthenaArray<Real> e1;
  AthenaArray<Real> e2;
  AthenaArray<Real> e3;

  // BoundaryVariable:
  int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) override;
  int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) override;

  // BoundaryCommunication:
  void SetupPersistentMPI() override;
  // void StartReceivingForInit(bool cons_and_field) override;
  // void ClearBoundaryForInit(bool cons_and_field) override;
  void StartReceiving(BoundaryCommSubset phase) override;
  void ClearBoundary(BoundaryCommSubset phase) override;

  // BoundaryBuffer:
  // 1x LoadField*() don't use: int nl, int nu
  void SendBoundaryBuffers() override;
  bool ReceiveBoundaryBuffers() override;
  void ReceiveAndSetBoundariesWithWait() override;
  void SetBoundaries() override;
  void SendFluxCorrection() override;
  bool ReceiveFluxCorrection() override;
  // originally: SendEMFCorrection(), ReceiveEMFCorrection()

  // Shearingbox Field
  // void LoadFieldShearing(FaceField &src, Real *buf, int nb);
  // void SendFieldShearingboxBoundaryBuffersForInit(FaceField &src, bool cons);
  // void SendFieldShearingboxBoundaryBuffers(FaceField &src, bool cons);
  // void SetFieldShearingboxBoundarySameLevel(FaceField &dst, Real *buf, const int nb);
  // bool ReceiveFieldShearingboxBoundaryBuffers(FaceField &dst);
  // void RemapFluxField(const int k, const int jinner, const int jouter, const int i,
  //                     const Real eps, const AthenaArray<Real> &U,
  //                     AthenaArray<Real> &Flux);
  // // Shearingbox EMF
  // void LoadEMFShearing(EdgeField &src, Real *buf, const int nb);
  // void SendEMFShearingboxBoundaryCorrectionForInit();
  // void SendEMFShearingboxBoundaryCorrection();
  // void SetEMFShearingboxBoundarySameLevel(EdgeField &dst, Real *buf, const int nb);
  // bool ReceiveEMFShearingboxBoundaryCorrection();
  // void RemapEMFShearingboxBoundary();
  // void ClearEMFShearing(EdgeField &work);
  // void RemapFluxEMF(const int k, const int jinner, const int jouter, const Real eps,
  //                   const AthenaArray<Real> &U, AthenaArray<Real> &Flux);

  // BoundaryPhysics:
  void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;

  void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;
  void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int ngh) override;

  void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                         int il, int iu, int jl,
                         int ju, int kl, int ku, int ngh) override;
  void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                         int il, int iu, int jl,
                         int ju, int kl, int ku, int ngh) override;
  //protected:

 private:
  // standard Field and emf BV private variables
  // KGF: currently declared in base BoundaryValues class:
  //BoundaryData bd_fc_;
  //BoundaryData bd_fc_flcor_; // bd_emfcor_;
  BoundaryStatus *flux_north_flag_;
  BoundaryStatus *flux_south_flag_;
  Real **flux_north_send_, **flux_north_recv_;
  Real **flux_south_send_, **flux_south_recv_;

  // original bvals_fc.cpp functions never took "bool *flip" as a function parameter
  // because "flip_across_pole_field" was hardcoded in 3x SetFieldFrom*() fns
  const bool *flip_across_pole_;

  bool edge_flag_[12];
  int nedge_fine_[12];

  // ready to recv flux from same level and apply correction? false= 2nd pass for fine lvl
  bool recv_flx_same_lvl_;
  // KGF: formerly "firsttime_". The variable switch is used in only 2x functions:
  // ReceiveEMFCorrection() and StartReceivingAll()

#ifdef MPI_PARALLEL
  int fc_phys_id_, fc_flx_phys_id_, fc_flx_pole_phys_id_;
  MPI_Request *req_flux_north_send_, *req_flux_north_recv_;
  MPI_Request *req_flux_south_send_, *req_flux_south_recv_;
#endif

  // BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  // 4x Send/Receive/Set-FieldBoundaryBuffers() don't use: HydroBoundaryQuantity type
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  // cbuf parameter is unique to CC variable and the 2x Coarser load/set fns.
  // needed for switching HydroBoundaryQuantity::cons and HydroBoundaryQuantity::prim
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;
  // 3x SetFieldBoundary*() don't use: int nl, int nu (like Load); also not "bool flip"
  void PolarBoundarySingleAzimuthalBlock() override;

  // Face-centered/Field/EMF unique class methods:
  // KGF: should we require CellCenteredBoundaryVariable class to provide similar
  // encapsulated functions (or no-op fns) for load/set fluxes?

  // called in SetBoundaries() and ReceiveAndSetBoundariesWithWait()
  void PolarFieldBoundaryAverage(); // formerly PolarAxisFieldAverage()

  // all 3x only called in SendFluxCorrection()
  int LoadFluxBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadFluxBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadFluxBoundaryBufferToPolar(Real *buf, const PolarNeighborBlock &nb);

  // all 6x only called in ReceiveFluxCorrection()
  void SetFluxBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetFluxBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetFluxBoundaryFromPolar(Real **buf_list, int num_bufs, bool north);

  void ClearCoarseFluxBoundary();
  void AverageFluxBoundary();  // average flux from fine and equal lvls
  void PolarFluxBoundarySingleAzimuthalBlock();
  // called in SetupPersistentMPI()
  void CountFineEdges();

  void CopyPolarBufferSameProcess(const PolarNeighborBlock& nb, int ssize,
                                  int polar_block_index, bool is_north);
  // Shearingbox Field
  //   BoundaryStatus shbox_inner_field_flag_[4], shbox_outer_field_flag_[4];
  //   FaceField shboxvar_inner_field_, shboxvar_outer_field_;
  //   FaceField flx_inner_field_, flx_outer_field_;
  //   int  send_innersize_field_[4], recv_innersize_field_[4];
  //   Real *send_innerbuf_field_[4], *recv_innerbuf_field_[4];
  //   int  send_outersize_field_[4], recv_outersize_field_[4];
  //   Real *send_outerbuf_field_[4], *recv_outerbuf_field_[4];
  // #ifdef MPI_PARALLEL
  //   MPI_Request rq_innersend_field_[4], rq_innerrecv_field_[4];
  //   MPI_Request rq_outersend_field_[4], rq_outerrecv_field_[4];
  // #endif
  //   // Shearing box EMF correction
  //   BoundaryStatus shbox_inner_emf_flag_[5], shbox_outer_emf_flag_[5];
  //   EdgeField shboxvar_inner_emf_, shboxvar_outer_emf_;
  //   EdgeField shboxmap_inner_emf_, shboxmap_outer_emf_;
  //   EdgeField flx_inner_emf_, flx_outer_emf_;
  //   int  send_innersize_emf_[4], recv_innersize_emf_[4];
  //   Real *send_innerbuf_emf_[4], *recv_innerbuf_emf_[4];
  //   int  send_outersize_emf_[4], recv_outersize_emf_[4];
  //   Real *send_outerbuf_emf_[4], *recv_outerbuf_emf_[4];
  // #ifdef MPI_PARALLEL
  //   MPI_Request rq_innersend_emf_[4],  rq_innerrecv_emf_[4];
  //   MPI_Request rq_outersend_emf_[4],  rq_outerrecv_emf_[4];
  // #endif
};

#endif // BVALS_FC_BVALS_FC_HPP_
