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
  // Allow functions to access most any variable to allow for fully general BC dependence
  friend class Hydro;
  friend class Field;
 public:
  FaceCenteredBoundaryVariable();
  ~FaceCenteredBoundaryVariable();

  FaceField &var_fc;
  FaceField &src, &dst;

  // BoundaryCommunication:
  void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) override;
  void DestroyBoundaryData(BoundaryData &bd) override;
  void Initialize(void) override;
  void StartReceivingForInit(bool cons_and_field) override;
  void ClearBoundaryForInit(bool cons_and_field) override;
  void StartReceivingAll(const Real time) override;
  void ClearBoundaryAll(void) override;

  // BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf,
                                  const NeighborBlock& nb) override;
  // 1x LoadField*() don't use: int nl, int nu
  void SendBoundaryBuffers(void) override;
  bool ReceiveBoundaryBuffers(void) override;
  void ReceiveAndSetBoundariesWithWait(void) override;
  void SetBoundaries(void) override;
  // 4x Send/Receive/Set-FieldBoundaryBuffers() don't use: enum CCBoundaryType type
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  // cbuf parameter is unique to CC variable and the 2x Coarser load/set fns.
  // needed for switching HYDRO_CONS and HYDRO_PRIM
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;
  // 3x SetFieldBoundary*() don't use: int nl, int nu (like Load); also not "bool flip"
  void SendFluxCorrection(enum FluxCorrectionType type) override;
  bool ReceiveFluxCorrection(enum FluxCorrectionType type) override;
  // originally: SendEMFCorrection(void), ReceiveEMFCorrection(void)

  void PolarBoundarySingleAzimuthalBlock(void);

  // Face-centered/Field/EMF unique methods:
  void PolarBoundaryAverageField(); // formerly PolarAxisFieldAverage()

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb);

  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north);

  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  void PolarBoundarySingleAzimuthalBlockEMF(void);

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
  // void SendEMFShearingboxBoundaryCorrectionForInit(void);
  // void SendEMFShearingboxBoundaryCorrection(void);
  // void SetEMFShearingboxBoundarySameLevel(EdgeField &dst, Real *buf, const int nb);
  // bool ReceiveEMFShearingboxBoundaryCorrection(void);
  // void RemapEMFShearingboxBoundary(void);
  // void ClearEMFShearing(EdgeField &work);
  // void RemapFluxEMF(const int k, const int jinner, const int jouter, const Real eps,
  //                   const AthenaArray<Real> &U, AthenaArray<Real> &Flux);

  // BoundaryPhysics:
  void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);

  void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);
  void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh);

  void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                         int il, int iu, int jl,
                         int ju, int kl, int ku, int nu, int ngh);
  void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                         int il, int iu, int jl,
                         int ju, int kl, int ku, int nu, int ngh);
  //protected:

 private:
  // standard Field and emf BV private variables
  BoundaryData bd_fc_, bd_fc_flcor_; // bd_emfcor_;
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;

  // original bvals_fc.cpp functions never took "bool *flip" as a function parameter
  // because "flip_across_pole_field" was hardcoded in 3x SetFieldFrom*() fns
  bool *flip_across_pole;

#ifdef MPI_PARALLEL
  MPI_Request *req_emf_north_send_, *req_emf_north_recv_;
  MPI_Request *req_emf_south_send_, *req_emf_south_recv_;
#endif

  // Shearingbox Field
  //   enum BoundaryStatus shbox_inner_field_flag_[4], shbox_outer_field_flag_[4];
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
  //   enum BoundaryStatus shbox_inner_emf_flag_[5], shbox_outer_emf_flag_[5];
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
