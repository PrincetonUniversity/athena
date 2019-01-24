#ifndef BVALS_FC_BVALS_FC_HPP_
#define BVALS_FC_BVALS_FC_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.hpp
//  \brief

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
  virtual ~FaceCenteredBoundaryVariable();

  FaceField &var_fc;

  // BoundaryCommunication pure virtual functions:
  virtual void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) override;
  virtual void DestroyBoundaryData(BoundaryData &bd) override;
  virtual void Initialize(void) override;
  virtual void StartReceivingForInit(bool cons_and_field) override;
  virtual void ClearBoundaryForInit(bool cons_and_field) override;
  virtual void StartReceivingAll(const Real time) override;
  virtual void ClearBoundaryAll(void) override;

  // BoundaryBuffer pure virtual functions:
  virtual int LoadBoundaryBufferSameLevel(
       int nl, int nu, Real *buf, const NeighborBlock& nb) override;
  // 1x LoadField*() don't use: int nl, int nu
  virtual void SendBoundaryBuffers( enum CCBoundaryType type) override;
  virtual bool ReceiveBoundaryBuffers(enum CCBoundaryType type) override;
  virtual void ReceiveAndSetBoundariesWithWait(
                                               enum CCBoundaryType type) override;
  virtual void SetBoundaries( enum CCBoundaryType type) override;
  // 4x Send/Receive/Set-FieldBoundaryBuffers() don't use: enum CCBoundaryType type
  virtual void SetBoundarySameLevel(
       int nl, int nu, Real *buf,
      const NeighborBlock& nb, bool *flip) override;
  virtual int LoadBoundaryBufferToCoarser(
       int nl, int nu, Real *buf, AthenaArray<Real> &cbuf,
      const NeighborBlock& nb) override;
  // cbuf parameter is unique to CC variable and the 2x Coarser load/set fns.
  // needed for switching HYDRO_CONS and HYDRO_PRIM
  virtual int LoadBoundaryBufferToFiner(
       int nl, int nu, Real *buf, const NeighborBlock& nb) override;
  virtual void SetBoundaryFromCoarser(
      int nl, int nu, Real *buf, AthenaArray<Real> &cbuf,
      const NeighborBlock& nb, bool *flip) override;
  virtual void SetBoundaryFromFiner(
      int nl, int nu,
      Real *buf, const NeighborBlock& nb, bool *flip) override;
  // 3x SetFieldBoundary*() don't use: int nl, int nu (like Load), but also "bool flip"
  virtual void SendFluxCorrection(enum FluxCorrectionType type) override;
  virtual bool ReceiveFluxCorrection(enum FluxCorrectionType type) override;
  // originally: SendEMFCorrection(void), ReceiveEMFCorrection(void)

  virtual void PolarBoundarySingleAzimuthalBlockField(FaceField &dst);

  // Face-centered/Field/EMF unique methods:
  void PolarBoundaryAverageField(FaceField &dst); // formerly PolarAxisFieldAverage()

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

  //-------------------- prototypes for all BC functions ---------------------------------
  virtual void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);

  virtual void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);
  virtual void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                              FaceField &b, int il, int iu, int jl, int ju,
                              int kl, int ku, int nu, int ngh);

  virtual void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 FaceField &b, int il, int iu, int jl,
                                 int ju, int kl, int ku, int nu, int ngh);
  virtual void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                                 FaceField &b, int il, int iu, int jl,
                                 int ju, int kl, int ku, int nu, int ngh);
  //protected:

 private:
  // standard Field and emf BV private variables
  BoundaryData bd_field_, bd_emfcor_;
  enum BoundaryStatus *emf_north_flag_;
  enum BoundaryStatus *emf_south_flag_;
  Real **emf_north_send_, **emf_north_recv_;
  Real **emf_south_send_, **emf_south_recv_;

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
