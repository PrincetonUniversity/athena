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
  // Allow functions to access most any variable to allow for fully general BC dependence
  friend class Hydro;
  friend class Field;
 public:
  CellCenteredBoundaryVariable();
  ~CellCenteredBoundaryVariable();

  AthenaArray<Real> &var_cc;

  // Probably should pull these variables out of function signatures, since FaceCentered
  // does not use them, only all CellCenteredBoundaryVariable instances (not specific to
  // Hydro, unlike enum CCBoundaryType and AthenaArray<Real> coarse_buf
  // int nl, nu;

  // what about bool *flip? CC or Hydro-specific?

  // Current Facade class BoundaryValues calls in time_integrator.cpp:
  // StartReceivingAll(time);
  // ClearBoundaryAll();
  // SendFluxCorrection(FLUX_HYDRO);
  // SendEMFCorrection();
  // ReceiveFluxCorrection(FLUX_HYDRO)
  // ReceiveEMFCorrection()
  // SendCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  // SendFieldBoundaryBuffers(pmb->pfield->b);
  // ReceiveCellCenteredBoundaryBuffers(HYDRO_CONS);
  // ReceiveFieldBoundaryBuffers();
  // SetCellCenteredBoundaries(pmb->phydro->u, HYDRO_CONS);
  // SetFieldBoundaries(pmb->pfield->b);
  // +8x shearing box-specific functions

  // - Replace all of these pbval->fn() calls with phbval->fn() or pfbval->fn()
  // - Create a unique function for HydroBoundaryVariable to change:
  // enum CCBoundaryType, AthenaArray<Real> coarse_buf, &src, &dst

  // HYDRO_PRIM is passed only in 2x lines in mesh.cpp:
  // SendCellCenteredBoundaryBuffers(pmb->phydro->w, HYDRO_PRIM);
  // ReceiveAndSetCellCenteredBoundariesWithWait(pmb->phydro->w, HYDRO_PRIM);

  // BoundaryCommunication:
  void InitBoundaryData(BoundaryData &bd, enum BoundaryType type) override;
  void DestroyBoundaryData(BoundaryData &bd) override;
  void Initialize(void) override;
  void StartReceivingForInit(bool cons_and_field) override;
  void ClearBoundaryForInit(bool cons_and_field) override;
  void StartReceivingAll(const Real time) override;
  void ClearBoundaryAll(void) override;

  // BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(int nl, int nu, Real *buf,
                                  const NeighborBlock& nb) override;
  void SendBoundaryBuffers(void) override;
  bool ReceiveBoundaryBuffers(void) override;
  void ReceiveAndSetBoundariesWithWait(void) override;
  void SetBoundaries(void) override;
  void SetBoundarySameLevel(int nl, int nu, Real *buf,
                            const NeighborBlock& nb,
                            bool *flip) override;
  int LoadBoundaryBufferToCoarser(int nl, int nu, Real *buf,  // coarse_buf
                                  const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(int nl, int nu, Real *buf,
                                const NeighborBlock& nb) override;
  void SetBoundaryFromCoarser(int nl, int nu, Real *buf, // coarse_buf
                              const NeighborBlock& nb,
                              bool *flip) override;
  void SetBoundaryFromFiner(int nl, int nu, Real *buf,
                            const NeighborBlock& nb,
                            bool *flip) override;

  void SendFluxCorrection(enum FluxCorrectionType type) override;
  bool ReceiveFluxCorrection(enum FluxCorrectionType type) override;
  // TODO(felker): FLUX_HYDRO=0 is the only defined FluxCorrectionType enum in athena.hpp
  // TODO(felker): handle the 6x unique Field-related flux correction functions
  // Cell-centered flux correction functions are much simpler than Field counterpart
  // In addition to 2x simple Send/Recv EMFCorrection() functions, there are:
  // - 6x Load/Set EMF (not correction). No Load to finer, to Set to coarser, but
  //   Set/LoadEMFBoundaryPolarBuffer()
  // - AverageEMFBoundary(), ClearCoarseEMFBoundary(),
  //                         PolarBoundarySingleAzimuthalBlockEMF()

  // optional: compare to PolarBoundarySingleAzimuthalBlockField(),
  //                      PolarBoundarySingleAzimuthalBlockEMF()
  // what about PolarBoundaryAverageField()?
  void PolarBoundarySingleAzimuthalBlock(int nl, int nu);

  // Shearingbox Hydro
  // void LoadHydroShearing(AthenaArray<Real> &src, Real *buf, int nb);
  // void SendHydroShearingboxBoundaryBuffersForInit(AthenaArray<Real> &src, bool cons);
  // void SendHydroShearingboxBoundaryBuffers(AthenaArray<Real> &src, bool cons);

  // void SetHydroShearingboxBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
  //                                           const int nb);
  // bool ReceiveHydroShearingboxBoundaryBuffers(AthenaArray<Real> &dst);
  // void FindShearBlock(const Real time);
  // void RemapFlux(const int n, const int k, const int jinner, const int jouter,
  //                const int i, const Real eps, const AthenaArray<Real> &U,
  //                AthenaArray<Real> &Flux);

  // BoundaryPhysics:
  void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;

  void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;
  void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                      int il, int iu, int jl, int ju,
                      int kl, int ku, int nu, int ngh) override;

  void PolarWedgeInnerX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                         int il, int iu, int jl,
                         int ju, int kl, int ku, int nu, int ngh) override;
  void PolarWedgeOuterX2(MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
                         int il, int iu, int jl,
                         int ju, int kl, int ku, int nu, int ngh) override;
  //protected:

 private:
  // standard cell-centered and flux BV private variables
  BoundaryData bd_cc_, bd_cc_flcor_;

  // Shearingbox Hydro
  //   enum BoundaryStatus shbox_inner_hydro_flag_[4], shbox_outer_hydro_flag_[4];
  //   // working arrays of remapped quantities
  //   AthenaArray<Real>  shboxvar_inner_hydro_, shboxvar_outer_hydro_;
  //   // Hydro flux from conservative remapping
  //   AthenaArray<Real>  flx_inner_hydro_, flx_outer_hydro_;
  //   int  send_innersize_hydro_[4], recv_innersize_hydro_[4]; // buffer sizes
  //   Real *send_innerbuf_hydro_[4], *recv_innerbuf_hydro_[4]; // send and recv buffers
  //   int  send_outersize_hydro_[4], recv_outersize_hydro_[4]; // buffer sizes
  //   Real *send_outerbuf_hydro_[4], *recv_outerbuf_hydro_[4]; // send and recv buffers
  // #ifdef MPI_PARALLEL
  //   // MPI request for send and recv msgs
  //   MPI_Request rq_innersend_hydro_[4], rq_innerrecv_hydro_[4];
  //   MPI_Request rq_outersend_hydro_[4], rq_outerrecv_hydro_[4];
  // #endif
};

#endif // BVALS_CC_BVALS_CC_HPP_