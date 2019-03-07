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
  CellCenteredBoundaryVariable(MeshBlock *pmb,
                               AthenaArray<Real> &var, AthenaArray<Real> *var_flux);
  ~CellCenteredBoundaryVariable();

  AthenaArray<Real> var_cc;
  // AthenaArray<Real> src, dst;

  // KGF: considered moving variable to derived HydroBoundaryVariable class
  AthenaArray<Real> coarse_buf;  // FaceCentered functions just use "pmr->coarse_b_.x1f"

  AthenaArray<Real> x1flux, x2flux, x3flux;

  // what about bool *flip? CC or Hydro-specific? Assuming CC now, see protected: vars

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

  // BoundaryVariable:
  int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) override;
  int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) override;

  // BoundaryCommunication:
  void Initialize() override;
  void StartReceivingForInit(bool cons_and_field) override;
  void ClearBoundaryForInit(bool cons_and_field) override;
  void StartReceivingAll(const Real time) override;
  void ClearBoundaryAll() override;

  // BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  void SendBoundaryBuffers() override;
  bool ReceiveBoundaryBuffers() override;
  void ReceiveAndSetBoundariesWithWait() override;
  void SetBoundaries() override;
  // "bool *flip" is passed in 3x Set...From*(), computed by enum switch in wrapper
  // function SetCellCenteredBoundaries(): nullptr (grav?) vs. flip_across_pole_hyd
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;
  // coarse_buf:
  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;
  // coarse_buf:
  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;

  void SendFluxCorrection() override;
  bool ReceiveFluxCorrection() override;
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
  // what about PolarBoundaryAverageField()? -- make analog no-ops for cell-centered var:
  // void PolarBoundaryAverage()
  // and for EMF:
  // void PolarBoundarySingleAzimuthalBlockFluxCorrection()
  // void PolarBoundaryAverageFluxCorrection()
  void PolarBoundarySingleAzimuthalBlock() override;

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

 protected:
  // CellCenteredBoundaryVariable is assumed 4D
  // (3D) nl=nu=0 for gravity
  // nl=0, nu=NHYDRO-1 for Hydro
  int nl_, nu_;
  const bool *flip_across_pole_;

 private:
  // standard cell-centered and flux BV private variables
  // KGF: currently declared in base BoundaryValues class
  //BoundaryData bd_cc_;
  //BoundaryData bd_cc_flcor_;

  // Pulling these variables out of function signatures, since FaceCentered
  // does not use them, only all CellCenteredBoundaryVariable instances (not specific to
  // Hydro, unlike enum CCBoundaryType and AthenaArray<Real> coarse_buf)

#ifdef MPI_PARALLEL
  int cc_phys_id_, cc_flx_phys_id_;
#endif

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
