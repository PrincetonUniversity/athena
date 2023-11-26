#ifndef BVALS_CC_BVALS_CC_HPP_
#define BVALS_CC_BVALS_CC_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone!@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.hpp
//! \brief handle boundaries for any AthenaArray type variable that represents a physical
//!        quantity indexed along / located around cell-centers

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
//! \brief

class CellCenteredBoundaryVariable : public BoundaryVariable {
 public:
  CellCenteredBoundaryVariable(MeshBlock *pmb,
                               AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
                               AthenaArray<Real> *var_flux, bool fflux);
    //override function for arrays need different initialization of nu_
  CellCenteredBoundaryVariable(MeshBlock *pmb,
                               AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
                               AthenaArray<Real> *var_flux, bool fflux, int flag);
  ~CellCenteredBoundaryVariable();

  //! \note
  //! may want to rebind var_cc to u,u1,u2,w,w1, etc. registers for time integrator logic.
  //! Also, derived class HydroBoundaryVariable needs to keep switching var and coarse_var
  //! arrays between primitive and conserved variables ---> ptr members, not references
  AthenaArray<Real> *var_cc;
  AthenaArray<Real> *coarse_buf;  //!< may pass nullptr if mesh refinement is unsupported

  //!@{
  //! \note
  //! currently, no need to ever switch flux[] ---> keep as reference members (not ptrs)
  //! flux[3] w/ 3x empty AthenaArrays may be passed if mesh refinement is unsupported,
  //! but nullptr is not allowed
  AthenaArray<Real> &x1flux, &x2flux, &x3flux;
  //!@}


  //! maximum number of reserved unique "physics ID" component of MPI tag bitfield
  //! \note
  //! (CellCenteredBoundaryVariable only actually uses 1x if multilevel==false, no shear)
  //! must correspond to the # of "int *phys_id_" private members, below.
  //! Convert to array?
  static constexpr int max_phys_id = 4;

  //!@{
  //! BoundaryVariable:
  int ComputeVariableBufferSize(const NeighborIndexes& ni, int cng) override;
  int ComputeFluxCorrectionBufferSize(const NeighborIndexes& ni, int cng) override;
  //!@}

  //!@{
  //! BoundaryCommunication:
  void SetupPersistentMPI() override;
  void StartReceiving(BoundaryCommSubset phase) override;
  void ClearBoundary(BoundaryCommSubset phase) override;
  void StartReceivingShear(BoundaryCommSubset phase) override;
  //!@}

  //!@{
  //! BoundaryBuffer:
  void SendFluxCorrection() override;
  bool ReceiveFluxCorrection() override;
  //!@}

  //!@{
  //! Shearing box
  void SendShearingBoxBoundaryBuffers();
  bool ReceiveShearingBoxBoundaryBuffers();
  void SetShearingBoxBoundaryBuffers();
  void SendFluxShearingBoxBoundaryBuffers();
  bool ReceiveFluxShearingBoxBoundaryBuffers();
  void SetFluxShearingBoxBoundaryBuffers();
  //!@}

  // for Multigrid
  void ExpandPhysicalBoundaries();

  //!@{
  //! BoundaryPhysics:
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

  void PolarWedgeInnerX2(Real time, Real dt,
                         int il, int iu, int jl, int kl, int ku, int ngh) override;
  void PolarWedgeOuterX2(Real time, Real dt,
                         int il, int iu, int ju, int kl, int ku, int ngh) override;
  //!@}

 protected:
  int nl_, nu_;
  const bool *flip_across_pole_;

  //! shearing box:
  //! working arrays of remapped quantities
  AthenaArray<Real>  shear_cc_[2];

 private:
  //!@{
  //! BoundaryBuffer:
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;

  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;
  //!@}

  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;

  virtual int LoadFluxBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadFluxBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);

  void SetFluxBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetFluxBoundaryFromFiner(Real *buf, const NeighborBlock& nb);

  void PolarBoundarySingleAzimuthalBlock() override;

#ifdef MPI_PARALLEL
  int cc_phys_id_, cc_flx_phys_id_;
#endif
  // shearing box:
  //! flux from conservative remapping
  AthenaArray<Real> pbuf;
  // KGF: these should probably be combined into a struct or array with send/recv switch
  int shear_send_count_cc_[2][4], shear_recv_count_cc_[2][4]; // buffer sizes

#ifdef MPI_PARALLEL
  int shear_cc_phys_id_;
#endif

  void LoadShearingBoxBoundarySameLevel(AthenaArray<Real> &src, Real *buf, int nb);
  virtual void ShearQuantities(AthenaArray<Real> &shear_cc_, bool upper) {}
  void SetShearingBoxBoundarySameLevel(AthenaArray<Real> &src,
                                       Real *buf, const int nb);

  AthenaArray<Real> shear_var_flx_[2];
  AthenaArray<Real> shear_map_flx_[2];
  int shear_send_count_flx_[2][3], shear_recv_count_flx_[2][3];
#ifdef MPI_PARALLEL
  int shear_flx_phys_id_;
#endif

  void LoadFluxShearingBoxBoundarySameLevel(AthenaArray<Real> &src,
                                            Real *buf, int nb);
  void SetFluxShearingBoxBoundarySameLevel(AthenaArray<Real> &src,
                                           Real *buf, const int nb);
  friend class ParticleMeshBoundaryVariable;
  friend class RadBoundaryVariable;
};

#endif // BVALS_CC_BVALS_CC_HPP_
