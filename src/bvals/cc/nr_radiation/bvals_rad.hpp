#ifndef BVALS_CC_NR_RADIATION_BVALS_RAD_HPP_
#define BVALS_CC_NR_RADIATION_BVALS_RAD_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_rad.hpp
//  \brief
//========================================================================================
// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "../bvals_cc.hpp"


//----------------------------------------------------------------------------------------
//! \class CellCenteredBoundaryVariable
//  \brief

class RadBoundaryVariable : public CellCenteredBoundaryVariable {
 public:
  RadBoundaryVariable(MeshBlock *pmb,
                      AthenaArray<Real> *var_rad, AthenaArray<Real> *coarse_var,
                      AthenaArray<Real> *var_flux);
  virtual ~RadBoundaryVariable() = default;

  // functions unique implementation to radiation class
  void SendFluxCorrection() override;
  bool ReceiveFluxCorrection() override;


  void SetBoundaries() override;


  // function for shearing box
  void AddRadShearForInit();
  void ShearQuantities(AthenaArray<Real> &shear_cc_, bool upper) override;
  void SetShearingBoxBoundaryBuffers();
  void SetFluxShearingBoxBoundaryBuffers();
  void SendFluxShearingBoxBoundaryBuffers();
  bool ReceiveFluxShearingBoxBoundaryBuffers();
  bool ReceiveShearingBoxBoundaryBuffers();
  void SendShearingBoxBoundaryBuffers();

  // BoundaryPhysics: need to rotate the intensity
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

  void VacuumInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void VacuumOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void VacuumInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void VacuumOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void VacuumInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void VacuumOuterX3(Real time, Real dt,
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

  // rotate the specific intensity by some amount
  // these functions are typically called with periodic boundary together
  void RotateHPi_InnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh);

  void RotateHPi_OuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh);

  void RotateHPi_InnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh);

  void RotateHPi_OuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh);

  void RotatePi_InnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh);

  void RotatePi_OuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh);

  void PolarWedgeInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;

  void PolarWedgeOuterX2(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;

 private:
  AthenaArray<Real>  azimuthal_shift_rad_;

  // override function for flux correction
  void SetFluxBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetFluxBoundarySameLevel(Real *buf, const NeighborBlock& nb);

  int LoadFluxBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadFluxBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;

  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;

  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;

  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;


  void PolarBoundarySingleAzimuthalBlock() override;

  // shearing box functions
  void LoadShearingBoxBoundarySameLevel(AthenaArray<Real> &src, Real *buf, int nb);
  void SetShearingBoxBoundarySameLevel(AthenaArray<Real> &src,
                                       Real *buf, const int nb);
  void LoadFluxShearingBoxBoundarySameLevel(AthenaArray<Real> &src,
                                            Real *buf, int nb);
  void SetFluxShearingBoxBoundarySameLevel(AthenaArray<Real> &src,
                                           Real *buf, const int nb);

  AthenaArray<Real> ir_cm_, ir_lab_, pflux_; // co-moving frame specific inteisites
};

#endif //  BVALS_CC_NR_RADIATION_BVALS_RAD_HPP_
