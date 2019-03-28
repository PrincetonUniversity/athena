#ifndef BVALS_CC_HYDRO_BVALS_HYDRO_HPP_
#define BVALS_CC_HYDRO_BVALS_HYDRO_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_hydro.hpp
//  \brief

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \class CellCenteredBoundaryVariable
//  \brief

class HydroBoundaryVariable : public CellCenteredBoundaryVariable {
 public:
  HydroBoundaryVariable(MeshBlock *pmb,
                        AthenaArray<Real> *var_hydro, AthenaArray<Real> *coarse_var,
                        AthenaArray<Real> *var_flux,
                        HydroBoundaryQuantity hydro_type);
                                                // AthenaArray<Real> &prim);
  virtual ~HydroBoundaryVariable() = default;

  // switch between Hydro class members "u" and "w" (or "u" and "u1", ...)
  void SwapHydroQuantity(AthenaArray<Real> &var_hydro, HydroBoundaryQuantity hydro_type);
  void SelectCoarseBuffer(HydroBoundaryQuantity hydro_type);

  // BoundaryPhysics: need to flip sign of velocity vectors for Reflect*()
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
  //protected:
 private:
  // HydroBoundaryQuantity::prim is passed only in 2x lines in mesh.cpp:
  // SendCellCenteredBoundaryBuffers(pmb->phydro->w, HydroBoundaryQuantity::prim);
  // ReceiveAndSetCellCenteredBoundariesWithWait(pmb->phydro->w,
  //                                             HydroBoundaryQuantity::prim);

  // Hydro is a unique cell-centered variable because of the relationship between
  // HydroBoundaryQuantity::cons u and HydroBoundaryQuantity::prim w.
  HydroBoundaryQuantity hydro_type_;
};

#endif // BVALS_CC_HYDRO_BVALS_HYDRO_HPP_
