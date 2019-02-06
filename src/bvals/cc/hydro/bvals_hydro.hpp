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

// Athena++ classes headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../bvals.hpp"
#include "../../bvals_interfaces.hpp"
#include "../bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \class CellCenteredBoundaryVariable
//  \brief

class HydroBoundaryVariable : public CellCenteredBoundaryVariable {
 public:
  HydroBoundaryVariable(MeshBlock *pmb, enum BoundaryType type,
                        AthenaArray<Real> &cons, enum HydroBoundaryType hydro_type);
                                                // AthenaArray<Real> &prim);
  ~HydroBoundaryVariable();

  // switch between Hydro class members "u" and "w"
  void SelectCoarseBuffer(enum HydroBoundaryType hydro_type);

  // BoundaryPhysics: need to flip sign of velocity vectors for Reflect*()
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
  //protected:
 private:
  // HYDRO_PRIM is passed only in 2x lines in mesh.cpp:
  // SendCellCenteredBoundaryBuffers(pmb->phydro->w, HYDRO_PRIM);
  // ReceiveAndSetCellCenteredBoundariesWithWait(pmb->phydro->w, HYDRO_PRIM);

  // Hydro is a unique cell-centered variable because of the relationship between
  // HYDRO_CONS u and HYDRO_PRIM w.
  enum HydroBoundaryType hydro_type_;
};

#endif // BVALS_CC_HYDRO_BVALS_HYDRO_HPP_
