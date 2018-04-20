#ifndef MESH_MESH_REFINEMENT_HPP_
#define MESH_MESH_REFINEMENT_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh_refinement.hpp
//  \brief defines MeshRefinement class used for static/adaptive mesh refinement

// Athena++ classes headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;
class ParameterInput;
class Coordinates;
struct FaceField;
class BoundaryValues;

//----------------------------------------------------------------------------------------
//! \class MeshRefinement
//  \brief

class MeshRefinement {
  friend class BoundaryValues;
  friend class MeshBlock;
  friend class Mesh;
public:
  MeshRefinement(MeshBlock *pmb, ParameterInput *pin);
  ~MeshRefinement();

  // functions
  void RestrictCellCenteredValues(const AthenaArray<Real> &fine,
                                  AthenaArray<Real> &coarse, int sn, int en,
                                  int csi, int cei, int csj, int cej, int csk, int cek);
  void RestrictFieldX1(const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
                       int csi, int cei, int csj, int cej, int csk, int cek);
  void RestrictFieldX2(const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
                       int csi, int cei, int csj, int cej, int csk, int cek);
  void RestrictFieldX3(const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
                       int csi, int cei, int csj, int cej, int csk, int cek);
  void ProlongateCellCenteredValues(const AthenaArray<Real> &coarse,
                                    AthenaArray<Real> &fine, int sn, int en,
                                    int si, int ei, int sj, int ej, int sk, int ek);
  void ProlongateSharedFieldX1(const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
                               int si, int ei, int sj, int ej, int sk, int ek);
  void ProlongateSharedFieldX2(const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
                               int si, int ei, int sj, int ej, int sk, int ek);
  void ProlongateSharedFieldX3(const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
                               int si, int ei, int sj, int ej, int sk, int ek);
  void ProlongateInternalField(FaceField &fine,
                               int si, int ei, int sj, int ej, int sk, int ek);
  void CheckRefinementCondition(void);

private:
  // data
  MeshBlock *pmy_block_;
  Coordinates *pcoarsec;
  AthenaArray<Real> coarse_cons_, coarse_prim_, coarse_bcc_;
  FaceField coarse_b_;
  AthenaArray<Real> fvol_[2][2], sarea_x1_[2][2], sarea_x2_[2][3], sarea_x3_[3][2];
  int refine_flag_, neighbor_rflag_, deref_count_, deref_threshold_;

  // functions
  AMRFlagFunc_t AMRFlag_;
};

#endif // MESH_MESH_REFINEMENT_HPP_
