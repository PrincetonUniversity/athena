#ifndef MESH_MESH_REFINEMENT_HPP_
#define MESH_MESH_REFINEMENT_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh_refinement.hpp
//  \brief defines MeshRefinement class used for static/adaptive mesh refinement

// C headers

// C++ headers
#include <tuple>
#include <vector>

// Athena++ headers
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
class FaceCenteredBoundaryVariable;
class HydroBoundaryVariable;

//----------------------------------------------------------------------------------------
//! \class MeshRefinement
//  \brief

class MeshRefinement {
  // KGF: needs to access pcoarsec in ProlongateBoundaries() for passing to BoundaryFunc()
  friend class BoundaryValues;

  // friend class FaceCenteredBoundaryVariable;
  // friend class HydroBoundaryVariable;
  // friend class MeshBlock;

  // KGF: needs to access refine_flag_ in Mesh::AdaptiveMeshRefinement(). Make var public?
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
  void CheckRefinementCondition();

  // setter functions for "enrolling" variable arrays in refinement via Mesh::AMR()
  // and/or in BoundaryValues::ProlongateBoundaries() (for SMR and AMR)
  void AddToRefinement(AthenaArray<Real> *pvar_cc, AthenaArray<Real> *pcoarse_cc);
  void AddToRefinement(FaceField *pvar_fc, FaceField *pcoarse_fc);

  // KGF: for switching first entry in pvars_cc_ to/from:
  // (w, coarse_prim) vs. (u, coarse_cons_)
  void SetHydroRefinement(HydroBoundaryQuantity hydro_type);

 private:
  // data
  MeshBlock *pmy_block_;
  Coordinates *pcoarsec;

  AthenaArray<Real> fvol_[2][2], sarea_x1_[2][2], sarea_x2_[2][3], sarea_x3_[3][2];
  int refine_flag_, neighbor_rflag_, deref_count_, deref_threshold_;

  // functions
  AMRFlagFunc AMRFlag_; // duplicate of member of Mesh class

  // KGF: tuples of references to AMR-enrolled (quantity, coarse_quantity) arrays
  // KGF: rename to something more descriptive, e.g. "refined_cc_vars"
  std::vector<std::tuple<AthenaArray<Real> *, AthenaArray<Real> *>> pvars_cc_;
  std::vector<std::tuple<FaceField *, FaceField *>> pvars_fc_;
};

#endif // MESH_MESH_REFINEMENT_HPP_
