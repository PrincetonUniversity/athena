#ifndef FIELD_FIELD_HPP_
#define FIELD_FIELD_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field.hpp
//! \brief defines Field class which implements data and functions for E/B fields

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/fc/bvals_fc.hpp"
#include "../coordinates/coordinates.hpp"
#include "field_diffusion/field_diffusion.hpp"

class MeshBlock;
class ParameterInput;
class Hydro;

//========================================================================================
//! \class Field
//! \brief electric and magnetic field data and functions

class Field {
  friend class Hydro;
 public:
  Field(MeshBlock *pmb, ParameterInput *pin);

  MeshBlock* pmy_block;  //!> ptr to MeshBlock containing this Field

  // face-centered magnetic fields
  FaceField b;       //!> time-integrator memory register #1
  FaceField b1;      //!> time-integrator memory register #2
  FaceField b2;      //!> time-integrator memory register #3
  FaceField b0, ct_update; // rkl2 STS memory registers
  // (no more than MAX_NREGISTER allowed)

  // cell-centered magnetic fields
  AthenaArray<Real> bcc;  //!> time-integrator memory register #1

  //-------- begin fourth-order MHD (UCT4)
  // In the 4th-order scheme, "b" and "bcc" hold face-AVERAGED and cell-AVERAGED fields;
  // the following hold the point-valued (face-/cell-centered) counterparts:
  FaceField b_fc;                 //!> fourth-order approx. to face-centered fields
  AthenaArray<Real> bcc_center;   //!> fourth-order approx. to cell-centered fields

  // fourth-order UCT reconstructions of transverse B and velocity at cell corners.
  // Naming (LDZ04): E/W = L1/R1, N/S = L2/R2 one-sided limits w.r.t. the corner
  AthenaArray<Real> by_W; // R1
  AthenaArray<Real> by_E; // L1
  AthenaArray<Real> bx_S; // R2
  AthenaArray<Real> bx_N; // L2
  // 3D states
  AthenaArray<Real> bz_R1, bz_L1;
  AthenaArray<Real> bz_R2, bz_L2;
  AthenaArray<Real> by_R3, by_L3;
  AthenaArray<Real> bx_R3, bx_L3;
  // 2D corner velocity states
  AthenaArray<Real> v_NE; // L2 L1
  AthenaArray<Real> v_SE; // R2 L1
  AthenaArray<Real> v_NW; // L2 R1
  AthenaArray<Real> v_SW; // R2 R1
  // 3D corner velocity states
  AthenaArray<Real> v_R3R2, v_R3L2, v_L3R2, v_L3L2;
  AthenaArray<Real> v_R3R1, v_R3L1, v_L3R1, v_L3L1;
  // max/min fast magnetosonic wavespeeds at interfaces (LDZ04 eq. 55)
  AthenaArray<Real> alpha_plus_x1_, alpha_minus_x1_;
  AthenaArray<Real> alpha_plus_x2_, alpha_minus_x2_;
  AthenaArray<Real> alpha_plus_x3_, alpha_minus_x3_;

  // 1D pencil scratch for the corner reconstruction sweeps
  AthenaArray<Real> by_W_, by_E_, bx_S_, bx_N_;
  AthenaArray<Real> bz_R1_, bz_L1_, bz_R2_, bz_L2_;
  AthenaArray<Real> by_R3_, by_L3_, bx_R3_, bx_L3_;
  AthenaArray<Real> v_NE_, v_SE_, v_NW_, v_SW_;
  AthenaArray<Real> v_R3R2_, v_R3L2_, v_L3R2_, v_L3L2_;
  AthenaArray<Real> v_R3R1_, v_R3L1_, v_L3R1_, v_L3L1_;
  AthenaArray<Real> vl_temp2_, vr_temp2_;
  AthenaArray<Real> v_NEb_, v_NWb_, bx_Nb_;
  //-------- end fourth-order MHD

  EdgeField e;    //!> edge-centered electric fields used in CT
  FaceField wght; //!> weights used to integrate E to corner using GS algorithm
  AthenaArray<Real> e2_x1f, e3_x1f; // electric fields at x1-face from Riemann solver
  AthenaArray<Real> e1_x2f, e3_x2f; // electric fields at x2-face from Riemann solver
  AthenaArray<Real> e1_x3f, e2_x3f; // electric fields at x3-face from Riemann solver

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_bcc_;
  int refinement_idx{-1};
  FaceField coarse_b_;

  FaceCenteredBoundaryVariable fbvar;
  FieldDiffusion fdif;

  void CalculateCellCenteredField(
      const FaceField &bf, AthenaArray<Real> &bc,
      Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);
  void CT(const Real wght, FaceField &b_out);
  void CT_STS(const Real wght, int stage, FaceField &b_out, FaceField &ct_update_out);
  void ComputeCornerE(AthenaArray<Real> &w, AthenaArray<Real> &bcc);
  void ComputeCornerE_STS();

  //------- begin fourth-order MHD functions
  void ComputeCornerE_UCT4();
  void CalculateCellCenteredFieldFourth(const FaceField &bf_center,
                                        AthenaArray<Real> &bc_center,
                                        Coordinates *pco, int il, int iu, int jl,
                                        int ju, int kl, int ku);
  void CellCenteredToAveragedField(const AthenaArray<Real> &bc_center,
                                   AthenaArray<Real> &bc, Coordinates *pco,
                                   int il, int iu, int jl, int ju, int kl, int ku);
  void FaceAveragedToCellAveragedField(const FaceField &bf, FaceField &bf_center,
                                       AthenaArray<Real> &bc,
                                       AthenaArray<Real> &bc_center, Coordinates *pco,
                                       int il, int iu, int jl, int ju, int kl, int ku);
  void CalculateFaceCenteredField(const FaceField &bf, FaceField &bf_center,
                                  Coordinates *pco, int il, int iu, int jl, int ju,
                                  int kl, int ku);
  //------- end fourth-order MHD functions

 private:
  // scratch space used to compute fluxes
  AthenaArray<Real> cc_e_;
  AthenaArray<Real> face_area_, edge_length_, edge_length_p1_;
  AthenaArray<Real> g_, gi_;  // only used in GR

  // fourth-order MHD: 4D scratch arrays
  AthenaArray<Real> scr1_nkji_cc_;
  AthenaArray<Real> scr1_kji_x1fc_, scr2_kji_x2fc_, scr3_kji_x3fc_;
};
#endif // FIELD_FIELD_HPP_
