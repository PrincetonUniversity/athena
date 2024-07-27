#ifndef COORDINATES_COORDINATES_HPP_
#define COORDINATES_COORDINATES_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file coordinates.hpp
//! \brief defines abstract base and derived classes for coordinates.
//!
//! This classe provides data and functions to compute/store coordinate positions and
//! spacing, as well as geometrical factors (areas, volumes, coordinate source terms)
//! for various coordinate systems.
//!
//! Change in version 2024: now we have only Coordinates class and no derived class.
//! One of the coordinate source files is selected by the configuration script,
//! and default functions are implemented with "__attribute__((weak))".

// C headers

// C++ headers
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"

// forward declarations
class MeshBlock;
class ParameterInput;

//----------------------------------------------------------------------------------------
//! \class Coordinates
//! \brief abstract base class for all coordinate derived classes

class Coordinates {
 public:
  friend class HydroSourceTerms;
  Coordinates(MeshBlock *pmb, ParameterInput *pin, bool flag = false);
  ~Coordinates();

  void Initialize(ParameterInput *pin);

  // data
  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Coordinates
  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f;    // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v;    // volume spacing and positions
  AthenaArray<Real> x1s2, x1s3, x2s1, x2s3, x3s1, x3s2; // area averaged positions for AMR
  // geometry coefficients (only used in SphericalPolar, Cylindrical, Cartesian)
  AthenaArray<Real> h2f, dh2fd1, h31f, h32f, dh31fd1, dh32fd2;
  AthenaArray<Real> h2v, dh2vd1, h31v, h32v, dh31vd1, dh32vd2;

  // functions...
  // ...to compute length of edges
  void Edge1Length(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &len);
  void Edge2Length(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &len);
  void Edge3Length(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &len);
  Real GetEdge1Length(const int k, const int j, const int i);
  Real GetEdge2Length(const int k, const int j, const int i);
  Real GetEdge3Length(const int k, const int j, const int i);
  // ...to compute length connecting cell centers (for non-ideal MHD)
  void VolCenter1Length(const int k, const int j, const int il, const int iu,
                        AthenaArray<Real> &len);
  void VolCenter2Length(const int k, const int j, const int il, const int iu,
                        AthenaArray<Real> &len);
  void VolCenter3Length(const int k, const int j, const int il, const int iu,
                        AthenaArray<Real> &len);
  // ...to compute physical width at cell center
  void CenterWidth1(const int k, const int j, const int il, const int iu,
                    AthenaArray<Real> &dx1);
  void CenterWidth2(const int k, const int j, const int il, const int iu,
                    AthenaArray<Real> &dx2);
  void CenterWidth3(const int k, const int j, const int il, const int iu,
                    AthenaArray<Real> &dx3);

  // ...to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area);
  void Face2Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area);
  Real GetFace1Area(const int k, const int j, const int i);
  Real GetFace2Area(const int k, const int j, const int i);
  Real GetFace3Area(const int k, const int j, const int i);
  // ...to compute area of faces joined by cell centers (for non-ideal MHD)
  void VolCenterFace1Area(const int k, const int j, const int il, const int iu,
                                  AthenaArray<Real> &area);
  void VolCenterFace2Area(const int k, const int j, const int il, const int iu,
                                  AthenaArray<Real> &area);
  void VolCenterFace3Area(const int k, const int j, const int il, const int iu,
                                  AthenaArray<Real> &area);

  // ...to compute Laplacian of quantities in the coord system and orthogonal subspaces
  void Laplacian(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                 const int il, const int iu, const int jl, const int ju,
                 const int kl, const int ku, const int nl, const int nu);
  void LaplacianX1(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                   const int n, const int k, const int j, const int il, const int iu);
  void LaplacianX1All(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                      const int nl, const int nu, const int kl, const int ku,
                      const int jl, const int ju, const int il, const int iu);
  void LaplacianX2(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                   const int n, const int k, const int j, const int il, const int iu);
  void LaplacianX2All(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                      const int nl, const int nu, const int kl, const int ku,
                      const int jl, const int ju, const int il, const int iu);
  void LaplacianX3(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                   const int n, const int k, const int j, const int il, const int iu);
  void LaplacianX3All(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                      const int nl, const int nu, const int kl, const int ku,
                      const int jl, const int ju, const int il, const int iu);

  // ...to compute volume of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
                  AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  void AddCoordTermsDivergence(const Real dt, const AthenaArray<Real> *flux,
                             const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
                             AthenaArray<Real> &u);

  // coordinate source term for cosmic rays
  void AddCRCoordTermsDivergence(const AthenaArray<Real> &u_cr,
                                 AthenaArray<Real> &coord_src);

  // Coordinate (Geometric) source term for Grad Pc
  void CRGradPcCoordTermsDivergence(const AthenaArray<Real> &u_cr,
                                    AthenaArray<Real> &grad_pc);

  void AddCoordTermsDivergence_STS(const Real dt, int stage,
                                   const AthenaArray<Real> *flux,
                                   AthenaArray<Real> &u, AthenaArray<Real> &flux_div);


    // ...define the direction of radiation angle with respect to the local axis
  void AxisDirection(int *axisx, int *axisy, int *axisz);
  void ConvertAngle(MeshBlock *pmb, const int nang, AthenaArray<Real> &mu);

// Functions for angular flux
  void ZetaArea(NRRadiation *prad, AthenaArray<Real> &area);
  void PsiArea(NRRadiation *prad, AthenaArray<Real> &area);
  void AngularVol(NRRadiation *prad, AthenaArray<Real> &vol);
  void GetGeometryZeta(NRRadiation *prad, const int k, const int j,
                       const int i, AthenaArray<Real> &g_zeta);
  void GetGeometryPsi(NRRadiation *prad, const int k, const int j,
                      const int i, const int n_zeta, AthenaArray<Real> &g_psi);
  // function overwirte in case nzeta = 0
  void GetGeometryPsi(NRRadiation *prad, const int k, const int j,
                      const int i, AthenaArray<Real> &g_psi);

  // ...to determine if index is a pole
  bool IsPole(int j);


  // In GR, functions...
  // ...to return private variables
  Real GetMass() const {return bh_mass_;}
  Real GetSpin() const {return bh_spin_;}

  // ...to compute metric
  void Metric(Real x1, Real x2, Real x3, ParameterInput *pin, AthenaArray<Real> &g,
       AthenaArray<Real> &g_inv, AthenaArray<Real> &dg_dx1, AthenaArray<Real> &dg_dx2,
       AthenaArray<Real> &dg_dx3);
  void CellMetric(const int k, const int j, const int il, const int iu,
                  AthenaArray<Real> &g, AthenaArray<Real> &gi);
  void Face1Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  // ...to transform primitives to locally flat space
  void PrimToLocal1(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
                    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal2(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
                    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal3(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
                    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);

  // ...to transform fluxes in locally flat space to global frame
  void FluxToGlobal1(const int k, const int j, const int il, const int iu,
             const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
             AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez);
  void FluxToGlobal2(const int k, const int j, const int il, const int iu,
             const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
             AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez);
  void FluxToGlobal3(const int k, const int j, const int il, const int iu,
             const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
             AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez);

  // ...to raise (lower) covariant (contravariant) components of a vector
  void RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
                       int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3);
  void LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
                       Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3);

 protected:
  bool coarse_flag;  // true if this coordinate object is parent (coarse) mesh in AMR
  Mesh *pm;
  int il, iu, jl, ju, kl, ku, ng;  // limits of indices of arrays (normal or coarse)
  int nc1, nc2, nc3;               // # cells in each dir of arrays (normal or coarse)
  // Scratch arrays for coordinate factors
  // Format: coord_<type>[<direction>]_<index>[<count>]_
  //   type: vol[ume], area, etc.
  //   direction: 1/2/3 depending on which face, edge, etc. is in play
  //   index: i/j/k indicating which coordinates index array
  //   count: 1/2/... in case multiple arrays are needed for different terms
  AthenaArray<Real> coord_vol_i_, coord_vol_i1_, coord_vol_i2_;
  AthenaArray<Real> coord_vol_j_, coord_vol_j1_, coord_vol_j2_;
  AthenaArray<Real> coord_vol_k1_;
  AthenaArray<Real> coord_vol_kji_;
  AthenaArray<Real> coord_area1_i_, coord_area1_i1_;
  AthenaArray<Real> coord_area1_j_, coord_area1_j1_, coord_area1_j2_;
  AthenaArray<Real> coord_area1_k1_;
  AthenaArray<Real> coord_area1_kji_;
  AthenaArray<Real> coord_area2_i_, coord_area2_i1_, coord_area2_i2_;
  AthenaArray<Real> coord_area2_j_, coord_area2_j1_, coord_area2_j2_;
  AthenaArray<Real> coord_area2_k1_;
  AthenaArray<Real> coord_area2_kji_;
  AthenaArray<Real> coord_area3_i_, coord_area3_i1_, coord_area3_i2_;
  AthenaArray<Real> coord_area3_j1_, coord_area3_j2_;
  AthenaArray<Real> coord_area3_kji_;
  AthenaArray<Real> coord_area1vc_i_,coord_area1vc_j_; //nonidealmhd additions
  AthenaArray<Real> coord_area2vc_i_,coord_area2vc_j_; //nonidealmhd additions
  AthenaArray<Real> coord_area3vc_i_; //nonidealmhd addition
  AthenaArray<Real> coord_len1_i1_, coord_len1_i2_;
  AthenaArray<Real> coord_len1_j1_, coord_len1_j2_;
  AthenaArray<Real> coord_len1_kji_;
  AthenaArray<Real> coord_len2_i1_;
  AthenaArray<Real> coord_len2_j1_, coord_len2_j2_;
  AthenaArray<Real> coord_len2_kji_;
  AthenaArray<Real> coord_len3_i1_;
  AthenaArray<Real> coord_len3_j1_, coord_len3_j2_;
  AthenaArray<Real> coord_len3_k1_;
  AthenaArray<Real> coord_len3_kji_;
  AthenaArray<Real> coord_width1_i1_;
  AthenaArray<Real> coord_width1_kji_;
  AthenaArray<Real> coord_width2_i1_;
  AthenaArray<Real> coord_width2_j1_;
  AthenaArray<Real> coord_width2_kji_;
  AthenaArray<Real> coord_width3_j1_, coord_width3_j2_, coord_width3_j3_;
  AthenaArray<Real> coord_width3_k1_;
  AthenaArray<Real> coord_width3_ji1_;
  AthenaArray<Real> coord_width3_kji_;
  AthenaArray<Real> coord_src_j1_, coord_src_j2_;
  AthenaArray<Real> coord_src_kji_;
  AthenaArray<Real> coord_src1_i_;
  AthenaArray<Real> coord_src1_j_;
  AthenaArray<Real> coord_src2_i_;
  AthenaArray<Real> coord_src2_j_;
  AthenaArray<Real> coord_src3_j_;

  // Scratch arrays for physical source terms
  AthenaArray<Real> phy_src1_i_, phy_src2_i_;

  // GR-specific scratch arrays
  AthenaArray<Real> metric_cell_i1_, metric_cell_i2_;
  AthenaArray<Real> metric_cell_j1_, metric_cell_j2_;
  AthenaArray<Real> metric_cell_kji_;
  AthenaArray<Real> metric_face1_i1_, metric_face1_i2_;
  AthenaArray<Real> metric_face1_j1_, metric_face1_j2_;
  AthenaArray<Real> metric_face1_kji_;
  AthenaArray<Real> metric_face2_i1_, metric_face2_i2_;
  AthenaArray<Real> metric_face2_j1_, metric_face2_j2_;
  AthenaArray<Real> metric_face2_kji_;
  AthenaArray<Real> metric_face3_i1_, metric_face3_i2_;
  AthenaArray<Real> metric_face3_j1_, metric_face3_j2_;
  AthenaArray<Real> metric_face3_kji_;
  AthenaArray<Real> trans_face1_i1_, trans_face1_i2_;
  AthenaArray<Real> trans_face1_j1_;
  AthenaArray<Real> trans_face1_ji1_, trans_face1_ji2_, trans_face1_ji3_,
    trans_face1_ji4_, trans_face1_ji5_, trans_face1_ji6_, trans_face1_ji7_;
  AthenaArray<Real> trans_face1_kji_;
  AthenaArray<Real> trans_face2_i1_, trans_face2_i2_;
  AthenaArray<Real> trans_face2_j1_;
  AthenaArray<Real> trans_face2_ji1_, trans_face2_ji2_, trans_face2_ji3_,
    trans_face2_ji4_, trans_face2_ji5_, trans_face2_ji6_;
  AthenaArray<Real> trans_face2_kji_;
  AthenaArray<Real> trans_face3_i1_, trans_face3_i2_;
  AthenaArray<Real> trans_face3_j1_;
  AthenaArray<Real> trans_face3_ji1_, trans_face3_ji2_, trans_face3_ji3_,
    trans_face3_ji4_, trans_face3_ji5_, trans_face3_ji6_;
  AthenaArray<Real> trans_face3_kji_;
  AthenaArray<Real> g_, gi_;

  // GR-specific variables
  Real bh_mass_;
  Real bh_spin_;
};

#endif // COORDINATES_COORDINATES_HPP_
