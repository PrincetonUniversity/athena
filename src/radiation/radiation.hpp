#ifndef RADIATION_RADIATION_HPP_
#define RADIATION_RADIATION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

// Forward declarations
class MeshBlock;
class ParameterInput;

//----------------------------------------------------------------------------------------
// Radiation class
// Notes:
//   designed for general relativity

class Radiation {

public:

  // Constructor and destructor
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();

  // Object pointers
  MeshBlock* pmy_block;  // pointer to containing MeshBlock

  // Flags
  bool source_terms_defined;

  // Parameters
  int nzeta;     // number of polar radiation angles in active zone
  int npsi;      // number of azimuthal radiation angles in active zone
  int nang;      // total number of radiation angles, including ghost zones
  int nang_zf;   // total number of radiation angles when zeta is on faces
  int nang_pf;   // total number of radiation angles when psi is on faces
  int nang_zpf;  // total number of radiation angles when zeta and psi are on faces
  int zs, ze;    // start and end zeta-indices
  int ps, pe;    // start and end psi-indices
  int is, ie;    // start and end x1-indices
  int js, je;    // start and end x2-indices
  int ks, ke;    // start and end x3-indices

  // Data arrays
  AthenaArray<Real> zetaf;        // face-centered polar radiation angles
  AthenaArray<Real> zetav;        // volume-centered polar radiation angles
  AthenaArray<Real> dzetaf;       // face-to-face polar radiation angle differences
  AthenaArray<Real> psif;         // face-centered azimuthal radiation angles
  AthenaArray<Real> psiv;         // volume-centered azimuthal radiation angles
  AthenaArray<Real> dpsif;        // face-to-face azimuthal radiation angle differences
  AthenaArray<Real> zeta_length;  // angular length at constant psi
  AthenaArray<Real> psi_length;   // angular length at constant zeta
  AthenaArray<Real> solid_angle;  // angular area of cell
  AthenaArray<Real> prim;         // primitive intensity I
  AthenaArray<Real> prim1;        // primitive intensity I, for substeps
  AthenaArray<Real> cons;         // conserved intensity n^0 I
  AthenaArray<Real> cons1;        // conserved intensity n^0 I, for substeps
  AthenaArray<Real> cons2;        // conserved intensity n^0 I, for substeps
  AthenaArray<Real> flux_x[3];    // spatial fluxes of intensity n^i I
  AthenaArray<Real> flux_a[2];    // angular fluxes of intensity n^a I

  // Function pointers
  RadSrcTermFunc_t UserSourceTerm;

  // Task list functions
  void WeightedAveCons(AthenaArray<Real> &cons_out, AthenaArray<Real> &cons_in_1,
      AthenaArray<Real> &cons_in_2, const Real weights[3]);
  void CalculateFluxes(AthenaArray<Real> &prim_in, int order);
  void AddFluxDivergenceToAverage(AthenaArray<Real> &prim_in, const Real weight,
      AthenaArray<Real> &cons_out);
  void PrimitiveToConserved(const AthenaArray<Real> &prim_in, AthenaArray<Real> &cons_out,
      Coordinates *pcoord, int il, int iu, int jl, int ju, int kl, int ku);
  void ConservedToPrimitive(AthenaArray<Real> &cons_in, AthenaArray<Real> &prim_out,
      Coordinates *pcoord, int il, int iu, int jl, int ju, int kl, int ku);
  void AddSourceTerms(const Real time, const Real dt, const AthenaArray<Real> &prim_in,
      AthenaArray<Real> &cons_out);

  // Other functions
  int AngleInd(int l, int m, bool zeta_face = false, bool psi_face = false);
  void CalculateBeamSource(Real pos_1, Real pos_2, Real pos_3, Real width, Real dir_1,
      Real dir_2, Real dir_3, Real spread, AthenaArray<Real> &prim_vals,
      AthenaArray<Real> &cons_vals, bool spherical = false);
  void SetMoments(AthenaArray<Real> &moments);

private:

  // Data arrays
  AthenaArray<Real> n0_;        // n^0 at cell and angle centers
  AthenaArray<Real> n1_;        // n^1 at x^1-faces and angle centers
  AthenaArray<Real> n2_;        // n^2 at x^2-faces and angle centers
  AthenaArray<Real> n3_;        // n^3 at x^3-faces and angle centers
  AthenaArray<Real> na0_;       // -n^ah n^bh omega^0h_{ah,bh} at cell and angle centers
  AthenaArray<Real> na1_;       // n^zeta at cell centers and zeta-faces
  AthenaArray<Real> na2_;       // n^psi at cell centers and psi-faces
  AthenaArray<Real> prim_l_;    // left reconstructed state
  AthenaArray<Real> prim_r_;    // right reconstructed state
  AthenaArray<Real> area_l_;    // left face areas
  AthenaArray<Real> area_r_;    // right face areas
  AthenaArray<Real> vol_;       // cell volumes
  AthenaArray<Real> flux_div_;  // flux divergences in spatial coordinates
};

#endif // RADIATION_RADIATION_HPP_
