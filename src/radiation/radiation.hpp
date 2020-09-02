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
#include "../athena.hpp"                  // Real, indices, function prototypes
#include "../athena_arrays.hpp"           // AthenaArray
#include "../bvals/cc/rad/bvals_rad.hpp"  // RadBoundaryVariable
#include "../mesh/mesh.hpp"               // MeshBlock
#include "../parameter_input.hpp"         // ParameterInput

//----------------------------------------------------------------------------------------
// Radiation class
// Notes:
//   Currently designed for general relativity.

class Radiation {

public:

  // Constructor and destructor
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();

  // Object and function pointers
  MeshBlock* pmy_block;
  RadSrcTermFunc UserSourceTerm;
  OpacityFunc UpdateOpacity;

  // Flags
  bool coupled_to_matter;
  bool source_terms_defined;
  bool moment_fix;

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

  // Physical constants in CGS
  Real c_cgs = 2.99792458e10;                  // speed of light in cm/s
  Real m_p_cgs = 1.67262192369e-24;            // proton mass in g
  Real k_b_cgs = 1.380649e-16;                 // Boltzmann constant in erg/K
  Real sigma_sb_cgs = 5.670374419e-5;          // Stefan-Boltzmann constant in
                                               // erg/(cm^2*s*K^4)
  Real arad_cgs = 4.0 * sigma_sb_cgs / c_cgs;  // radiation constant in erg/(cm^3*K^4)

  // Physical constants in code units
  Real arad;  // radiation constant in code units

  // User-specified units
  Real density_cgs;  // code unit of density in g/cm^3
  Real mol_weight;   // molecular weight of gas in proton masses

  // Data arrays
  AthenaArray<Real> zetaf;           // face-centered polar radiation angles
  AthenaArray<Real> zetav;           // volume-centered polar radiation angles
  AthenaArray<Real> dzetaf;          // face-to-face polar radiation angle differences
  AthenaArray<Real> psif;            // face-centered azimuthal radiation angles
  AthenaArray<Real> psiv;            // volume-centered azimuthal radiation angles
  AthenaArray<Real> dpsif;           // face-to-face azimuthal radiation angle differences
  AthenaArray<Real> zeta_length;     // angular length at constant psi
  AthenaArray<Real> psi_length;      // angular length at constant zeta
  AthenaArray<Real> solid_angle;     // angular area of cell
  AthenaArray<Real> prim;            // primitive intensity I
  AthenaArray<Real> prim1;           // primitive intensity I, for substeps
  AthenaArray<Real> cons;            // conserved intensity n^0 n_0 I
  AthenaArray<Real> cons1;           // conserved intensity n^0 n_0 I, for substeps
  AthenaArray<Real> cons2;           // conserved intensity n^0 n_0 I, for substeps
  AthenaArray<Real> flux_x[3];       // spatial fluxes of intensity n^i n_0 I
  AthenaArray<Real> flux_a[2];       // angular fluxes of intensity n^a n_0 I
  AthenaArray<Real> coarse_prim;     // prolongation/restriction buffer
  AthenaArray<Real> coarse_cons;     // prolongation/restriction buffer
  AthenaArray<Real> opacity;         // opacity array
  AthenaArray<Real> moments_coord;   // contravariant stress tensor in coordinate frame
  AthenaArray<Real> moments_tetrad;  // contravariant stress tensor in tetrad frame
  AthenaArray<Real> moments_fluid;   // contravariant stress tensor in fluid frame

  // Boundary communication
  RadBoundaryVariable rbvar;
  int refinement_idx{-1};

  // Task list functions - core
  void WeightedAve(AthenaArray<Real> &cons_out, AthenaArray<Real> &cons_in_1,
      AthenaArray<Real> &cons_in_2, const Real weights[3]);
  void PrimitiveToConserved(const AthenaArray<Real> &prim_in, AthenaArray<Real> &cons_out,
      Coordinates *pcoord, int il, int iu, int jl, int ju, int kl, int ku);
  void ConservedToPrimitive(AthenaArray<Real> &cons_in, AthenaArray<Real> &prim_out,
      int il, int iu, int jl, int ju, int kl, int ku);
  void ConservedToPrimitiveWithMoments(AthenaArray<Real> &cons_in,
      AthenaArray<Real> &prim_out, const AthenaArray<Real> &prim_hydro,
      Coordinates *pcoord, int il, int iu, int jl, int ju, int kl, int ku);
  void AddSourceTerms(const Real time, const Real dt, const AthenaArray<Real> &prim_rad,
      const AthenaArray<Real> &prim_hydro, AthenaArray<Real> &cons_rad,
      AthenaArray<Real> &cons_hydro);

  // Task list functions - fluxes (defined in rad_fluxes.cpp)
  void CalculateFluxes(AthenaArray<Real> &prim_rad, const AthenaArray<Real> &prim_hydro,
      int order);
  void AddFluxDivergenceToAverage(AthenaArray<Real> &prim_in, const Real weight,
      AthenaArray<Real> &cons_out);

  // Reconstruction functions (defined in rad_reconstruction.cpp)
  void RadiationDonorCellX1(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationDonorCellX2(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationDonorCellX3(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationDonorCellA1(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationDonorCellA2(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationPiecewiseLinearX1(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationPiecewiseLinearX2(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationPiecewiseLinearX3(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationPiecewiseLinearA1(const AthenaArray<Real> &intensity, int k, int j);
  void RadiationPiecewiseLinearA2(const AthenaArray<Real> &intensity, int k, int j);

  // Setup functions for problem generators (defined in rad_setup.cpp)
  void CalculateBeamSource(Real pos_1, Real pos_2, Real pos_3, Real width, Real dir_1,
      Real dir_2, Real dir_3, Real spread, Real dii_dt, AthenaArray<Real> &dcons_dt,
      bool cylindrical = false, bool spherical = false);
  void CalculateConstantRadiation(Real energy, Real u1, Real u2, Real u3,
      AthenaArray<Real> &cons_out);
  void CalculateRadiationInCellM1(Real energy, Real u1, Real u2, Real u3, int k, int j,
      int i, const AthenaArray<Real> &g, AthenaArray<Real> &cons_out);
  void CalculateRadiationInCellLinear(Real ee_f, Real ff1_f, Real ff2_f, Real ff3_f,
      Real uu1, Real uu2, Real uu3, int k, int j, int i, const AthenaArray<Real> &g,
      AthenaArray<Real> &cons_out);

  // Other functions
  int AngleInd(int l, int m, bool zeta_face = false, bool psi_face = false);
  void SetMoments(const AthenaArray<Real> &prim_hydro, Coordinates *pcoord, int il,
      int iu, int jl, int ju, int kl, int ku);
  void EnrollOpacityFunction(OpacityFunc MyOpacityFunction);

private:

  // Data arrays
  AthenaArray<Real> nh_cc_;          // n^\hat{mu} at angle centers
  AthenaArray<Real> nmu_;            // n^mu at cell and angle centers
  AthenaArray<Real> n0_n_mu_;        // n^0 n_mu at cell and angle centers
  AthenaArray<Real> n1_n_0_;         // n^1 n_0 at x^1-faces and angle centers
  AthenaArray<Real> n2_n_0_;         // n^2 n_0 at x^2-faces and angle centers
  AthenaArray<Real> n3_n_0_;         // n^3 n_0 at x^3-faces and angle centers
  AthenaArray<Real> na1_n_0_;        // n^zeta n_0 at cell centers and zeta-faces
  AthenaArray<Real> na2_n_0_;        // n^psi n_0 at cell centers and psi-faces
  AthenaArray<Real> rad_l_;          // left reconstructed radiation state
  AthenaArray<Real> rad_r_;          // right reconstructed radiation state
  AthenaArray<Real> area_l_;         // left face areas
  AthenaArray<Real> area_r_;         // right face areas
  AthenaArray<Real> vol_;            // cell volumes
  AthenaArray<Real> flux_div_;       // flux divergences in spatial coordinates
  AthenaArray<Real> g_, gi_;         // metric and inverse
  AthenaArray<Real> norm_to_tet_;    // transformation from normal to tetrad frame
  AthenaArray<Real> moments_old_;    // moments of radiation field before fluid coupling
  AthenaArray<Real> moments_new_;    // moments of radiation field after fluid coupling
  AthenaArray<Real> u_tet_;          // fluid 4-velocity in tetrad frame
  AthenaArray<Real> coefficients_;   // quartic coefficients for implicit update
  AthenaArray<bool> bad_cell_;       // flag indicating problem with coupling
  AthenaArray<Real> tt_plus_;        // gas temperature after coupling
  AthenaArray<Real> ee_f_minus_;     // fluid-frame radiation energy before coupling
  AthenaArray<Real> ee_f_plus_;      // fluid-frame radiation energy after coupling
};

#endif // RADIATION_RADIATION_HPP_
