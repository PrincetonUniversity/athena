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
#include "../mesh.hpp"                    // MeshBlock
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
  bool using_planck_mean;

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
  Real m_p;       // proton mass in code units
  Real k_b;       // Boltzmann constant in code units
  Real sigma_sb;  // Stefan-Boltzmann constant in code units
  Real arad;      // radiation constant in code units

  // User-specified units
  Real length_cgs;   // code unit of length in cm
  Real density_cgs;  // code unit of density in g/cm^3
  Real mol_weight;   // molecular weight of gas in proton masses

  // Inferred units
  Real time_cgs;         // code unit of time in s
  Real pressure_cgs;     // code unit of pressure (and energy density) in dyne/cm^2
  Real temperature_cgs;  // code unit of temperature in K
  Real intensity_cgs;    // code unit of intensity I in erg/(cm^2*s*sr)
  Real opacity_cgs;      // code unit of opacity in cm^2/g

  // Emission, absorption, and scattering coefficients in CGS
  Real kappa_cgs;  // opacity in cm^2/g

  // Emission, absorption, and scattering coefficients in code units
  Real kappa;  // opacity in code units

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
  AthenaArray<Real> cons;         // conserved intensity n^0 n_0 I
  AthenaArray<Real> cons1;        // conserved intensity n^0 n_0 I, for substeps
  AthenaArray<Real> cons2;        // conserved intensity n^0 n_0 I, for substeps
  AthenaArray<Real> flux_x[3];    // spatial fluxes of intensity n^i n_0 I
  AthenaArray<Real> flux_a[2];    // angular fluxes of intensity n^a n_0 I
  AthenaArray<Real> coarse_prim;  // prolongation/restriction buffer
  AthenaArray<Real> coarse_cons;  // prolongation/restriction buffer
  AthenaArray<Real> opacity;      // opacity array

  // Boundary communication
  RadBoundaryVariable rbvar;
  int refinement_idx{-1};

  // Task list functions
  void WeightedAve(AthenaArray<Real> &cons_out, AthenaArray<Real> &cons_in_1,
      AthenaArray<Real> &cons_in_2, const Real weights[3]);
  void CalculateFluxes(AthenaArray<Real> &prim_in, int order);
  void AddFluxDivergenceToAverage(AthenaArray<Real> &prim_in, const Real weight,
      AthenaArray<Real> &cons_out);
  void PrimitiveToConserved(const AthenaArray<Real> &prim_in, AthenaArray<Real> &cons_out,
      Coordinates *pcoord, int il, int iu, int jl, int ju, int kl, int ku);
  void ConservedToPrimitive(AthenaArray<Real> &cons_in, AthenaArray<Real> &prim_out,
      Coordinates *pcoord, int il, int iu, int jl, int ju, int kl, int ku);
  void AddSourceTerms(const Real time, const Real dt, const AthenaArray<Real> &prim_rad,
      const AthenaArray<Real> &prim_hydro, AthenaArray<Real> &cons_rad,
      AthenaArray<Real> &cons_hydro);

  // Fluid coupling functions
  void EnrollOpacityFunction(OpacityFunc MyOpacityFunction);
  void Coupling(const AthenaArray<Real> &prim_hydro, const AthenaArray<Real> &normal,
      const AthenaArray<Real> &omega, Real dt, int k, int j,
      AthenaArray<Real> &intensity);

  // Other functions
  int AngleInd(int l, int m, bool zeta_face = false, bool psi_face = false);
  void CalculateBeamSource(Real pos_1, Real pos_2, Real pos_3, Real width, Real dir_1,
      Real dir_2, Real dir_3, Real spread, Real dii_dt, AthenaArray<Real> &dcons_dt,
      bool cylindrical = false, bool spherical = false);
  void SetMoments(AthenaArray<Real> &moments);

private:

  // Data arrays
  AthenaArray<Real> n0_n_mu_;        // n^0 n_mu at cell and angle centers
  AthenaArray<Real> n1_n_0_;         // n^1 n_0 at x^1-faces and angle centers
  AthenaArray<Real> n2_n_0_;         // n^2 n_0 at x^2-faces and angle centers
  AthenaArray<Real> n3_n_0_;         // n^3 n_0 at x^3-faces and angle centers
  AthenaArray<Real> na1_n_0_;        // n^zeta n_0 at cell centers and zeta-faces
  AthenaArray<Real> na2_n_0_;        // n^psi n_0 at cell centers and psi-faces
  AthenaArray<Real> prim_l_;         // left reconstructed state
  AthenaArray<Real> prim_r_;         // right reconstructed state
  AthenaArray<Real> area_l_;         // left face areas
  AthenaArray<Real> area_r_;         // right face areas
  AthenaArray<Real> vol_;            // cell volumes
  AthenaArray<Real> flux_div_;       // flux divergences in spatial coordinates
  AthenaArray<Real> intensity_scr_;  // intensity in a single cell
  AthenaArray<Real> tran_coef_;      // transformation coefficient in a single cell
  AthenaArray<Real> weight_;         // weight (fractional solid angle) in a single cell
  AthenaArray<Real> vncsigma2_;      // source terms in a single cell
  AthenaArray<Real> g_, gi_;         // metric and inverse
  AthenaArray<Real> norm_to_tet_;    // transformation from normal to tetrad frame
  AthenaArray<Real> u_tet_;          // fluid 4-velocity in tetrad frame
  AthenaArray<Real> dtau_;           // timestep in fluid frame
  AthenaArray<Real> weight_sum_;     // sum of solid angles, used for normalizing
  AthenaArray<Real> n_cm_;           // unit null direction in comoving fluid frame
  AthenaArray<Real> omega_cm_;       // solid angle in comoving fluid frame
  AthenaArray<Real> intensity_cm_;   // intensity I in comoving fluid frame
  AthenaArray<Real> moments_old_;    // moments of radiation field before fluid coupling
  AthenaArray<Real> moments_new_;    // moments of radiation field after fluid coupling
};

#endif // RADIATION_RADIATION_HPP_
