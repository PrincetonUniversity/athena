// General relativistic Fishbone-Moncrief torus generator

// Primary header
#include "../mesh/mesh.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // abs(), cos(), exp(), log(), NAN, pow(), sin(), sqrt()
#include <iostream>   // endl
#include <limits>     // numeric_limits::max()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str(), string

// Athena headers
#include "../athena.hpp"                   // macros, enums, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro

// Declarations
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void InflowBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
static Real CalculateLFromRPeak(Real r);
static Real CalculateRPeakFromL(Real l_target);
static Real LogHAux(Real r, Real sin_theta);
static void CalculateVelocityInTorus(Real r, Real sin_theta, Real *pu0, Real *pu3);
static Real CalculateBetaMin();
static bool CalculateBeta(Real r_m, Real r_c, Real r_p, Real theta_m, Real theta_c,
    Real theta_p, Real *pbeta);
static bool CalculateBetaFromA(Real r_m, Real r_c, Real r_p, Real theta_m, Real theta_c,
    Real theta_p, Real a_cm, Real a_cp, Real a_mc, Real a_pc, Real *pbeta);
static Real CalculateMagneticPressure(Real bb1, Real bb2, Real bb3, Real r,
    Real sin_theta, Real cos_theta);

// Global variables
static Real m, a;                                // black hole parameters
static Real gamma_adi, k_adi;                    // hydro parameters
static Real r_edge, r_peak, l, rho_max;          // fixed torus parameters
static Real log_h_edge, log_h_peak;              // calculated torus parameters
static Real pgas_over_rho_peak, rho_peak;        // more calculated torus parameters
static Real rho_min, rho_pow, u_min, u_pow;      // background parameters
static Real potential_cutoff;                    // sets region of torus to magnetize
static Real potential_r_pow, potential_rho_pow;  // set how vector potential scales
static Real beta_min;                            // min ratio of gas to mag pressure
static std::string field_config;                 // type of magnetic field
static int sample_n_r, sample_n_theta;           // number of cells in sample grid
static Real sample_r_rat;                        // sample grid geometric spacing ratio
static Real sample_cutoff;                       // density cutoff for sample grid
static Real x1_min, x1_max, x2_min, x2_max;      // limits in chosen coordinate system
static Real r_min, r_max, theta_min, theta_max;  // limits in r,theta
static Real pert_amp, pert_kr, pert_kz;          // parameters for initial perturbations
static AthenaArray<Real> g, gi;                  // metric and its inverse

//--------------------------------------------------------------------------------------

// Function for setting up arrays to handle user work
// Inputs:
//   pin: input parameters (unused)
// Outputs: (none)
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Read problem-specific properties from input file
  rho_min = pin->GetReal("hydro", "rho_min");
  rho_pow = pin->GetReal("hydro", "rho_pow");
  u_min = pin->GetReal("hydro", "u_min");
  u_pow = pin->GetReal("hydro", "u_pow");
  k_adi = pin->GetReal("problem", "k_adi");
  r_edge = pin->GetReal("problem", "r_edge");
  r_peak = pin->GetReal("problem", "r_peak");
  l = pin->GetReal("problem", "l");
  rho_max = pin->GetReal("problem", "rho_max");
  if (MAGNETIC_FIELDS_ENABLED)
  {
    potential_cutoff = pin->GetReal("problem", "potential_cutoff");
    potential_r_pow = pin->GetReal("problem", "potential_r_pow");
    potential_rho_pow = pin->GetReal("problem", "potential_rho_pow");
    beta_min = pin->GetReal("problem", "beta_min");
    field_config = pin->GetString("problem", "field_config");
    sample_n_r = pin->GetInteger("problem", "sample_n_r");
    sample_n_theta = pin->GetInteger("problem", "sample_n_theta");
    sample_r_rat = pin->GetReal("problem", "sample_r_rat");
    sample_cutoff = pin->GetReal("problem", "sample_cutoff");
    x1_min = pin->GetReal("mesh", "x1min");
    x1_max = pin->GetReal("mesh", "x1max");
    x2_min = pin->GetReal("mesh", "x2min");
    x2_max = pin->GetReal("mesh", "x2max");
  }
  pert_amp = pin->GetOrAddReal("problem", "pert_amp", 0.0);
  pert_kr = pin->GetOrAddReal("problem", "pert_kr", 0.0);
  pert_kz = pin->GetOrAddReal("problem", "pert_kz", 0.0);

  // Prepare arrays if needed for extra outputs
  if (NIFOV == 1 or NIFOV == 5 or
      (MAGNETIC_FIELDS_ENABLED and (NIFOV == 2 or NIFOV == 10)))
  {
    g.NewAthenaArray(NMETRIC, mesh_size.nx1/nrbx1+NGHOST);
    gi.NewAthenaArray(NMETRIC, mesh_size.nx1/nrbx1+NGHOST);
  }

  // Enroll boundary functions
  EnrollUserBoundaryFunction(INNER_X1, InflowBoundary);
  EnrollUserBoundaryFunction(OUTER_X1, FixedBoundary);
  return;
}

//--------------------------------------------------------------------------------------

// Function for freeing arrays needed for user work
// Inputs:
//   pin: parameters
// Outputs: (none)
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  if (NIFOV == 1 or NIFOV == 5 or
      (MAGNETIC_FIELDS_ENABLED and (NIFOV == 2 or NIFOV == 10)))
  {
    g.DeleteAthenaArray();
    gi.DeleteAthenaArray();
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   initializes Fishbone-Moncrief torus
//     sets both primitive and conserved variables
//   defines and enrolls fixed r- and theta-direction boundary conditions
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//              Fishbone 1977, ApJ 215 323 (F)
//   assumes x3 is axisymmetric direction
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Get mass and spin of black hole
  m = pcoord->GetMass();
  a = pcoord->GetSpin();

  // Get ratio of specific heats
  gamma_adi = peos->GetGamma();

  // Reset whichever of l,r_peak is not specified
  if (r_peak >= 0.0)
    l = CalculateLFromRPeak(r_peak);
  else
    r_peak = CalculateRPeakFromL(l);

  // Prepare scratch arrays
  AthenaArray<bool> in_torus;
  in_torus.NewAthenaArray(ju+1, iu+1);
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);

  // Initialize primitive values
  log_h_edge = LogHAux(r_edge, 1.0);
  log_h_peak = LogHAux(r_peak, 1.0) - log_h_edge;
  pgas_over_rho_peak = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_peak)-1.0);
  rho_peak = std::pow(pgas_over_rho_peak/k_adi, 1.0/(gamma_adi-1.0)) / rho_max;
  for (int j = jl; j <= ju; ++j)
  {
    pcoord->CellMetric(kl, j, il, iu, g, gi);
    for (int i = il; i <= iu; ++i)
    {
      // Get Boyer-Lindquist coordinates of cell
      Real r, theta, phi;
      pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
          pcoord->x3v(kl), &r, &theta, &phi);
      Real sin_theta = std::sin(theta);

      // Determine if we are in the torus
      Real log_h;
      in_torus(j,i) = false;
      if (r >= r_edge)
      {
        log_h = LogHAux(r, sin_theta) - log_h_edge;  // (FM 3.6)
        if (log_h >= 0.0)
          in_torus(j,i) = true;
      }

      // Calculate primitives depending on location
      Real rho, pgas, uu1, uu2, uu3;
      if (in_torus(j,i))
      {
        Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
        rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
        pgas = pgas_over_rho * rho;
        Real u0, u3;
        CalculateVelocityInTorus(r, sin_theta, &u0, &u3);
        Real u0_pref, u1_pref, u2_pref, u3_pref;
        pcoord->TransformVectorCell(u0, 0.0, 0.0, u3, kl, j, i, &u0_pref, &u1_pref,
            &u2_pref, &u3_pref);
        uu1 = u1_pref - gi(I01,i)/gi(I00,i) * u0_pref;
        uu2 = u2_pref - gi(I02,i)/gi(I00,i) * u0_pref;
        uu3 = u3_pref - gi(I03,i)/gi(I00,i) * u0_pref;
      }
      else
      {
        rho = rho_min * std::pow(r, rho_pow);
        Real u = u_min * std::pow(r, u_pow);
        pgas = (gamma_adi-1.0) * u;
        uu1 = 0.0;
        uu2 = 0.0;
        uu3 = 0.0;
      }

      // Set primitive values, including cylindrically symmetric radial velocity
      // perturbations
      Real rr = r * std::abs(sin_theta);
      Real z = r * std::cos(theta);
      Real amp_rel = pert_amp * std::sin(pert_kr*rr) * std::cos(pert_kz*z);
      Real amp_abs = amp_rel * uu3;
      Real pert_uur = rr/r * amp_abs;
      Real pert_uutheta = std::cos(theta)/r * amp_abs;
      for (int k = kl; k <= ku; ++k)
      {
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = uu1 + pert_uur;
        phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = uu2 + pert_uutheta;
        phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = uu3;
      }
    }
  }

  // Free scratch arrays
  in_torus.DeleteAthenaArray();
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Initialize magnetic fields
  if (MAGNETIC_FIELDS_ENABLED)
  {
    // Determine limits of sample grid
    Real r1, r2, r3, r4, theta1, theta2, theta3, theta4, temp;
    pcoord->GetBoyerLindquistCoordinates(x1_min, x2_min, pcoord->x3v(kl), &r1, &theta1,
        &temp);
    pcoord->GetBoyerLindquistCoordinates(x1_max, x2_min, pcoord->x3v(kl), &r2, &theta2,
        &temp);
    pcoord->GetBoyerLindquistCoordinates(x1_min, x2_max, pcoord->x3v(kl), &r3, &theta3,
        &temp);
    pcoord->GetBoyerLindquistCoordinates(x1_max, x2_max, pcoord->x3v(kl), &r4, &theta4,
        &temp);
    r_min = std::min(std::min(r1, r2), std::min(r3, r4));
    r_max = std::max(std::max(r1, r2), std::max(r3, r4));
    theta_min = std::min(std::min(theta1, theta2), std::min(theta3, theta4));
    theta_max = std::max(std::max(theta1, theta2), std::max(theta3, theta4));

    // Prepare 2D arrays of vector potential values
    AthenaArray<Real> a_phi_edges, a_phi_cells;
    a_phi_edges.NewAthenaArray(ju+2, iu+2);
    a_phi_cells.NewAthenaArray(ju+1, iu+1);
    Real normalization;

    // Set vector potential in normal case
    if (field_config.compare("normal") == 0)
    {
      // Set edge-centered vector potential values
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2f(j),
              pcoord->x3v(kl), &r, &theta, &phi);
          if (r >= r_edge)
          {
            Real log_h = LogHAux(r, std::sin(theta)) - log_h_edge;  // (FM 3.6)
            if (log_h >= 0.0)
            {
              Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
              Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
              Real rho_cutoff = std::max(rho-potential_cutoff, 0.0);
              a_phi_edges(j,i) = std::pow(r,potential_r_pow)
                  * std::pow(rho_cutoff,potential_rho_pow);
            }
          }
        }

      // Set cell-centered vector potential values
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
              pcoord->x3v(kl), &r, &theta, &phi);
          if (r >= r_edge)
          {
            Real log_h = LogHAux(r, std::sin(theta)) - log_h_edge;  // (FM 3.6)
            if (log_h >= 0.0)
            {
              Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
              Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
              Real rho_cutoff = std::max(rho-potential_cutoff, 0.0);
              a_phi_cells(j,i) = std::pow(r,potential_r_pow)
                  * std::pow(rho_cutoff,potential_rho_pow);
            }
          }
        }

      // Calculate magnetic field normalization
      if (beta_min < 0.0)
        normalization = 0.0;
      else
      {
        Real beta_min_actual = CalculateBetaMin();
        normalization = std::sqrt(beta_min_actual/beta_min);
      }
    }

    // Set vector potential in renormalized case
    else if (field_config.compare("renorm") == 0)
    {
      // Prepare global 2D sample arrays for integrating
      AthenaArray<Real> r_face, r_cell, theta_face, theta_cell;
      r_face.NewAthenaArray(sample_n_r+1);
      r_cell.NewAthenaArray(sample_n_r);
      theta_face.NewAthenaArray(sample_n_theta/2+1);
      theta_cell.NewAthenaArray(sample_n_theta/2);
      AthenaArray<Real> a_phi_global_edges, a_phi_global_cells;
      a_phi_global_edges.NewAthenaArray(sample_n_theta/2+1, sample_n_r+1);
      a_phi_global_cells.NewAthenaArray(sample_n_theta/2, sample_n_r);
      AthenaArray<Real> bbr_r_faces, bbr_theta_faces;
      bbr_r_faces.NewAthenaArray(sample_n_theta/2, sample_n_r+1);
      bbr_theta_faces.NewAthenaArray(sample_n_theta/2+1, sample_n_r);

      // Calculate r values
      Real delta_r;
      for (int i = 0; i < sample_n_r; ++i)
      {
        if (i == 0)
        {
          r_face(i) = r_min;
          Real ratio_power = 1.0;
          Real ratio_sum = 1.0;
          for (int ii = 1; ii < sample_n_r; ++ii)
          {
            ratio_power *= sample_r_rat;
            ratio_sum += ratio_power;
          }
          delta_r = (r_max-r_min) / ratio_sum;
        }
        else
          delta_r *= sample_r_rat;
        r_face(i+1) = r_face(i) + delta_r;
        r_cell(i) = 0.5 * (r_face(i) + r_face(i+1));
      }

      // Calculate theta values
      for (int j = 0; j < sample_n_theta/2; ++j)
      {
        if (j == 0)
          theta_face(j) = theta_min;
        theta_face(j+1) = theta_min
            + static_cast<Real>(j+1)/static_cast<Real>(sample_n_theta)
            * (theta_max-theta_min);
        theta_cell(j) = 0.5 * (theta_face(j) + theta_face(j+1));
      }

      // Calculate edge-centered A_phi based on radius and density
      for (int j = 0; j < sample_n_theta/2+1; ++j)
        for (int i = 0; i < sample_n_r+1; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(r_face(i), theta_face(j),
              pcoord->x3v(kl), &r, &theta, &phi);
          Real rho = 0.0;
          if (r >= r_edge)
          {
            Real log_h = LogHAux(r, std::sin(theta)) - log_h_edge;  // (FM 3.6)
            if (log_h >= 0.0)
            {
              Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
              rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
            }
          }
          a_phi_global_edges(j,i) =
              std::pow(r,potential_r_pow) * std::pow(rho,potential_rho_pow);
        }

      // Calculate cell-centered A_phi based on radius and density
      for (int j = 0; j < sample_n_theta/2; ++j)
        for (int i = 0; i < sample_n_r; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(r_cell(i), theta_cell(j),
              pcoord->x3v(kl), &r, &theta, &phi);
          Real rho = 0.0;
          if (r >= r_edge)
          {
            Real log_h = LogHAux(r, std::sin(theta)) - log_h_edge;  // (FM 3.6)
            if (log_h >= 0.0)
            {
              Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
              rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
            }
          }
          a_phi_global_cells(j,i) =
              std::pow(r,potential_r_pow) * std::pow(rho,potential_rho_pow);
        }

      // Calculate r-face-centered B^r based on edge-centered A_phi and pgas
      for (int j = 0; j < sample_n_theta/2; ++j)
        for (int i = 1; i < sample_n_r; ++i)
        {
          Real r_m = r_cell(i-1);
          Real r_c = r_face(i);
          Real r_p = r_cell(i);
          Real theta_m = theta_face(j);
          Real theta_c = theta_cell(j);
          Real theta_p = theta_face(j+1);
          Real cos_theta = std::cos(theta_c);
          Real det = (SQR(r_c) + SQR(a) * SQR(cos_theta)) * std::abs(std::sin(theta_c));
          Real a_phi_cm = a_phi_global_edges(j,i);
          Real a_phi_cp = a_phi_global_edges(j+1,i);
          Real a_phi_mc = a_phi_global_cells(j,i-1);
          Real a_phi_pc = a_phi_global_cells(j,i);
          Real bbr = 1.0/det * (a_phi_cp-a_phi_cm) / (theta_p-theta_m);
          Real beta;
          bool value_set = CalculateBetaFromA(r_m, r_c, r_p, theta_m, theta_p, theta_c,
              a_phi_cm, a_phi_cp, a_phi_mc, a_phi_pc, &beta);
          if (value_set)
            bbr_r_faces(j,i) = bbr * std::sqrt(beta);
        }

      // Calculate theta-face-centered B^r based on cell-centered A_phi and pgas
      for (int j = 1; j < sample_n_theta/2; ++j)
        for (int i = 0; i < sample_n_r; ++i)
        {
          Real r_m = r_face(i);
          Real r_c = r_cell(i);
          Real r_p = r_face(i+1);
          Real theta_m = theta_cell(j-1);
          Real theta_c = theta_face(j);
          Real theta_p = theta_cell(j);
          Real cos_theta = std::cos(theta_c);
          Real det = (SQR(r_c) + SQR(a) * SQR(cos_theta)) * std::abs(std::sin(theta_c));
          Real a_phi_cm = a_phi_global_cells(j-1,i);
          Real a_phi_cp = a_phi_global_cells(j,i);
          Real a_phi_mc = a_phi_global_edges(j,i);
          Real a_phi_pc = a_phi_global_edges(j,i+1);
          Real bbr = 1.0/det * (a_phi_cp-a_phi_cm) / (theta_p-theta_m);
          Real beta;
          bool value_set = CalculateBetaFromA(r_m, r_c, r_p, theta_m, theta_p, theta_c,
              a_phi_cm, a_phi_cp, a_phi_mc, a_phi_pc, &beta);
          if (value_set)
            bbr_theta_faces(j,i) = bbr * std::sqrt(beta);
        }

      // Calculate edge-centered A_phi based on r-face-centered B^r
      for (int i = 0; i < sample_n_r+1; ++i)
      {
        Real r = r_face(i);
        a_phi_global_edges(0,i) = 0.0;
        for (int j = 1; j < sample_n_theta/2+1; ++j)
        {
          Real theta_m = theta_face(j-1);
          Real theta_c = theta_cell(j-1);
          Real theta_p = theta_face(j);
          Real cos_theta = std::cos(theta_c);
          Real det = (SQR(r) + SQR(a) * SQR(cos_theta)) * std::abs(std::sin(theta_c));
          a_phi_global_edges(j,i) = a_phi_global_edges(j-1,i)
              + bbr_r_faces(j-1,i) * det * (theta_p-theta_m);
        }
      }

      // Calculate cell-centered A_phi based on theta-face-centered B^r
      for (int i = 0; i < sample_n_r; ++i)
      {
        Real r = r_cell(i);
        a_phi_global_cells(0,i) = 0.0;
        for (int j = 1; j < sample_n_theta/2; ++j)
        {
          Real theta_m = theta_cell(j-1);
          Real theta_c = theta_face(j);
          Real theta_p = theta_cell(j);
          Real cos_theta = std::cos(theta_c);
          Real det = (SQR(r) + SQR(a) * SQR(cos_theta)) * std::abs(std::sin(theta_c));
          a_phi_global_cells(j,i) = a_phi_global_cells(j-1,i)
              + bbr_theta_faces(j,i) * det * (theta_p-theta_m);
        }
      }

      // Calculate maximum of vector potential
      Real a_phi_max = 0.0;
      for (int i = 0; i < sample_n_r+1; ++i)
        for (int j = 0; j < sample_n_theta/2+1; ++j)
          a_phi_max = std::max(a_phi_max, a_phi_global_edges(j,i));

      // Floor edge-centered A_phi
      for (int i = 0; i < sample_n_r+1; ++i)
        for (int j = 0; j < sample_n_theta/2+1; ++j)
          a_phi_global_edges(j,i) = std::max(a_phi_global_edges(j,i),
              potential_cutoff*a_phi_max);

      // Floor cell-centered A_phi
      for (int i = 0; i < sample_n_r; ++i)
        for (int j = 0; j < sample_n_theta/2; ++j)
          a_phi_global_cells(j,i) = std::max(a_phi_global_cells(j,i),
              potential_cutoff*a_phi_max);

      // Calculate minimum value of beta over global sample grid
      Real beta_min_actual = std::numeric_limits<Real>::max();
      for (int j = 0; j < sample_n_theta/2; ++j)
        for (int i = 0; i < sample_n_r; ++i)
        {
          Real r_m = r_face(i);
          Real r_c = r_cell(i);
          Real r_p = r_face(i+1);
          Real theta_m = theta_face(j);
          Real theta_c = theta_cell(j);
          Real theta_p = theta_face(j+1);
          if (r_m < r_edge)
            continue;
          Real log_h = LogHAux(r_c, std::sin(theta_c)) - log_h_edge;
          Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
          Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
          if (rho < sample_cutoff)
            continue;
          Real a_phi_cm = 0.5 * (a_phi_global_edges(j,i) + a_phi_global_edges(j,i+1));
          Real a_phi_cp =
              0.5 * (a_phi_global_edges(j+1,i) + a_phi_global_edges(j+1,i+1));
          Real a_phi_mc = 0.5 * (a_phi_global_edges(j,i) + a_phi_global_edges(j+1,i));
          Real a_phi_pc =
              0.5 * (a_phi_global_edges(j,i+1) + a_phi_global_edges(j+1,i+1));
          Real beta;
          bool value_set = CalculateBetaFromA(r_m, r_c, r_p, theta_m, theta_p, theta_c,
              a_phi_cm, a_phi_cp, a_phi_mc, a_phi_pc, &beta);
          if (value_set)
            beta_min_actual = std::min(beta_min_actual, beta);
        }

      // Calculate magnetic field normalization
      if (beta_min < 0.0)
        normalization = 0.0;
      else
        normalization = std::sqrt(beta_min_actual/beta_min);

      // Interpolate edge-centered vector potential onto local grid
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          Real r_c, theta_c, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2f(j),
              pcoord->x3v(kl), &r_c, &theta_c, &phi);
          if (theta_c > PI/2.0)
            theta_c = PI - theta_c;
          int r_index;
          for (r_index = 0; r_index < sample_n_r; ++r_index)
            if (r_face(r_index+1) > r_c)
              break;
          Real r_frac;
          if (r_index == 0)
            r_frac = 0.0;
          else if (r_index == sample_n_r)
          {
            r_index = sample_n_r - 1;
            r_frac = 1.0;
          }
          else
            r_frac = (r_c-r_face(r_index)) / (r_face(r_index+1)-r_face(r_index));
          int theta_index;
          for (theta_index = 0; theta_index < sample_n_theta/2; ++theta_index)
            if (theta_face(theta_index+1) > theta_c)
              break;
          Real theta_frac;
          if (theta_index == 0)
            theta_frac = 0.0;
          else if (theta_index == sample_n_theta/2)
          {
            theta_index = sample_n_theta/2 - 1;
            theta_frac = 1.0;
          }
          else
            theta_frac = (theta_c-theta_face(theta_index))
                / (theta_face(theta_index+1)-theta_face(theta_index));
          Real a_mm = a_phi_global_edges(theta_index,r_index);
          Real a_mp = a_phi_global_edges(theta_index+1,r_index);
          Real a_pm = a_phi_global_edges(theta_index,r_index+1);
          Real a_pp = a_phi_global_edges(theta_index+1,r_index+1);
          a_phi_edges(j,i) = (1.0-r_frac-theta_frac+r_frac*theta_frac) * a_mm
                           + theta_frac*(1.0-r_frac) * a_mp
                           + r_frac*(1.0-theta_frac) * a_pm
                           + r_frac*theta_frac * a_pp;
        }

      // Interpolate cell-centered vector potential onto local grid
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
        {
          Real r_c, theta_c, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
              pcoord->x3v(kl), &r_c, &theta_c, &phi);
          if (theta_c > PI/2.0)
            theta_c = PI - theta_c;
          int r_index;
          for (r_index = 0; r_index < sample_n_r-1; ++r_index)
            if (r_cell(r_index+1) > r_c)
              break;
          Real r_frac;
          if (r_index == 0)
            r_frac = 0.0;
          else if (r_index == sample_n_r-1)
          {
            r_index = sample_n_r - 2;
            r_frac = 1.0;
          }
          else
            r_frac = (r_c-r_cell(r_index)) / (r_cell(r_index+1)-r_cell(r_index));
          int theta_index;
          for (theta_index = 0; theta_index < sample_n_theta/2-1; ++theta_index)
            if (theta_cell(theta_index+1) > theta_c)
              break;
          Real theta_frac;
          if (theta_index == 0)
            theta_frac = 0.0;
          else if (theta_index == sample_n_theta/2-1)
          {
            theta_index = sample_n_theta/2 - 2;
            theta_frac = 1.0;
          }
          else
            theta_frac = (theta_c-theta_cell(theta_index))
                / (theta_cell(theta_index+1)-theta_cell(theta_index));
          Real a_mm = a_phi_global_cells(theta_index,r_index);
          Real a_mp = a_phi_global_cells(theta_index+1,r_index);
          Real a_pm = a_phi_global_cells(theta_index,r_index+1);
          Real a_pp = a_phi_global_cells(theta_index+1,r_index+1);
          a_phi_cells(j,i) = (1.0-r_frac-theta_frac+r_frac*theta_frac) * a_mm
                           + theta_frac*(1.0-r_frac) * a_mp
                           + r_frac*(1.0-theta_frac) * a_pm
                           + r_frac*theta_frac * a_pp;
        }

      // Free global arrays
      r_face.DeleteAthenaArray();
      r_cell.DeleteAthenaArray();
      theta_face.DeleteAthenaArray();
      theta_cell.DeleteAthenaArray();
      a_phi_global_edges.DeleteAthenaArray();
      a_phi_global_cells.DeleteAthenaArray();
      bbr_r_faces.DeleteAthenaArray();
      bbr_theta_faces.DeleteAthenaArray();
    }

    // Handle unknown input
    else
    {
      std::stringstream msg;
      msg << "### FATAL ERROR in Problem Generator\n"
          << "field_config must be \"normal\" or \"renorm\"" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // Set magnetic fields according to vector potential
    // Note: This does very rough differencing for in-face fields on exterior faces of
    //    domain. This should not matter, as these will be identically 0 in nice
    //    coordinate systems or as long as the initial torus is within the domain.
    for (int k = kl; k <= ku+1; ++k)
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          // Set B^1
          if (j != ju+1 and k != ku+1)
          {
            Real r, theta, phi;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2v(j),
                pcoord->x3v(k), &r, &theta, &phi);
            Real r_1, theta_1, phi_1;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2f(j),
                pcoord->x3v(k), &r_1, &theta_1, &phi_1);
            Real r_2, theta_2, phi_2;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2f(j+1),
                pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            Real cos_theta = std::cos(theta);
            Real det = (SQR(r) + SQR(a) * SQR(cos_theta)) * std::abs(std::sin(theta));
            Real bbr =
                1.0/det * (a_phi_edges(j+1,i)-a_phi_edges(j,i)) / (theta_2-theta_1);
            Real a_phi_1, a_phi_2;
            if (i == il)
            {
              a_phi_1 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j+1,i));
              a_phi_2 = a_phi_cells(j,i);
              r_1 = r;
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
                  pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            else if (i == iu+1)
            {
              a_phi_1 = a_phi_cells(j,i-1);
              a_phi_2 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j+1,i));
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i-1), pcoord->x2v(j),
                  pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              r_2 = r;
            }
            else
            {
              a_phi_1 = a_phi_cells(j,i-1);
              a_phi_2 = a_phi_cells(j,i);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i-1), pcoord->x2v(j),
                  pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
                  pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            Real bbtheta = -1.0/det * (a_phi_2-a_phi_1) / (r_2-r_1);
            if (det == 0.0 or (bbr == 0.0 and bbtheta == 0.0))
              pfield->b.x1f(k,j,i) = 0.0;
            else
            {
              Real ut, uphi;
              Real sin_theta = std::sin(theta);
              CalculateVelocityInTorus(r, sin_theta, &ut, &uphi);
              Real sin_sq_theta = SQR(sin_theta);
              Real cos_sq_theta = 1.0 - sin_sq_theta;
              Real rho_sq = SQR(r) + SQR(a) * cos_sq_theta;
              Real bt = -2.0*m*a*r * SQR(sin_theta) / rho_sq * bbr * ut;
              Real br = 1.0/ut * bbr;
              Real btheta = 1.0/ut * bbtheta;
              Real u0, u1, u2, u3;
              pcoord->TransformVectorFace1(ut, 0.0, 0.0, uphi, k, j, i, &u0, &u1, &u2,
                  &u3);
              Real b0, b1, b2, b3;
              pcoord->TransformVectorFace1(bt, br, btheta, 0.0, k, j, i, &b0, &b1, &b2,
                  &b3);
              pfield->b.x1f(k,j,i) = (b1 * u0 - b0 * u1) * normalization;
            }
          }

          // Set B^2
          if (i != iu+1 and k != ku+1)
          {
            Real r, theta, phi;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2f(j),
                pcoord->x3v(k), &r, &theta, &phi);
            Real r_1, theta_1, phi_1;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2f(j),
                pcoord->x3v(k), &r_1, &theta_1, &phi_1);
            Real r_2, theta_2, phi_2;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i+1), pcoord->x2f(j),
                pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            Real cos_theta = std::cos(theta);
            Real det = (SQR(r) + SQR(a) * SQR(cos_theta)) * std::abs(std::sin(theta));
            Real bbtheta = -1.0/det * (a_phi_edges(j,i+1)-a_phi_edges(j,i)) / (r_2-r_1);
            Real a_phi_1, a_phi_2;
            if (j == jl)
            {
              a_phi_1 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j,i+1));
              a_phi_2 = a_phi_cells(j,i);
              theta_1 = theta;
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
                  pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            else if (j == ju+1)
            {
              a_phi_1 = a_phi_cells(j-1,i);
              a_phi_2 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j,i+1));
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j-1),
                  pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              theta_2 = theta;
            }
            else
            {
              a_phi_1 = a_phi_cells(j-1,i);
              a_phi_2 = a_phi_cells(j,i);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j-1),
                  pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
                  pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            Real bbr = 1.0/det * (a_phi_2 - a_phi_1) / (theta_2 - theta_1);
            if (det == 0.0 or (bbr == 0.0 and bbtheta == 0.0))
              pfield->b.x2f(k,j,i) = 0.0;
            else
            {
              Real ut, uphi;
              Real sin_theta = std::sin(theta);
              CalculateVelocityInTorus(r, sin_theta, &ut, &uphi);
              Real sin_sq_theta = SQR(sin_theta);
              Real cos_sq_theta = 1.0 - sin_sq_theta;
              Real rho_sq = SQR(r) + SQR(a) * cos_sq_theta;
              Real bt = -2.0*m*a*r * SQR(sin_theta) / rho_sq * bbr * ut;
              Real br = 1.0/ut * bbr;
              Real btheta = 1.0/ut * bbtheta;
              Real u0, u1, u2, u3;
              pcoord->TransformVectorFace2(ut, 0.0, 0.0, uphi, k, j, i, &u0, &u1, &u2,
                  &u3);
              Real b0, b1, b2, b3;
              pcoord->TransformVectorFace2(bt, br, btheta, 0.0, k, j, i, &b0, &b1, &b2,
                  &b3);
              pfield->b.x2f(k,j,i) = (b2 * u0 - b0 * u2) * normalization;
            }
          }

          // Set B^3
          if (i != iu+1 and j != ju+1)
            pfield->b.x3f(k,j,i) = 0.0;
        }
    a_phi_edges.DeleteAthenaArray();
    a_phi_cells.DeleteAthenaArray();
  }

  // Impose density and pressure floors
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
      {
        Real r, theta, phi;
        pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
            pcoord->x3v(kl), &r, &theta, &phi);
        Real &rho = phydro->w(IDN,k,j,i);
        Real &pgas = phydro->w(IEN,k,j,i);
        rho = std::max(rho, rho_min * std::pow(r, rho_pow));
        pgas = std::max(pgas, (gamma_adi-1.0) * u_min * std::pow(r, u_pow));
        phydro->w1(IDN,k,j,i) = rho;
        phydro->w1(IEN,k,j,i) = pgas;
      }

  // Calculate cell-centered magnetic field
  AthenaArray<Real> bb;
  if (MAGNETIC_FIELDS_ENABLED)
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, il, iu, jl, ju,
        kl, ku);
  else
    bb.NewAthenaArray(3, ku+1, ju+1, iu+1);

  // Initialize conserved values
  if (MAGNETIC_FIELDS_ENABLED)
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, il, iu, jl,
        ju, kl, ku);
  else
  {
    peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, il, iu, jl, ju, kl,
        ku);
    bb.DeleteAthenaArray();
  }

  // Call user work function to set output variables
  UserWorkInLoop();
  return;
}

//--------------------------------------------------------------------------------------

// Function responsible for storing useful quantities for output
// Inputs: (none)
// Outputs: (none)
// Notes:
//   writes to ifov array the following quantities:
//     gamma (normal frame Lorentz factor)
//     p_mag (magnetic pressure)
//     u^mu (coordinate 4-velocity components)
//     b^mu (coordinate 4-magnetic field components)
//   quantities written are specified by NIFOV:
//     0: (nothing)
//     1: gamma
//     2: gamma, p_mag
//     5: gamma, u^mu
//     10: gamma, p_mag, u^mu, b^mu
void MeshBlock::UserWorkInLoop()
{
  // Only proceed if appropriate number of extra output variables specified
  if (not (NIFOV == 1 or NIFOV == 5 or
      (MAGNETIC_FIELDS_ENABLED and (NIFOV == 2 or NIFOV == 10))))
    return;

  // Go through all cells
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
    {
      pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i)
      {
        // Calculate normal frame Lorentz factor
        Real uu1 = phydro->w(IM1,k,j,i);
        Real uu2 = phydro->w(IM2,k,j,i);
        Real uu3 = phydro->w(IM3,k,j,i);
        Real tmp = g(I11,i)*uu1*uu1 + 2.0*g(I12,i)*uu1*uu2 + 2.0*g(I13,i)*uu1*uu3
                 + g(I22,i)*uu2*uu2 + 2.0*g(I23,i)*uu2*uu3
                 + g(I33,i)*uu3*uu3;
        Real gamma = std::sqrt(1.0 + tmp);
        phydro->ifov(0,k,j,i) = gamma;
        if (NIFOV == 1)
          continue;

        // Calculate 4-velocity
        Real alpha = std::sqrt(-1.0/gi(I00,i));
        Real u0 = gamma/alpha;
        Real u1 = uu1 - alpha * gamma * gi(I01,i);
        Real u2 = uu2 - alpha * gamma * gi(I02,i);
        Real u3 = uu3 - alpha * gamma * gi(I03,i);
        if (NIFOV == 5 or NIFOV == 10)
        {
          int offset = (NIFOV == 5) ? 1 : 2;
          phydro->ifov(offset+0,k,j,i) = u0;
          phydro->ifov(offset+1,k,j,i) = u1;
          phydro->ifov(offset+2,k,j,i) = u2;
          phydro->ifov(offset+3,k,j,i) = u3;
        }
        if (NIFOV == 5)
          continue;
        Real u_0, u_1, u_2, u_3;
        pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2, &u_3);

        // Calculate 4-magnetic field
        Real bb1 = pfield->bcc(IB1,k,j,i);
        Real bb2 = pfield->bcc(IB2,k,j,i);
        Real bb3 = pfield->bcc(IB3,k,j,i);
        Real b0 = g(I01,i)*u0*bb1 + g(I02,i)*u0*bb2 + g(I03,i)*u0*bb3
                + g(I11,i)*u1*bb1 + g(I12,i)*u1*bb2 + g(I13,i)*u1*bb3
                + g(I12,i)*u2*bb1 + g(I22,i)*u2*bb2 + g(I23,i)*u2*bb3
                + g(I13,i)*u3*bb1 + g(I23,i)*u3*bb2 + g(I33,i)*u3*bb3;
        Real b1 = (bb1 + b0 * u1) / u0;
        Real b2 = (bb2 + b0 * u2) / u0;
        Real b3 = (bb3 + b0 * u3) / u0;
        if (NIFOV == 10)
        {
          int offset = 6;
          phydro->ifov(offset+0,k,j,i) = b0;
          phydro->ifov(offset+1,k,j,i) = b1;
          phydro->ifov(offset+2,k,j,i) = b2;
          phydro->ifov(offset+3,k,j,i) = b3;
        }
        Real b_0, b_1, b_2, b_3;
        pcoord->LowerVectorCell(b0, b1, b2, b3, k, j, i, &b_0, &b_1, &b_2, &b_3);

        // Calculate magnetic pressure
        Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;
        phydro->ifov(1,k,j,i) = b_sq/2.0;
      }
    }
  return;
}

//--------------------------------------------------------------------------------------

// Fixed boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
// Notes:
//   does nothing
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  return;
}

//--------------------------------------------------------------------------------------

// Inflow boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
void InflowBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // Set hydro variables
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is-NGHOST; i <= is-1; ++i)
      {
        prim(IDN,k,j,i) = prim(IDN,k,j,is);
        prim(IEN,k,j,i) = prim(IEN,k,j,is);
        prim(IM1,k,j,i) = std::min(prim(IM1,k,j,is), 0.0);
        prim(IM2,k,j,i) = prim(IM2,k,j,is);
        prim(IM3,k,j,i) = prim(IM3,k,j,is);
      }

  // Set radial magnetic field
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is-NGHOST; i <= is-1; ++i)
        bb.x1f(k,j,i) = bb.x1f(k,j,is);

  // Set polar magnetic field
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je+1; ++j)
      for (int i = is-NGHOST; i <= is-1; ++i)
        bb.x2f(k,j,i) = bb.x2f(k,j,is);

  // Set azimuthal magnetic field
  for (int k = ks; k <= ke+1; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is-NGHOST; i <= is-1; ++i)
        bb.x3f(k,j,i) = bb.x3f(k,j,is);
  return;
}

//--------------------------------------------------------------------------------------

// Function for calculating angular momentum variable l
// Inputs:
//   r: desired radius of pressure maximum
// Outputs:
//   returned value: l = u^t u_\phi such that pressure maximum occurs at r_peak
// Notes:
//   beware many different definitions of l abound
//     this is *not* -u_phi/u_t
//   Harm has a similar function: lfish_calc() in init.c
//     Harm's function assumes M = 1 and that corotation is desired
//     it is equivalent to this, though seeing this requires much manipulation
//   implements (3.8) from Fishbone & Moncrief 1976, ApJ 207 962
//   assumes corotation
//   TODO: add counterrotation option
//   see CalculateRPeakFromL()
static Real CalculateLFromRPeak(Real r)
{
  Real num = SQR(SQR(r)) + SQR(a*r) - 2.0*m*SQR(a)*r - a*(SQR(r)-SQR(a))*std::sqrt(m*r);
  Real denom = SQR(r) - 3.0*m*r + 2.0*a*std::sqrt(m*r);
  return 1.0/r * std::sqrt(m/r) * num/denom;
}

//--------------------------------------------------------------------------------------

// Function for calculating pressure maximum radius r_peak
// Inputs:
//   l_target: desired u^t u_\phi
// Outputs:
//   returned value: location of pressure maximum given l_target
// Notes:
//   beware many different definitions of l abound
//     this is *not* -u_phi/u_t
//   uses (3.8) from Fishbone & Moncrief 1976, ApJ 207 962
//   assumes corotation
//   TODO: add counterrotation option
//   uses bisection to find r such that formula for l agrees with given value
//   proceeds until either absolute tolerance is met
//   returns best value after max_iterations reached if tolerances not met
//   returns NAN in case of failure (e.g. root not bracketed)
//   see CalculateLFromRPeak()
static Real CalculateRPeakFromL(Real l_target)
{
  // Parameters
  const Real tol_r = 1.0e-10;      // absolute tolerance on abscissa r_peak
  const Real tol_l = 1.0e-10;      // absolute tolerance on ordinate l
  const int max_iterations = 100;  // maximum number of iterations before best res

  // Prepare initial values
  Real r_a = r_min;
  Real r_b = r_max;
  Real r_c = 0.5 * (r_min + r_max);
  Real l_a = CalculateLFromRPeak(r_a);
  Real l_b = CalculateLFromRPeak(r_b);
  Real l_c = CalculateLFromRPeak(r_c);
  if (not ((l_a < l_target and l_b > l_target) or (l_a > l_target and l_b < l_target)))
    return NAN;

  // Find root
  for (int n = 0; n < max_iterations; ++n)
  {
    if (std::abs(r_b-r_a) <= 2.0*tol_r or std::abs(l_c-l_target) <= tol_l)
      break;
    if ((l_a < l_target and l_c < l_target) or (l_a > l_target and l_c > l_target))
    {
      r_a = r_c;
      l_a = l_c;
    }
    else
    {
      r_b = r_c;
      l_b = l_c;
    }
    r_c = 0.5 * (r_min + r_max);
    l_c = CalculateLFromRPeak(r_c);
  }
  return r_c;
}

//--------------------------------------------------------------------------------------

// Function for helping to calculate enthalpy
// Inputs:
//   r: radial Boyer-Lindquist coordinate
//   sin_sq_theta: square of sine of polar Boyer-Lindquist coordinate
// Outputs:
//   returned value: log(h)
// Notes:
//   enthalpy defined here as h = p_gas/rho
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//   implements first half of (FM 3.6)
static Real LogHAux(Real r, Real sin_theta)
{
  Real sin_sq_theta = SQR(sin_theta);
  Real cos_sq_theta = 1.0 - sin_sq_theta;
  Real delta = SQR(r) - 2.0*m*r + SQR(a);                // \Delta
  Real sigma = SQR(r) + SQR(a)*cos_sq_theta;             // \Sigma
  Real aa = SQR(SQR(r)+SQR(a))
      - delta*SQR(a)*sin_sq_theta;                       // A
  Real exp_2nu = sigma * delta / aa;                     // \exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sin_sq_theta;             // \exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;                 // \exp(-2\chi) (cf. FM 2.15)
  Real omega = 2.0*m*a*r/aa;                             // \omega (FM 3.5)
  Real var_a = std::sqrt(1.0 + 4.0*SQR(l)*exp_neg2chi);
  Real var_b = 0.5 * std::log((1.0+var_a)
      / (sigma*delta/aa));
  Real var_c = -0.5 * var_a;
  Real var_d = -l * omega;
  return var_b + var_c + var_d;                          // (FM 3.4)
}

//--------------------------------------------------------------------------------------

// Function for computing 4-velocity components at a given position inside torus
// Inputs:
//   r: Boyer-Lindquist r
//   sin_theta: sine of Boyer-Lindquist theta
// Outputs:
//   pu0: u^t set (Boyer-Lindquist coordinates)
//   pu3: u^\phi set (Boyer-Lindquist coordinates)
// Notes:
//   The formula for u^3 as a function of u_{(\phi)} is tedious to derive,
//       but this matches the formula used in Harm (init.c).
static void CalculateVelocityInTorus(Real r, Real sin_theta, Real *pu0, Real *pu3)
{
  Real sin_sq_theta = SQR(sin_theta);
  Real cos_sq_theta = 1.0 - sin_sq_theta;
  Real delta = SQR(r) - 2.0*m*r + SQR(a);            // \Delta
  Real sigma = SQR(r) + SQR(a)*cos_sq_theta;         // \Sigma
  Real aa = SQR(SQR(r)+SQR(a))
      - delta*SQR(a)*sin_sq_theta;                   // A
  Real exp_2nu = sigma * delta / aa;                 // \exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sin_sq_theta;         // \exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;             // \exp(-2\chi) (cf. FM 2.15)
  Real u_phi_proj_a = 1.0 + 4.0*SQR(l)*exp_neg2chi;
  Real u_phi_proj_b = -1.0
      + std::sqrt(u_phi_proj_a);
  Real u_phi_proj = std::sqrt(0.5 * u_phi_proj_b);   // (FM 3.3)
  Real u3_a = (1.0+SQR(u_phi_proj))
      / (aa*sigma*delta);
  Real u3_b = 2.0*m*a*r * std::sqrt(u3_a);
  Real u3_c = std::sqrt(sigma/aa) / sin_theta;
  Real u3 = u3_b + u3_c * u_phi_proj;
  Real g_00 = -(1.0 - 2.0*m*r/sigma);
  Real g_03 = -2.0*m*a*r/sigma * sin_sq_theta;
  Real g_33 = (sigma + (1.0 + 2.0*m*r/sigma)
      * SQR(a) * sin_sq_theta) * sin_sq_theta;
  Real u0_a = (SQR(g_03) - g_00*g_33) * SQR(u3);
  Real u0_b = std::sqrt(u0_a - g_00);
  Real u0 = -1.0/g_00 * (g_03*u3 + u0_b);
  *pu0 = u0;
  *pu3 = u3;
  return;
}

//--------------------------------------------------------------------------------------

// Function for finding approximate minimum value of plasma beta expected
// Inputs: (none)
// Outputs:
//   returned value: minimum beta found by sampling grid covering whole domain
// Notes:
//   constructs grid over entire mesh, not just block
//   grid is not necessarily the same as used for the problem proper
//   calculation is done entirely in Boyer-Lindquist coordinates
static Real CalculateBetaMin()
{
  // Prepare container to hold minimum
  Real beta_min_actual = std::numeric_limits<Real>::max();

  // Go through sample grid in theta
  for (int j = 0; j < sample_n_theta; ++j)
  {
    // Calculate theta values
    Real theta_m = theta_min + static_cast<Real>(j)/static_cast<Real>(sample_n_theta)
        * (theta_max-theta_min);
    Real theta_p = theta_min + static_cast<Real>(j+1)/static_cast<Real>(sample_n_theta)
        * (theta_max-theta_min);
    Real theta_c = 0.5 * (theta_m + theta_p);
    Real sin_theta_c = std::sin(theta_c);
    Real cos_theta_c = std::cos(theta_c);

    // Go through sample grid in r
    Real r_m, r_p, delta_r;
    for (int i = 0; i < sample_n_r; ++i)
    {
      // Calculate r values
      if (i == 0)
      {
        r_m = r_min;
        Real ratio_power = 1.0;
        Real ratio_sum = 1.0;
        for (int ii = 1; ii < sample_n_r; ++ii)
        {
          ratio_power *= sample_r_rat;
          ratio_sum += ratio_power;
        }
        delta_r = (r_max-r_min) / ratio_sum;
      }
      else
      {
        r_m = r_p;
        delta_r *= sample_r_rat;
      }
      r_p = r_m + delta_r;
      Real r_c = 0.5 * (r_m + r_p);

      // Calculate beta
      Real beta;
      bool value_set = CalculateBeta(r_m, r_c, r_p, theta_m, theta_c, theta_p, &beta);
      if (value_set)
        beta_min_actual = std::min(beta_min_actual, beta);
    }
  }
  return beta_min_actual;
}

//--------------------------------------------------------------------------------------

// Function for calculating beta from four nearby points
// Inputs:
//   r_m,r_c,r_p: inner, center, and outer radii
//   theta_m,theta_c,theta_p: upper, center, and lower polar angles
// Outputs:
//   pbeta: value set to plasma beta at cell center
//   returned value: true if pbeta points to meaningful number (inside torus)
// Notes:
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
static bool CalculateBeta(Real r_m, Real r_c, Real r_p, Real theta_m, Real theta_c,
    Real theta_p, Real *pbeta)
{
  // Calculate trigonometric functions of theta
  Real sin_theta_m = std::sin(theta_m);
  Real sin_theta_p = std::sin(theta_p);
  Real sin_theta_c = std::sin(theta_c);
  Real cos_theta_c = std::cos(theta_c);

  // Determine if we are in the torus (FM 3.6)
  if (r_m < r_edge)
    return false;
  Real log_h_cc = LogHAux(r_c, sin_theta_c) - log_h_edge;
  Real log_h_cm = LogHAux(r_c, sin_theta_m) - log_h_edge;
  Real log_h_cp = LogHAux(r_c, sin_theta_p) - log_h_edge;
  Real log_h_mc = LogHAux(r_m, sin_theta_c) - log_h_edge;
  Real log_h_pc = LogHAux(r_p, sin_theta_c) - log_h_edge;
  if (log_h_cc < 0.0 or log_h_cm < 0.0 or log_h_cp < 0.0 or log_h_mc < 0.0
      or log_h_pc < 0.0)
    return false;

  // Calculate primitives depending on location
  Real pgas_over_rho_cc = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_cc)-1.0);
  Real pgas_over_rho_cm = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_cm)-1.0);
  Real pgas_over_rho_cp = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_cp)-1.0);
  Real pgas_over_rho_mc = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_mc)-1.0);
  Real pgas_over_rho_pc = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_pc)-1.0);
  Real rho_cc = std::pow(pgas_over_rho_cc/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  if (rho_cc < sample_cutoff)
    return false;
  Real rho_cm = std::pow(pgas_over_rho_cm/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real rho_cp = std::pow(pgas_over_rho_cp/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real rho_mc = std::pow(pgas_over_rho_mc/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real rho_pc = std::pow(pgas_over_rho_pc/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real pgas = pgas_over_rho_cc * rho_cc;
  Real u1 = 0.0;
  Real u2 = 0.0;
  Real u0, u3;
  CalculateVelocityInTorus(r_c, sin_theta_c, &u0, &u3);

  // Calculate phi-component of vector potential
  Real rho_cm_cutoff = std::max(rho_cm-potential_cutoff, 0.0);
  Real rho_cp_cutoff = std::max(rho_cp-potential_cutoff, 0.0);
  Real rho_mc_cutoff = std::max(rho_mc-potential_cutoff, 0.0);
  Real rho_pc_cutoff = std::max(rho_pc-potential_cutoff, 0.0);
  Real a_cm = std::pow(r_c,potential_r_pow) * std::pow(rho_cm_cutoff,potential_rho_pow);
  Real a_cp = std::pow(r_c,potential_r_pow) * std::pow(rho_cp_cutoff,potential_rho_pow);
  Real a_mc = std::pow(r_m,potential_r_pow) * std::pow(rho_mc_cutoff,potential_rho_pow);
  Real a_pc = std::pow(r_p,potential_r_pow) * std::pow(rho_pc_cutoff,potential_rho_pow);
  if (a_cm == 0.0 or a_cp == 0.0 or a_mc == 0.0 or a_pc == 0.0)
    return false;

  // Calculate cell-centered 3-magnetic field
  Real det = (SQR(r_c) + SQR(a) * SQR(cos_theta_c)) * std::abs(sin_theta_c);
  Real bb1 = 1.0/det * (a_cp-a_cm) / (theta_p-theta_m);
  Real bb2 = -1.0/det * (a_pc-a_mc) / (r_p-r_m);
  Real bb3 = 0.0;

  // Calculate beta
  Real pmag = CalculateMagneticPressure(bb1, bb2, bb3, r_c, sin_theta_c, cos_theta_c);
  *pbeta = pgas/pmag;
  return true;
}

//--------------------------------------------------------------------------------------

// Function for calculating beta given vector potential
// Inputs:
//   r_m,r_c,r_p: inner, center, and outer radii
//   theta_m,theta_c,theta_p: upper, center, and lower polar angles
//   a_cm,a_cp,a_mc,a_pc: A_phi offset by theta (down,up) and r (down,up)
// Outputs:
//   pbeta: value set to plasma beta at cell center
//   returned value: true if pbeta points to meaningful number (inside torus)
// Notes:
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
static bool CalculateBetaFromA(Real r_m, Real r_c, Real r_p, Real theta_m, Real theta_c,
    Real theta_p, Real a_cm, Real a_cp, Real a_mc, Real a_pc, Real *pbeta)
{
  // Calculate trigonometric functions of theta
  Real sin_theta_c = std::sin(theta_c);
  Real cos_theta_c = std::cos(theta_c);

  // Determine if we are in the torus (FM 3.6)
  if (r_m < r_edge)
    return false;
  Real log_h = LogHAux(r_c, sin_theta_c) - log_h_edge;
  if (log_h < 0.0)
    return false;

  // Calculate primitives
  Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
  Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real pgas = pgas_over_rho * rho;

  // Check A_phi
  if (a_cm == 0.0 or a_cp == 0.0 or a_mc == 0.0 or a_pc == 0.0)
    return false;

  // Calculate 3-magnetic field
  Real det = (SQR(r_c) + SQR(a) * SQR(cos_theta_c)) * std::abs(sin_theta_c);
  Real bb1 = 1.0/det * (a_cp-a_cm) / (theta_p-theta_m);
  Real bb2 = -1.0/det * (a_pc-a_mc) / (r_p-r_m);
  Real bb3 = 0.0;

  // Calculate beta
  Real pmag = CalculateMagneticPressure(bb1, bb2, bb3, r_c, sin_theta_c, cos_theta_c);
  *pbeta = pgas/pmag;
  return true;
}

//--------------------------------------------------------------------------------------

// Function to calculate b^lambda b_lambda / 2
// Inputs:
//   bb1,bb2,bb3: components of 3-magnetic field in Boyer-Lindquist coordinates
//   r: Boyer-Lindquist radius
//   sin_theta,cos_theta: sine and cosine of Boyer-Lindquist theta
// Outputs:
//   returned value: magnetic pressure
static Real CalculateMagneticPressure(Real bb1, Real bb2, Real bb3, Real r,
    Real sin_theta, Real cos_theta)
{
  // Calculate Boyer-Lindquist metric
  Real delta = SQR(r) - 2.0*m*r + SQR(a);
  Real sigma = SQR(r) + SQR(a) * SQR(cos_theta);
  Real g_00 = -(1.0 - 2.0*m*r/sigma);
  Real g_01 = 0.0;
  Real g_02 = 0.0;
  Real g_03 = -2.0*m*a*r/sigma * SQR(sin_theta);
  Real g_11 = sigma/delta;
  Real g_12 = 0.0;
  Real g_13 = 0.0;
  Real g_22 = sigma;
  Real g_23 = 0.0;
  Real g_33 = (SQR(r) + SQR(a) + 2.0*m*SQR(a)*r/sigma * SQR(sin_theta))
      * SQR(sin_theta);
  Real g_10 = g_01;
  Real g_20 = g_02;
  Real g_21 = g_12;
  Real g_30 = g_03;
  Real g_31 = g_13;
  Real g_32 = g_23;

  // Calculate 4-velocity
  Real u1 = 0.0;
  Real u2 = 0.0;
  Real u0, u3;
  CalculateVelocityInTorus(r, sin_theta, &u0, &u3);

  // Calculate 4-magnetic field
  Real b0 = bb1 * (g_10*u0 + g_11*u1 + g_12*u2 + g_13*u3)
          + bb2 * (g_20*u0 + g_21*u1 + g_22*u2 + g_23*u3)
          + bb3 * (g_30*u0 + g_31*u1 + g_32*u2 + g_33*u3);
  Real b1 = 1.0/u0 * (bb1 + b0 * u1);
  Real b2 = 1.0/u0 * (bb2 + b0 * u2);
  Real b3 = 1.0/u0 * (bb3 + b0 * u3);

  // Calculate magnetic pressure
  Real b_sq = g_00*b0*b0 + g_01*b0*b1 + g_02*b0*b2 + g_03*b0*b3
            + g_10*b1*b0 + g_11*b1*b1 + g_12*b1*b2 + g_13*b1*b3
            + g_20*b2*b0 + g_21*b2*b1 + g_22*b2*b2 + g_23*b2*b3
            + g_30*b3*b0 + g_31*b3*b1 + g_32*b3*b2 + g_33*b3*b3;
  return 0.5*b_sq;
}
