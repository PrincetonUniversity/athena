// General relativistic Fishbone-Moncrief torus generator

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // exp(), log(), NAN, pow(), sin(), sqrt()
#include <iostream>   // endl
#include <limits>     // numeric_limits::max()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str(), string

// Athena headers
#include "../athena.hpp"                   // macros, enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues, FaceField
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"

// Declarations
static Real calculate_l_from_r_peak(Real r);
static Real calculate_r_peak_from_l(Real l_target, Real r_min, Real r_max);
static Real log_h_aux(Real r, Real sin_theta);
static void calculate_velocity_in_torus(Real r, Real sin_theta, Real *pu0, Real *pu3);
static Real calculate_beta_min(Real r_min, Real r_max, Real theta_min, Real theta_max,
    int n_r, int n_theta, Real r_ratio);
static bool CalculateBeta(Real r_m, Real r_c, Real r_p, Real theta_m, Real theta_c,
    Real theta_p, Real *pbb1, Real *pbb2, Real *pbb3, Real *pbeta);
static Real CalculateMagneticPressure(Real bb1, Real bb2, Real bb3, Real r,
    Real sin_theta, Real cos_theta);
static Real CalculatePotentialOffset(Real r, Real theta, Real r_min, Real r_max,
    int n_r, int n_theta, Real r_rat);

// Global variables
static Real m, a;                            // black hole parameters
static Real gamma_adi, k_adi;                // hydro parameters
static Real r_edge, r_peak, l, rho_max;      // fixed torus parameters
static Real log_h_edge, log_h_peak;          // calculated torus parameters
static Real pgas_over_rho_peak, rho_peak;    // more calculated torus parameters
static Real rho_min, rho_pow, u_min, u_pow;  // background parameters
static Real potential_cutoff;                // sets region of torus to magnetize
static Real beta_min;                        // min ratio of gas to magnetic pressure
static std::string field_shape;              // type of field to use

// Function for initializing global mesh properties
void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  return;
}

// Function for cleaning up global mesh properties
void Mesh::TerminateUserMeshProperties(void)
{
  return;
}

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
  gamma_adi = phydro->peos->GetGamma();

  // Read other properties
  rho_min = pin->GetReal("hydro", "rho_min");
  rho_pow = pin->GetReal("hydro", "rho_pow");
  u_min = pin->GetReal("hydro", "u_min");
  u_pow = pin->GetReal("hydro", "u_pow");
  k_adi = pin->GetReal("problem", "k_adi");
  r_edge = pin->GetReal("problem", "r_edge");
  r_peak = pin->GetReal("problem", "r_peak");
  l = pin->GetReal("problem", "l");
  rho_max = pin->GetReal("problem", "rho_max");
  potential_cutoff = pin->GetReal("problem", "cutoff");
  beta_min = pin->GetReal("problem", "beta_min");

  // Reset whichever of l,r_peak is not specified
  if (r_peak >= 0.0)
    l = calculate_l_from_r_peak(r_peak);
  else
    r_peak = calculate_r_peak_from_l(l, pin->GetReal("mesh", "x1min"),
        pin->GetReal("mesh", "x1max"));

  // Prepare scratch arrays
  AthenaArray<bool> in_torus;
  in_torus.NewAthenaArray(ju+1, iu+1);
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, iu+1);
  gi.NewAthenaArray(NMETRIC, iu+1);

  // Initialize primitive values
  log_h_edge = log_h_aux(r_edge, 1.0);
  log_h_peak = log_h_aux(r_peak, 1.0) - log_h_edge;
  pgas_over_rho_peak = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h_peak)-1.0);
  rho_peak = std::pow(pgas_over_rho_peak/k_adi, 1.0/(gamma_adi-1.0));
  for (int j = jl; j <= ju; ++j)
  {
    pcoord->CellMetric(kl, j, il, iu, g, gi);
    for (int i = il; i <= iu; ++i)
    {
      // Get Boyer-Lindquist coordinates of cell
      Real r, theta, phi;
      pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
          pcoord->x2v(j), pcoord->x3v(kl), &r, &theta, &phi);
      Real sin_theta = std::sin(theta);

      // Determine if we are in the torus
      Real log_h;
      in_torus(j,i) = false;
      if (r >= r_edge)
      {
        log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
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
        calculate_velocity_in_torus(r, sin_theta, &u0, &u3);
        Real u0_pref, u1_pref, u2_pref, u3_pref;
        pcoord->TransformVectorCell(u0, 0.0, 0.0, u3, kl, j, i,
            &u0_pref, &u1_pref, &u2_pref, &u3_pref);
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

      // Set primitive values
      for (int k = kl; k <= ku; ++k)
      {
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IEN,k,j,i) = phydro->w1(IEN,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = uu1;
        phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = uu2;
        phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = uu3;
      }
    }
  }

  // Initialize magnetic fields
  if (MAGNETIC_FIELDS_ENABLED)
  {
    // Determine field shape
    field_shape = pin->GetOrAddString("problem", "field_shape", "standard");
    if (not (field_shape == "standard" or field_shape == "offset"))
    {
      std::stringstream msg;
      msg << "### FATAL ERROR in Problem Generator" << std::endl
          << "must select valid field_shape" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // Prepare 2D arrays of vector potential values
    AthenaArray<Real> a_phi_cells, a_phi_edges;
    a_phi_cells.NewAthenaArray(ju+1, iu+1);
    a_phi_edges.NewAthenaArray(ju+2, iu+2);

    // Set vector potential in 2D slice in standard case (proportional to rho)
    if (field_shape == "standard")
    {
      // Set cell-centered values
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
              pcoord->x2v(j), pcoord->x3v(kl), &r, &theta, &phi);
          Real sin_theta = std::sin(theta);
          if (r >= r_edge)
          {
            Real log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
            if (log_h >= 0.0)
            {
              Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
              Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
              a_phi_cells(j,i) = std::max(rho-potential_cutoff, 0.0);
            }
          }
        }

      // Set edge-centered values
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i),
              pcoord->x2f(j), pcoord->x3v(kl), &r, &theta, &phi);
          if (r >= r_edge)
          {
            Real sin_theta = std::sin(theta);
            Real log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
            if (log_h >= 0.0)
            {
              Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
              Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
              a_phi_edges(j,i) = std::max(rho-potential_cutoff, 0.0);
            }
          }
        }
    }

    // Extract parameters for sampling grid
    int sample_n_r = pin->GetInteger("problem", "sample_n_r");
    int sample_n_theta = pin->GetInteger("problem", "sample_n_theta");
    sample_n_theta += sample_n_theta % 2;
    Real sample_r_rat = pin->GetReal("problem", "sample_r_rat");
    Real x1_min = pin->GetReal("mesh", "x1min");
    Real x1_max = pin->GetReal("mesh", "x1max");
    Real x2_min = pin->GetReal("mesh", "x2min");
    Real x2_max = pin->GetReal("mesh", "x2max");
    Real r1, r2, r3, r4, theta1, theta2, theta3, theta4, temp;
    pcoord->GetBoyerLindquistCoordinates(x1_min, x2_min, pcoord->x3v(kl), &r1, &theta1,
        &temp);
    pcoord->GetBoyerLindquistCoordinates(x1_max, x2_min, pcoord->x3v(kl), &r2, &theta2,
        &temp);
    pcoord->GetBoyerLindquistCoordinates(x1_min, x2_max, pcoord->x3v(kl), &r3, &theta3,
        &temp);
    pcoord->GetBoyerLindquistCoordinates(x1_max, x2_max, pcoord->x3v(kl), &r4, &theta4,
        &temp);
    Real r_min = std::min(std::min(r1, r2), std::min(r3, r4));
    Real r_max = std::max(std::max(r1, r2), std::max(r3, r4));
    Real theta_min = std::min(std::min(theta1, theta2), std::min(theta3, theta4));
    Real theta_max = std::max(std::max(theta1, theta2), std::max(theta3, theta4));

    // Set vector potential in 2D slice in offset case (based on r^5*rho^2)
    if (field_shape == "offset")
    {
      // Set cell-centered values
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i), pcoord->x2v(j),
              pcoord->x3v(kl), &r, &theta, &phi);
          a_phi_cells(j,i) = CalculatePotentialOffset(r, theta, r_min, r_max,
              sample_n_r, sample_n_theta, sample_r_rat);
        }

      // Set edge-centered values
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          Real r, theta, phi;
          pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i), pcoord->x2f(j),
              pcoord->x3v(kl), &r, &theta, &phi);
          a_phi_edges(j,i) = CalculatePotentialOffset(r, theta, r_min, r_max,
              sample_n_r, sample_n_theta, sample_r_rat);
        }
    }

    // Calculate magnetic field normalization
    r_peak = calculate_r_peak_from_l(l, r_min, r_max);
    Real normalization;
    if (beta_min < 0.0)
      normalization = 0.0;
    else
    {
      Real beta_min_actual = calculate_beta_min(r_min, r_max, theta_min, theta_max,
          sample_n_r, sample_n_theta, sample_r_rat);
      normalization = std::sqrt(beta_min_actual/beta_min);
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
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i),
                pcoord->x2v(j), pcoord->x3v(k), &r, &theta, &phi);
            Real r_1, theta_1, phi_1;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i),
                pcoord->x2f(j), pcoord->x3v(k), &r_1, &theta_1, &phi_1);
            Real r_2, theta_2, phi_2;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i),
                pcoord->x2f(j+1), pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            Real bbr = -(a_phi_edges(j+1,i) - a_phi_edges(j,i)) / (theta_2 - theta_1);
            Real a_phi_1, a_phi_2;
            if (i == il)
            {
              a_phi_1 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j+1,i));
              a_phi_2 = a_phi_cells(j,i);
              r_1 = r;
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                  pcoord->x2v(j), pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            else if (i == iu+1)
            {
              a_phi_1 = a_phi_cells(j,i-1);
              a_phi_2 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j+1,i));
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i-1),
                  pcoord->x2v(j), pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              r_2 = r;
            }
            else
            {
              a_phi_1 = a_phi_cells(j,i-1);
              a_phi_2 = a_phi_cells(j,i);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i-1),
                  pcoord->x2v(j), pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                  pcoord->x2v(j), pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            Real bbtheta = (a_phi_2 - a_phi_1) / (r_2 - r_1);
            if (bbr == 0.0 and bbtheta == 0.0)
              pfield->b.x1f(k,j,i) = 0.0;
            else
            {
              Real ut, uphi;
              Real sin_theta = std::sin(theta);
              calculate_velocity_in_torus(r, sin_theta, &ut, &uphi);
              Real sin_sq_theta = SQR(sin_theta);
              Real cos_sq_theta = 1.0 - sin_sq_theta;
              Real rho_sq = SQR(r) + SQR(a) * cos_sq_theta;
              Real bt = -2.0*m*a*r * SQR(sin_theta) / rho_sq * bbr * ut;
              Real br = 1.0/ut * bbr;
              Real btheta = 1.0/ut * bbtheta;
              Real u0, u1, u2, u3;
              pcoord->TransformVectorFace1(ut, 0.0, 0.0, uphi, k, j, i,
                  &u0, &u1, &u2, &u3);
              Real b0, b1, b2, b3;
              pcoord->TransformVectorFace1(bt, br, btheta, 0.0, k, j, i,
                  &b0, &b1, &b2, &b3);
              pfield->b.x1f(k,j,i) = (b1 * u0 - b0 * u1) * normalization;
            }
          }

          // Set B^2
          if (i != iu+1 and k != ku+1)
          {
            Real r, theta, phi;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                pcoord->x2f(j), pcoord->x3v(k), &r, &theta, &phi);
            Real r_1, theta_1, phi_1;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i),
                pcoord->x2f(j), pcoord->x3v(k), &r_1, &theta_1, &phi_1);
            Real r_2, theta_2, phi_2;
            pcoord->GetBoyerLindquistCoordinates(pcoord->x1f(i+1),
                pcoord->x2f(j), pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            Real bbtheta = (a_phi_edges(j,i+1) - a_phi_edges(j,i)) / (r_2 - r_1);
            Real a_phi_1, a_phi_2;
            if (j == jl)
            {
              a_phi_1 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j,i+1));
              a_phi_2 = a_phi_cells(j,i);
              theta_1 = theta;
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                  pcoord->x2v(j), pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            else if (j == ju+1)
            {
              a_phi_1 = a_phi_cells(j-1,i);
              a_phi_2 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j,i+1));
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                  pcoord->x2v(j-1), pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              theta_2 = theta;
            }
            else
            {
              a_phi_1 = a_phi_cells(j-1,i);
              a_phi_2 = a_phi_cells(j,i);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                  pcoord->x2v(j-1), pcoord->x3v(k), &r_1, &theta_1, &phi_1);
              pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
                  pcoord->x2v(j), pcoord->x3v(k), &r_2, &theta_2, &phi_2);
            }
            Real bbr = -(a_phi_2 - a_phi_1) / (theta_2 - theta_1);
            if (bbr == 0.0 and bbtheta == 0.0)
              pfield->b.x2f(k,j,i) = 0.0;
            else
            {
              Real ut, uphi;
              Real sin_theta = std::sin(theta);
              calculate_velocity_in_torus(r, sin_theta, &ut, &uphi);
              Real sin_sq_theta = SQR(sin_theta);
              Real cos_sq_theta = 1.0 - sin_sq_theta;
              Real rho_sq = SQR(r) + SQR(a) * cos_sq_theta;
              Real bt = -2.0*m*a*r * SQR(sin_theta) / rho_sq * bbr * ut;
              Real br = 1.0/ut * bbr;
              Real btheta = 1.0/ut * bbtheta;
              Real u0, u1, u2, u3;
              pcoord->TransformVectorFace2(ut, 0.0, 0.0, uphi, k, j, i,
                  &u0, &u1, &u2, &u3);
              Real b0, b1, b2, b3;
              pcoord->TransformVectorFace2(bt, br, btheta, 0.0, k, j, i,
                  &b0, &b1, &b2, &b3);
              pfield->b.x2f(k,j,i) = (b2 * u0 - b0 * u2) * normalization;
            }
          }

          // Set B^3
          if (i != iu+1 and j != ju+1)
            pfield->b.x3f(k,j,i) = 0.0;
        }
    a_phi_cells.DeleteAthenaArray();
    a_phi_edges.DeleteAthenaArray();
  }

  // Impose density and pressure floors
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
      {
        Real r, theta, phi;
        pcoord->GetBoyerLindquistCoordinates(pcoord->x1v(i),
            pcoord->x2v(j), pcoord->x3v(kl), &r, &theta, &phi);
        Real &rho = phydro->w(IDN,k,j,i);
        Real &pgas = phydro->w(IEN,k,j,i);
        rho = std::max(rho, rho_min * std::pow(r, rho_pow));
        pgas = std::max(pgas, (gamma_adi-1.0) * u_min * std::pow(r, u_pow));
        phydro->w1(IDN,k,j,i) = rho;
        phydro->w1(IEN,k,j,i) = pgas;
      }

  // Calculate cell-centered magnetic field
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ku+1, ju+1, iu+1);
  if (MAGNETIC_FIELDS_ENABLED)
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
        {
          // Extract face-centered magnetic field
          const Real &bbf1m = pfield->b.x1f(k,j,i);
          const Real &bbf1p = pfield->b.x1f(k,j,i+1);
          const Real &bbf2m = pfield->b.x2f(k,j,i);
          const Real &bbf2p = pfield->b.x2f(k,j+1,i);
          const Real &bbf3m = pfield->b.x3f(k,j,i);
          const Real &bbf3p = pfield->b.x3f(k+1,j,i);

          // Calculate cell-centered magnetic field
          Real tmp = (pcoord->x1v(i) - pcoord->x1f(i)) / pcoord->dx1f(i);
          bb(IB1,k,j,i) = (1.0-tmp) * bbf1m + tmp * bbf1p;
          tmp = (pcoord->x2v(j) - pcoord->x2f(j)) / pcoord->dx2f(j);
          bb(IB2,k,j,i) = (1.0-tmp) * bbf2m + tmp * bbf2p;
          tmp = (pcoord->x3v(k) - pcoord->x3f(k)) / pcoord->dx3f(k);
          bb(IB3,k,j,i) = (1.0-tmp) * bbf3m + tmp * bbf3p;
       }

  // Initialize conserved values
  phydro->peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord,
                                          il, iu, jl, ju, kl, ku);

  // Free scratch arrays
  in_torus.DeleteAthenaArray();
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();
  bb.DeleteAthenaArray();
  return;
}


// User-defined work function called every time step
void MeshBlock::UserWorkInLoop(void)
{
  return;
}


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
//   see calculate_r_peak_from_l()
static Real calculate_l_from_r_peak(Real r)
{
  Real num = SQR(SQR(r)) + SQR(a*r) - 2.0*m*SQR(a)*r - a*(SQR(r)-SQR(a))*std::sqrt(m*r);
  Real denom = SQR(r) - 3.0*m*r + 2.0*a*std::sqrt(m*r);
  return 1.0/r * std::sqrt(m/r) * num/denom;
}

// Function for calculating pressure maximum radius r_peak
// Inputs:
//   l_target: desired u^t u_\phi
//   r_min,r_max: bounds on r_peak that bracket root
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
//   see calculate_l_from_r_peak()
static Real calculate_r_peak_from_l(Real l_target, Real r_min, Real r_max)
{
  // Parameters
  const Real tol_r = 1.0e-10;      // absolute tolerance on abscissa r_peak
  const Real tol_l = 1.0e-10;      // absolute tolerance on ordinate l
  const int max_iterations = 100;  // maximum number of iterations before best res

  // Prepare initial values
  Real a = r_min;
  Real b = r_max;
  Real c = 0.5 * (r_min + r_max);
  Real l_a = calculate_l_from_r_peak(a);
  Real l_b = calculate_l_from_r_peak(b);
  Real l_c = calculate_l_from_r_peak(c);
  if (not ((l_a < l_target and l_b > l_target) or (l_a > l_target and l_b < l_target)))
    return NAN;

  // Find root
  for (int n = 0; n < max_iterations; ++n)
  {
    if (std::abs(b-a) <= 2.0*tol_r or std::abs(l_c-l_target) <= tol_l)
      break;
    if ((l_a < l_target and l_c < l_target) or (l_a > l_target and l_c > l_target))
    {
      a = c;
      l_a = l_c;
    }
    else
    {
      b = c;
      l_b = l_c;
    }
    c = 0.5 * (r_min + r_max);
    l_c = calculate_l_from_r_peak(c);
  }
  return c;
}

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
static Real log_h_aux(Real r, Real sin_theta)
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
static void calculate_velocity_in_torus(Real r, Real sin_theta, Real *pu0, Real *pu3)
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

// Function for finding approximate minimum value of plasma beta expected
// Inputs:
//   r_min,r_max: bounds on Boyer-Lindquist radial extent of grid
//   theta_min,theta_max: bounds on Boyer-Lindquist polar extent of grid
//   n_r,n_theta: numbers of cells in radial and polar directions in sample grid
//   r_ratio: geometric ratio of adjacent radial widths in sample grid
// Outputs:
//   returned value: minimum beta found by sampling grid covering whole domain
// Notes:
//   constructs grid over entire mesh, not just block
//   grid is not necessarily the same as used for the problem proper
//   calculation is done entirely in Boyer-Lindquist coordinates
static Real calculate_beta_min(Real r_min, Real r_max, Real theta_min, Real theta_max,
    int n_r, int n_theta, Real r_ratio)
{
  // Prepare container to hold minimum
  Real beta_min_actual = std::numeric_limits<Real>::max();

  // Go through sample grid in theta
  for (int j = 0; j < n_theta; ++j)
  {
    // Calculate theta values
    Real theta_m = theta_min
        + static_cast<Real>(j)/static_cast<Real>(n_theta) * (theta_max-theta_min);
    Real theta_p = theta_min
        + static_cast<Real>(j+1)/static_cast<Real>(n_theta) * (theta_max-theta_min);
    Real theta_c = 0.5 * (theta_m + theta_p);
    Real sin_theta_c = std::sin(theta_c);
    Real cos_theta_c = std::cos(theta_c);

    // Go through sample grid in r
    for (int i = 0; i < n_r; ++i)
    {
      // Calculate r values
      Real r_m, r_p, delta_r;
      if (i == 0)
      {
        r_m = r_min;
        Real ratio_power = 1.0;
        Real ratio_sum = 1.0;
        for (int ii = 1; ii < n_r; ++ii)
        {
          ratio_power *= r_ratio;
          ratio_sum += ratio_power;
        }
        delta_r = (r_max-r_min) / ratio_sum;
      }
      else
      {
        r_m = r_p;
        delta_r *= r_ratio;
      }
      r_p = r_m + delta_r;
      Real r_c = 0.5 * (r_m + r_p);

      // Calculate beta in standard case
      if (field_shape == "standard")
      {
        Real beta, bbr, bbtheta, bbphi;
        bool values_set = CalculateBeta(r_m, r_c, r_p, theta_m, theta_c, theta_p, &bbr,
            &bbtheta, &bbphi, &beta);
        if (values_set)
          beta_min_actual = std::min(beta_min_actual, beta);
      }

      // Calculate beta in offset case
      if (field_shape == "offset")
      {
        // Calculate true vector potential at sample cell edges
        Real a_cm = CalculatePotentialOffset(r_c, theta_m, r_min, r_max, n_r, n_theta,
            r_ratio);
        Real a_cp = CalculatePotentialOffset(r_c, theta_p, r_min, r_max, n_r, n_theta,
            r_ratio);
        Real a_mc = CalculatePotentialOffset(r_m, theta_c, r_min, r_max, n_r, n_theta,
            r_ratio);
        Real a_pc = CalculatePotentialOffset(r_p, theta_c, r_min, r_max, n_r, n_theta,
            r_ratio);
        if (a_cm == 0.0 or a_cp == 0.0 or a_mc == 0.0 or a_pc == 0.0)
          continue;

        // Calculate magnetic field and pressure
        Real bbr = -(a_cp-a_cm) / (theta_p-theta_m);
        Real bbtheta = (a_pc-a_mc) / (r_p-r_m);
        Real pmag = CalculateMagneticPressure(bbr, bbtheta, 0.0, r_c, sin_theta_c,
            cos_theta_c);

        // Calculate gas pressure
        Real log_h = log_h_aux(r_c, sin_theta_c) - log_h_edge;
        Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1.0);
        Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
        Real pgas = pgas_over_rho * rho;

        // Reset minimum beta
        beta_min_actual = std::min(beta_min_actual, pgas/pmag);
      }
    }
  }
  return beta_min_actual;
}

// Function for calculating beta from four nearby points
// Inputs:
//   r_m,r_c,r_p: inner, center, and outer radii
//   theta_m,theta_c,theta_p: upper, center, and lower polar angles
// Outputs:
//   pbb1,pbb2,pbb3: values set to 3-magnetic field at cell center
//   pbeta: value set to plasma beta at cell center
// Notes:
//   in field_shape=="offset" case, only the original r^5*rho^2 potential is considered
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
static bool CalculateBeta(Real r_m, Real r_c, Real r_p, Real theta_m, Real theta_c,
    Real theta_p, Real *pbb1, Real *pbb2, Real *pbb3, Real *pbeta)
{
  // Calculate trigonometric functions of theta
  Real sin_theta_m = std::sin(theta_m);
  Real sin_theta_p = std::sin(theta_p);
  Real sin_theta_c = std::sin(theta_c);
  Real cos_theta_c = std::cos(theta_c);

  // Determine if we are in the torus (FM 3.6)
  if (r_m < r_edge)
    return false;
  Real log_h_cc = log_h_aux(r_c, sin_theta_c) - log_h_edge;
  Real log_h_cm = log_h_aux(r_c, sin_theta_m) - log_h_edge;
  Real log_h_cp = log_h_aux(r_c, sin_theta_p) - log_h_edge;
  Real log_h_mc = log_h_aux(r_m, sin_theta_c) - log_h_edge;
  Real log_h_pc = log_h_aux(r_p, sin_theta_c) - log_h_edge;
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
  Real rho_cm = std::pow(pgas_over_rho_cm/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real rho_cp = std::pow(pgas_over_rho_cp/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real rho_mc = std::pow(pgas_over_rho_mc/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real rho_pc = std::pow(pgas_over_rho_pc/k_adi, 1.0/(gamma_adi-1.0)) / rho_peak;
  Real pgas = pgas_over_rho_cc * rho_cc;
  Real u1 = 0.0;
  Real u2 = 0.0;
  Real u0, u3;
  calculate_velocity_in_torus(r_c, sin_theta_c, &u0, &u3);

  // Calculate phi-component of vector potential
  Real a_cm, a_cp, a_mc, a_pc;
  if (field_shape == "standard")
  {
    a_cm = std::max(rho_cm-potential_cutoff, 0.0);
    a_cp = std::max(rho_cp-potential_cutoff, 0.0);
    a_mc = std::max(rho_mc-potential_cutoff, 0.0);
    a_pc = std::max(rho_pc-potential_cutoff, 0.0);
  }
  if (field_shape == "offset")
  {
    a_cm = std::pow(r_c, 5) * SQR(rho_cm);
    a_cp = std::pow(r_c, 5) * SQR(rho_cp);
    a_mc = std::pow(r_m, 5) * SQR(rho_mc);
    a_pc = std::pow(r_p, 5) * SQR(rho_pc);
  }
  if (a_cm == 0.0 or a_cp == 0.0 or a_mc == 0.0 or a_pc == 0.0)
    return false;

  // Calculate cell-centered 3-magnetic field
  *pbb1 = -(a_cp-a_cm) / (theta_p-theta_m);
  *pbb2 = (a_pc-a_mc) / (r_p-r_m);
  *pbb3 = 0.0;

  // Calculate beta
  Real pmag = CalculateMagneticPressure(*pbb1, *pbb2, *pbb3, r_c, sin_theta_c,
      cos_theta_c);
  *pbeta = pgas / pmag;
  return true;
}

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
  calculate_velocity_in_torus(r, sin_theta, &u0, &u3);

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

// Function for calculating A_phi in offset potential
// Inputs:
//   r,theta: Boyer-Lindquist coordinates of point
//   r_min,r_max: inner and outer limits of sample grid
//   n_r,n_theta: numbers of cells in radial and polar directions in sample grid
//   r_ratio: geometric ratio of adjacent radial widths in sample grid
// Outputs:
//   returned value: A_phi at (r,theta)
// Notes:
//   procedure adapted from Tchekhovskoy, Narayan, & McKinney 2011, MNRAS 418 L79:
//     A'_phi = r^5 rho^2
//     B^r = -\partial(A'_phi) / \partial(theta)
//     B^r readjusted locally to correspond to beta = 1
//     A_phi = \int_0^theta B^r sin(theta') d theta'
//   integral done from nearest pole rather than North pole
//   uses sample grid over entire mesh, not just block
//   grid is not necessarily the same as used for the problem proper
static Real CalculatePotentialOffset(Real r, Real theta, Real r_min, Real r_max,
    int n_r, int n_theta, Real r_rat)
{
  // Assume no field beyond poles
  if (theta <= 0.0 or theta >= PI)
    return 0.0;
  Real delta_theta = PI / n_theta;

  // Find r-boundaries in sample grid
  Real r_m, r_p, delta_r;
  for (int i = 0; i < n_r; ++i)
  {
    if (i == 0)
    {
      r_m = r_min;
      Real ratio_power = 1.0;
      Real ratio_sum = 1.0;
      for (int ii = 1; ii < n_r; ++ii)
      {
        ratio_power *= r_rat;
        ratio_sum += ratio_power;
      }
      delta_r = (r_max-r_min) / ratio_sum;
    }
    else
    {
      r_m = r_p;
      delta_r *= r_rat;
    }
    r_p = r_m + delta_r;
    if (r_p >= r)
      break;
  }

  // Go through all cells at same radius from nearest pole to here
  Real a_phi = 0.0;
  for (int j = 0; j < n_theta/2; ++j)
  {
    // Calculate + and - theta values
    Real theta_m, theta_p;
    Real cell_fraction = 1.0;
    if (theta <= PI/2.0)
    {
      theta_m = static_cast<Real>(j)/static_cast<Real>(n_theta/2) * PI/2.0;
      theta_p = static_cast<Real>(j+1)/static_cast<Real>(n_theta/2) * PI/2.0;
      if (theta_m >= theta)
        break;
      if (theta_p > theta)
        cell_fraction = (theta-theta_m) / (theta_p-theta_m);
    }
    else
    {
      theta_p = PI - static_cast<Real>(j)/static_cast<Real>(n_theta/2) * PI/2.0;
      theta_m = PI - static_cast<Real>(j+1)/static_cast<Real>(n_theta/2) * PI/2.0;
      if (theta_p <= theta)
        break;
      if (theta_m < theta)
        cell_fraction = (theta-theta_m) / (theta_p-theta_m);
    }
    Real theta_c = 0.5*(theta_m+theta_p);

    // Calculate B^r and beta to add to integral
    Real beta, bbr, bbtheta, bbphi;
    bool values_set = CalculateBeta(r_m, r, r_p, theta_m, theta_c, theta_p, &bbr,
        &bbtheta, &bbphi, &beta);
    if (not values_set)
      continue;
    bbr *= std::sqrt(beta);
    a_phi += cell_fraction * bbr * std::sin(theta_c) * delta_theta;
  }
  return a_phi;
}
