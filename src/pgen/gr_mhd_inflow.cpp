//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_mhd_inflow.cpp
//  \brief Problem generator for magnetized equatorial inflow around Kerr black hole.

// C++ headers
#include <cmath>      // cos, sin, sqrt
#include <fstream>    // ifstream
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string, c_str()

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // macros, enums, FaceField
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro

// Configuration checking
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif

// Declarations
void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
static void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                         Real *ptheta, Real *pphi);
static void TransformVector(Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, Real r,
                     Real theta, Real phi, Real *pa0, Real *pa1, Real *pa2, Real *pa3);
static void CalculateFromTable(Real r, Real theta, Real *prho, Real *put, Real *pur,
                               Real *puphi, Real *pbt, Real *pbr, Real *pbphi);

// Global variables
static Real m;                           // mass M of black hole
static Real a;                           // spin of black hole (0 <= a < M)
static Real temperature;                 // temperature pgas/rho
static AthenaArray<Real> interp_values;  // table for analytic solution
static int num_lines;                    // number of lines in table

//----------------------------------------------------------------------------------------
// Function for initializing global mesh properties
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Read temperature
  temperature = pin->GetReal("problem", "temperature");

  // Allocate array for interpolation points
  num_lines = pin->GetInteger("problem", "num_data_lines");
  interp_values.NewAthenaArray(6, num_lines);

  // Read interpolation data from file
  std::string filename = pin->GetString("problem", "data_file");
  std::ifstream file(filename.c_str());
  if (not file.is_open()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator\n"
        << "file " << filename << " cannot be opened" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int n = 0; n < num_lines; ++n) {
    Real r, rho, ur, uphi, bbr, bbphi;
    file >> r >> rho >> ur >> uphi >> bbr >> bbphi;
    interp_values(0,n) = r;
    interp_values(1,n) = rho;
    interp_values(2,n) = ur;
    interp_values(3,n) = uphi;
    interp_values(4,n) = bbr;
    interp_values(5,n) = bbphi;
  }

  // Enroll fixed outer boundary function
  EnrollUserBoundaryFunction(OUTER_X1, FixedBoundary);
  return;
}

//----------------------------------------------------------------------------------------
// Function for cleaning up global mesh properties
// Inputs:
//   pin: parameters (unused)
// Outputs: (none)

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // Free interpolation table
  interp_values.DeleteAthenaArray();
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: pointer to runtime inputs
// Outputs: (none)
// Notes:
//   initializes equatorial inflow
//   see Gammie 1999, ApJ 522 L57
//       Gammie, McKinney, & Toth 2003, ApJ 589 444

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Get mass and spin of black hole
  m = pcoord->GetMass();
  a = pcoord->GetSpin();

  // Prepare variables to hold results from multiple-return functions
  Real r, theta, phi;   // Boyer-Lindquist (BL) coordinates
  Real rho;             // density
  Real ut, ur, uphi;    // BL u^\mu
  Real bt, br, bphi;    // BL b^\mu
  Real u0, u1, u2, u3;  // preferred coordinates u^\mu
  Real b0, b1, b2, b3;  // preferred coordinates b^\mu

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {

    // Initialize radial field components
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is-NGHOST; i <= ie+NGHOST+1; ++i) {
          Real x1 = pcoord->x1f(i);
          Real x2 = pcoord->x2v(j);
          Real x3 = pcoord->x3v(k);
          GetBoyerLindquistCoordinates(x1, x2, x3, &r, &theta, &phi);
          CalculateFromTable(r, theta, &rho, &ut, &ur, &uphi, &bt, &br, &bphi);
          TransformVector(ut, ur, 0.0, uphi, r, theta, phi, &u0, &u1, &u2, &u3);
          TransformVector(bt, br, 0.0, bphi, r, theta, phi, &b0, &b1, &b2, &b3);
          pfield->b.x1f(k,j,i) = b1*u0 - b0*u1;
        }
      }
    }

    // Initialize poloidal field components
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {
        for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2f(j);
          Real x3 = pcoord->x3v(k);
          GetBoyerLindquistCoordinates(x1, x2, x3, &r, &theta, &phi);
          CalculateFromTable(r, theta, &rho, &ut, &ur, &uphi, &bt, &br, &bphi);
          TransformVector(ut, ur, 0.0, uphi, r, theta, phi, &u0, &u1, &u2, &u3);
          TransformVector(bt, br, 0.0, bphi, r, theta, phi, &b0, &b1, &b2, &b3);
          pfield->b.x2f(k,j,i) = b2*u0 - b0*u2;
        }
      }
    }

    // Initialize azimuthal field components
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);
          Real x3 = pcoord->x3f(k);
          GetBoyerLindquistCoordinates(x1, x2, x3, &r, &theta, &phi);
          CalculateFromTable(r, theta, &rho, &ut, &ur, &uphi, &bt, &br, &bphi);
          TransformVector(ut, ur, 0.0, uphi, r, theta, phi, &u0, &u1, &u2, &u3);
          TransformVector(bt, br, 0.0, bphi, r, theta, phi, &b0, &b1, &b2, &b3);
          pfield->b.x3f(k,j,i) = b3*u0 - b0*u3;
        }
      }
    }
  }

  // Calculate cell-centered magnetic field
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ke+1, je+1, ie+NGHOST+1);
  if (MAGNETIC_FIELDS_ENABLED) {
    pfield->CalculateCellCenteredField(pfield->b, bb, pcoord, is-NGHOST, ie+NGHOST, js,
        je, ks, ke);
  }

  // Initialize primitive values
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC,ie+NGHOST+1);
  gi.NewAthenaArray(NMETRIC,ie+NGHOST+1);
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->CellMetric(k, j, is-NGHOST, ie+NGHOST, g, gi);
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        GetBoyerLindquistCoordinates(x1, x2, x3, &r, &theta, &phi);
        CalculateFromTable(r, theta, &rho, &ut, &ur, &uphi, &bt, &br, &bphi);
        TransformVector(ut, ur, 0.0, uphi, r, theta, phi, &u0, &u1, &u2, &u3);
        TransformVector(bt, br, 0.0, bphi, r, theta, phi, &b0, &b1, &b2, &b3);
        Real pgas = temperature * rho;
        Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
        Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
        Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = uu1;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uu2;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uu3;
      }
    }
  }
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Initialize conserved values
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is-NGHOST, ie+NGHOST, js,
      je, ks, ke);
  bb.DeleteAthenaArray();
  return;
}

//----------------------------------------------------------------------------------------
// Fixed boundary condition
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   time,dt: current time and timestep of simulation
//   is,ie,js,je,ks,ke: indices demarkating active region
// Outputs:
//   prim: primitives set in ghost zones
//   bb: face-centered magnetic field set in ghost zones
// Notes:
//   does nothing

void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  return;
}

//----------------------------------------------------------------------------------------
// Function for returning corresponding Boyer-Lindquist coordinates of point
// Inputs:
//   x1,x2,x3: global coordinates to be converted
// Outputs:
//   pr,ptheta,pphi: variables pointed to set to Boyer-Lindquist coordinates
// Notes:
//   conversion is trivial in all currently implemented coordinate systems

static void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3, Real *pr,
                                         Real *ptheta, Real *pphi) {
  if (COORDINATE_SYSTEM == "schwarzschild" or COORDINATE_SYSTEM == "kerr-schild") {
    *pr = x1;
    *ptheta = x2;
    *pphi = x3;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming 4-vector from Boyer-Lindquist to desired coordinates
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   r,theta,phi: Boyer-Lindquist coordinates of point
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in desired coordinates
// Notes:
//   Schwarzschild coordinates match Boyer-Lindquist when a = 0

static void TransformVector(Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, Real r,
                     Real theta, Real phi, Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  if (COORDINATE_SYSTEM == "schwarzschild") {
    *pa0 = a0_bl;
    *pa1 = a1_bl;
    *pa2 = a2_bl;
    *pa3 = a3_bl;
  } else if (COORDINATE_SYSTEM == "kerr-schild") {
    Real delta = SQR(r) - 2.0*m*r + SQR(a);
    *pa0 = a0_bl + 2.0*m*r/delta * a1_bl;
    *pa1 = a1_bl;
    *pa2 = a2_bl;
    *pa3 = a3_bl + a/delta * a1_bl;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating quantities based on table
// Inputs:
//   r,theta: Boyer-Lindquist radial and polar coordinates
// Outputs:
//   prho: value set to interpolated density
//   put,pur,puphi: values set to interpolated u^\mu in Boyer-Lindquist coordinates
//   pbt,pbr,pbphi: values set to interpolated b^\mu in Boyer-Lindquist coordinates

static void CalculateFromTable(Real r, Real theta, Real *prho, Real *put, Real *pur,
                               Real *puphi, Real *pbt, Real *pbr, Real *pbphi) {
  // Find location in interpolation table
  int n;
  Real fraction = 0.0;
  if (r < interp_values(0,0)) {
    n = 0;
    fraction = 0.0;
  } else if (r >= interp_values(0,num_lines-1)) {
    n = num_lines - 1;
    fraction = 1.0;
  } else {
    for (n = 0; n < num_lines-1; ++n) {
      if (r < interp_values(0,n+1)) {
        fraction = (r-interp_values(0,n)) / (interp_values(0,n+1)-interp_values(0,n));
        break;
      }
    }
  }

  // Interpolate to location based on table
  *prho = (1.0-fraction)*interp_values(1,n) + fraction*interp_values(1,n+1);
  Real ur = (1.0-fraction)*interp_values(2,n) + fraction*interp_values(2,n+1);
  Real uphi = (1.0-fraction)*interp_values(3,n) + fraction*interp_values(3,n+1);
  Real bbr = (1.0-fraction)*interp_values(4,n) + fraction*interp_values(4,n+1);
  Real bbphi = (1.0-fraction)*interp_values(5,n) + fraction*interp_values(5,n+1);

  // Calculate velocity
  Real sin = std::sin(theta);
  Real cos = std::cos(theta);
  Real delta = SQR(r) - 2.0*m*r + SQR(a);
  Real sigma = SQR(r) + SQR(a) * SQR(cos);
  Real g_tt = -(1.0 - 2.0*m*r/sigma);
  Real g_tphi = -2.0*m*a*r/sigma * SQR(sin);
  Real g_rr = sigma/delta;
  Real g_phiphi = (SQR(r) + SQR(a) + 2.0*m*SQR(a)*r/sigma * SQR(sin)) * SQR(sin);
  Real var_a = g_tt;
  Real var_b = 2.0*g_tphi*uphi;
  Real var_c = g_rr*SQR(ur) + g_phiphi*SQR(uphi) + 1.0;
  Real ut;
  if (var_a == 0.0) {
    ut = -var_c/var_b;
  } else {
    Real a1 = var_b/var_a;
    Real a0 = var_c/var_a;
    Real s2 = SQR(a1) - 4.0*a0;
    Real s = (s2 < 0.0) ? 0.0 : std::sqrt(s2);
    ut = (s2 >= 0.0 and a1 >= 0.0) ? -2.0*a0/(a1+s) : (-a1+s)/2.0;
  }
  *put = ut;
  *pur = ur;
  *puphi = uphi;

  // Calculate covariant magnetic field
  Real bt = g_tphi*ut*bbphi + g_rr*ur*bbr + g_phiphi*uphi*bbphi;
  *pbt = bt;
  *pbr = 1.0/ut * (bbr + bt * ur);
  *pbphi = 1.0/ut * (bbphi + bt * uphi);
  return;
}
