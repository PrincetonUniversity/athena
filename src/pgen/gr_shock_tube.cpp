//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_shock_tube.cpp
//  \brief Problem generator for shock tubes in special and general relativity.

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro

// Configuration checking
#if not RELATIVISTIC_DYNAMICS
#error "This problem generator must be used with relativity"
#endif

// Declarations
static void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
    Real *px, Real *py, Real *pz);
static void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3);

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters
// Outputs: (none)
// Notes:
//   sets conserved variables according to input primitives
//   assigns fields based on cell-center positions, rather than interface positions
//     this helps shock tube 2 from Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141
//     otherwise the middle interface would go to left variables, creating a
//         particularly troublesome jump leading to NaN's

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Read and set ratio of specific heats
  Real gamma_adi = peos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read and check shock direction and position
  int shock_dir = pin->GetInteger("problem", "shock_dir");
  Real shock_pos = pin->GetReal("problem", "xshock");
  Real min_bound, max_bound;
  std::stringstream msg;
  switch (shock_dir) {
    case 1:
      min_bound = pmy_mesh->mesh_size.x1min;
      max_bound = pmy_mesh->mesh_size.x1max;
      break;
    case 2:
      min_bound = pmy_mesh->mesh_size.x2min;
      max_bound = pmy_mesh->mesh_size.x2max;
      break;
    case 3:
      min_bound = pmy_mesh->mesh_size.x3min;
      max_bound = pmy_mesh->mesh_size.x3max;
      break;
    default:
      msg << "### FATAL ERROR in Problem Generator\n"
          << "shock_dir=" << shock_dir << " must be either 1, 2, or 3" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  if (shock_pos < min_bound || shock_pos > max_bound) {
    msg << "### FATAL ERROR in Problem Generator\n"
        << "xshock=" << shock_pos << " lies outside x" << shock_dir
            << " domain for shkdir=" << shock_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Read left state
  Real rho_left = pin->GetReal("problem", "dl");
  Real pgas_left = pin->GetReal("problem", "pl");
  Real vx_left = pin->GetReal("problem", "ul");
  Real vy_left = pin->GetReal("problem", "vl");
  Real vz_left = pin->GetReal("problem", "wl");
  Real bbx_left = 0.0, bby_left = 0.0, bbz_left = 0.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    bbx_left = pin->GetReal("problem", "bxl");
    bby_left = pin->GetReal("problem", "byl");
    bbz_left = pin->GetReal("problem", "bzl");
  }

  // Read right state
  Real rho_right = pin->GetReal("problem", "dr");
  Real pgas_right = pin->GetReal("problem", "pr");
  Real vx_right = pin->GetReal("problem", "ur");
  Real vy_right = pin->GetReal("problem", "vr");
  Real vz_right = pin->GetReal("problem", "wr");
  Real bbx_right = 0.0, bby_right = 0.0, bbz_right = 0.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    bbx_right = pin->GetReal("problem", "bxr");
    bby_right = pin->GetReal("problem", "byr");
    bbz_right = pin->GetReal("problem", "bzr");
  }

  // Prepare auxiliary arrays
  AthenaArray<Real> bb, g, gi;
  bb.NewAthenaArray(3, ke+1, je+1, ie+1);
  if (GENERAL_RELATIVITY) {
    g.NewAthenaArray(NMETRIC, ie+1);
    gi.NewAthenaArray(NMETRIC, ie+1);
  }

  // Initialize hydro variables
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      #if GENERAL_RELATIVITY
      {
        pcoord->CellMetric(k, j, is, ie, g, gi);
      }
      #endif  // GENERAL_RELATIVITY
      for (int i = is; i <= ie; ++i) {

        // Determine which variables to use
        Real rho = rho_right;
        Real pgas = pgas_right;
        Real vx = vx_right;
        Real vy = vy_right;
        Real vz = vz_right;
        Real bbx = bbx_right;
        Real bby = bby_right;
        Real bbz = bbz_right;
        bool left_side = false;
        switch(shock_dir) {
          case 1:
            left_side = pcoord->x1v(i) < shock_pos;
            break;
          case 2:
            left_side = pcoord->x2v(j) < shock_pos;
            break;
          case 3:
            left_side = pcoord->x3v(k) < shock_pos;
            break;
        }
        if (left_side) {
          rho = rho_left;
          pgas = pgas_left;
          vx = vx_left;
          vy = vy_left;
          vz = vz_left;
          bbx = bbx_left;
          bby = bby_left;
          bbz = bbz_left;
        }

        // Construct 4-vectors
        Real ut = std::sqrt(1.0 / (1.0 - (SQR(vx)+SQR(vy)+SQR(vz))));
        Real ux = ut * vx;
        Real uy = ut * vy;
        Real uz = ut * vz;
        Real bt = bbx*ux + bby*uy + bbz*uz;
        Real bx = (bbx + bt * ux) / ut;
        Real by = (bby + bt * uy) / ut;
        Real bz = (bbz + bt * uz) / ut;

        // Transform 4-vectors
        Real u0, u1, u2, u3;
        Real b0, b1, b2, b3;
        #if GENERAL_RELATIVITY
        {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);
          Real x3 = pcoord->x3v(k);
          Real t, x, y, z;
          GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
          TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
          TransformVector(bt, bx, by, bz, x, y, z, &b0, &b1, &b2, &b3);
        }
        #else  // SR
        {
          u0 = ut;
          u1 = ux;
          u2 = uy;
          u3 = uz;
          b0 = bt;
          b1 = bx;
          b2 = by;
          b3 = bz;
        }
        #endif  // GENERAL_RELATIVITY

        // Set primitives
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        if (GENERAL_RELATIVITY) {
          Real uu1 = u1 - gi(I01,i)/gi(I00,i) * u0;
          Real uu2 = u2 - gi(I02,i)/gi(I00,i) * u0;
          Real uu3 = u3 - gi(I03,i)/gi(I00,i) * u0;
          phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = uu1;
          phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = uu2;
          phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = uu3;
        } else {
          phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = u1 / u0;
          phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = u2 / u0;
          phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = u3 / u0;
        }

        // Set magnetic fields
        bb(IB1,k,j,i) = b1 * u0 - b0 * u1;
        bb(IB2,k,j,i) = b2 * u0 - b0 * u2;
        bb(IB3,k,j,i) = b3 * u0 - b0 * u3;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is, ie, js, je, ks, ke);

  // Delete auxiliary arrays
  bb.DeleteAthenaArray();
  if (GENERAL_RELATIVITY) {
    g.DeleteAthenaArray();
    gi.DeleteAthenaArray();
  }

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je+1; ++j) {
        for (int i = is; i <= ie+1; ++i) {

          // Determine which variables to use
          Real vx = vx_right;
          Real vy = vy_right;
          Real vz = vz_right;
          Real bbx = bbx_right;
          Real bby = bby_right;
          Real bbz = bbz_right;
          bool left_side = false;
          switch(shock_dir) {
            case 1:
              left_side = pcoord->x1v(i) < shock_pos;
              break;
            case 2:
              left_side = pcoord->x2v(j) < shock_pos;
              break;
            case 3:
              left_side = pcoord->x3v(k) < shock_pos;
              break;
          }
          if (left_side) {
            vx = vx_left;
            vy = vy_left;
            vz = vz_left;
            bbx = bbx_left;
            bby = bby_left;
            bbz = bbz_left;
          }

          // Construct 4-vectors
          Real ut = std::sqrt(1.0 / (1.0 - (SQR(vx)+SQR(vy)+SQR(vz))));
          Real ux = ut * vx;
          Real uy = ut * vy;
          Real uz = ut * vz;
          Real bt = bbx*ux + bby*uy + bbz*uz;
          Real bx = (bbx + bt * ux) / ut;
          Real by = (bby + bt * uy) / ut;
          Real bz = (bbz + bt * uz) / ut;

          // Set magnetic fields
          Real u0, u1, u2, u3;
          Real b0, b1, b2, b3;
          #if GENERAL_RELATIVITY
          {
            if (j != je+1 && k != ke+1) {
              Real x1 = pcoord->x1f(i);
              Real x2 = pcoord->x2v(j);
              Real x3 = pcoord->x3v(k);
              Real t, x, y, z;
              GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
              TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
              TransformVector(bt, bx, by, bz, x, y, z, &b0, &b1, &b2, &b3);
              pfield->b.x1f(k,j,i) = b1 * u0 - b0 * u1;
            }
            if (i != ie+1 && k != ke+1) {
              Real x1 = pcoord->x1v(i);
              Real x2 = pcoord->x2f(j);
              Real x3 = pcoord->x3v(k);
              Real t, x, y, z;
              GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
              TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
              TransformVector(bt, bx, by, bz, x, y, z, &b0, &b1, &b2, &b3);
              pfield->b.x2f(k,j,i) = b2 * u0 - b0 * u2;
            }
            if (i != ie+1 && j != je+1) {
              Real x1 = pcoord->x1v(i);
              Real x2 = pcoord->x2v(j);
              Real x3 = pcoord->x3f(k);
              Real t, x, y, z;
              GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
              TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
              TransformVector(bt, bx, by, bz, x, y, z, &b0, &b1, &b2, &b3);
              pfield->b.x3f(k,j,i) = b3 * u0 - b0 * u3;
            }
          }
          #else  // SR
          {
            if (j != je+1 && k != ke+1) {
              pfield->b.x1f(k,j,i) = bbx;
            }
            if (i != ie+1 && k != ke+1) {
              pfield->b.x2f(k,j,i) = bby;
            }
            if (i != ie+1 && j != je+1) {
              pfield->b.x3f(k,j,i) = bbz;
            }
          }
          #endif  // GENERAL_RELATIVITY
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for returning corresponding Minkowski coordinates of point
// Inputs:
//   x0,x1,x2,x3: global coordinates to be converted
// Outputs:
//   pt,px,py,pz: variables pointed to set to Minkowski coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

static void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
    Real *px, Real *py, Real *pz) {
  if (COORDINATE_SYSTEM == "minkowski") {
    *pt = x0;
    *px = x1;
    *py = x2;
    *pz = x3;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming 4-vector from Minkowski to desired coordinates
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   x,y,z: Minkowski coordinates of point
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in desired coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

static void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  if (COORDINATE_SYSTEM == "minkowski") {
    *pa0 = at;
    *pa1 = ax;
    *pa2 = ay;
    *pa3 = az;
  }
  return;
}
