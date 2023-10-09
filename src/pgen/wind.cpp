//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wind.cpp
//! \brief Problem generator for spherical wind problem.  Works in 2D Cylindrical,
//!        and 3D polar spherical coordinates. It is partially setup to  run in 
//!        cartesian coordinates, but this requires adding rotating system source terms 
//!        that will work in cartesian coordinates.


// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

// add suns gravity as source term
void SunGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

// user defined time step function
Real SetMinimumTimeStep(MeshBlock *pmb);

// User-defined boundary conditions
// radial in cylindrical coords
void CMEInnerX1(MeshBlock *pmb, Coordinates *pco,
                 AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void CMEOuterX1(MeshBlock *pmb, Coordinates *pco,
                 AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// helper functions
Real calc_v_radial_inner(ParameterInput *pin);
Real calc_GM_sun(ParameterInput *pin);
Real calc_b_radial_inner(ParameterInput *pin);
Real calc_b_azimuthal_inner(ParameterInput *pin, Real _b1, Real _v1_inner);
Real calc_n_inner(ParameterInput *pin, Real v1_inner);
Real calc_energy_inner(ParameterInput *pin, Real n_inner, Real gamma);
Real sign_radial_mag_field(Real theta, Real phi);

int RefinementConditionPressure(MeshBlock *pmb);
int RefinementConditionDensityJump(MeshBlock *pmb);

// namespaced variables, i.e. globals
// should be turned into a class with setters and getters
namespace {
  bool enable_pressure_refine;
  Real press_threshold;
  Real v1_inner, v2_inner, v3_inner;
  Real n_inner, e_inner;
  Real inner_radius, gamma_param;
  // Real r_measure;
  // Real m_p;
  Real b1, b2, b3;
  Real x_0, y_0, z_0;
  Real GM_norm;
  bool sun_gravity;
  bool enable_user_dt;
  Real user_dt;
  Real wave_rad;
  Real wave_mult;
  Real omega_sun;
  
  Real CME_start;
  Real CME_duration;
  
  Real CME_density;
  Real CME_velocity;
  Real CME_energy;

  Real CME_phi_0;
  Real CME_phi_extent;
  
} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Enroll user-defined physical source terms
  // external gravitational potential
  sun_gravity = pin->GetOrAddBoolean("problem", "enable_gravity",false);
  if (sun_gravity) {
    EnrollUserExplicitSourceFunction(SunGravity);
  }

  // enroll user defined timestep
  enable_user_dt = pin->GetOrAddBoolean("problem", "enable_user_dt",false);
  if (enable_user_dt) {
    EnrollUserTimeStepFunction(SetMinimumTimeStep);
  }

  if (adaptive) {
    enable_pressure_refine = pin->GetOrAddBoolean("problem", "enable_pressure_refine",false);
    if (enable_pressure_refine) {
      EnrollUserRefinementCondition(RefinementConditionPressure);
      press_threshold = pin->GetReal("problem","press_thresh");
    } else { // density ratio refinement
      EnrollUserRefinementCondition(RefinementConditionDensityJump);
    }
  }

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, CMEInnerX1);
  }

  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, CMEOuterX1);
  }

  ///////////////////
#if 1
  ///gamma_param = peos->GetGamma();
  gamma_param = pin->GetOrAddReal("problem", "gamma", 1.666666666667);
  // would need to change this if not cylindrical
  // get minimum x1 grid location get inner radius
  Real x1_min = pin->GetReal("mesh", "x1min");
  inner_radius = pin->GetOrAddReal("problem", "radius", x1_min);
  // r_measure = pin->GetOrAddReal("problem", "r_measure", 1.0);
  // m_p = pin->GetOrAddReal("problem", "m_p", 1.67E-27);
  // Real omega = pin->GetOrAddReal("problem", "omega_sun", 2.87e-6);
  //assume orbital_advection is enabled
  // Real omega_sun = pin->GetReal("orbital_advection", "Omega0");
  // TODO need to make it clear the coordinate systems and inputs
  v1_inner = calc_v_radial_inner(pin);
  v2_inner = 0.0; // create calc for v2 and v3 when orbital advection is off
  v3_inner = 0.0;
  n_inner = calc_n_inner(pin, v1_inner);
  e_inner = calc_energy_inner(pin, n_inner, gamma_param);

  // defaults
  b1 = calc_b_radial_inner(pin);
  b2 = calc_b_azimuthal_inner(pin, b1, v1_inner);
  b3 = 0.0;

  //output for each meshblock
  std::cout << "r_inner = " << inner_radius << std::endl;
  std::cout << "n_inner = " << n_inner << std::endl;
  std::cout << "v1_inner = " << v1_inner << std::endl;
  std::cout << "e_inner = " << e_inner << std::endl;
  std::cout << "b1 = " << b1 << std::endl;
  std::cout << "b3 = " << b2 << std::endl;
  std::cout << "--------------------------" << std::endl;
  

  // CME data
  Real t_o = pin->GetReal("problem", "t_o");
  Real vo = pin->GetReal("problem", "vo");

  //CME_start = pin->GetReal("problem", "CME_timestart_hrs")*3600/t_o;
  CME_start = pin->GetOrAddReal("problem", "CME_timestart_hrs",1e20)*3600/t_o;
  CME_duration = pin->GetOrAddReal("problem", "CME_timespan_hrs",0.)*3600/t_o;

  CME_density = pin->GetOrAddReal("problem", "CME_density",1.);
  CME_velocity = pin->GetOrAddReal("problem", "CME_velocity",0.)*1000/vo;
  CME_energy = e_inner*2;

  CME_phi_0 = (PI/180.0)*pin->GetOrAddReal("problem", "CME_phi_0",100.);
  CME_phi_extent = (PI/180.0)*pin->GetOrAddReal("problem", "CME_phi_extent",0.);

  // Real b0; //, angle;
  // if (MAGNETIC_FIELDS_ENABLED) {
  //   b0 = pin->GetReal("problem", "b0");
  //   angle = (PI/180.0)*pin->GetReal("problem", "angle");
  // }

  // get coordinates of center of sun, and convert to Cartesian if necessary
  // should all be zero for now
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x_0 = x1_0;
    y_0 = x2_0;
    z_0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x_0 = x1_0*std::cos(x2_0);
    y_0 = x1_0*std::sin(x2_0);
    z_0 = x3_0;
    b2 = calc_b_azimuthal_inner(pin, b1, v1_inner);
    b3 = 0.0; // z direction is zero
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x_0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y_0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z_0 = x1_0*std::cos(x2_0);
    b2 = 0.0; // polar is zero 
    b3 = calc_b_azimuthal_inner(pin, b1, v1_inner);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in wind.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // calculate normalized gravity potential
  if (sun_gravity) {
      GM_norm = calc_GM_sun(pin);
      // std::cout << GM_norm << " yo" << std::endl;
  }

  if (enable_user_dt) {
    user_dt = pin->GetOrAddReal("problem", "user_dt", 0.01);
  }

  // oscillation of current sheet
  wave_rad = pin->GetOrAddReal("problem", "wave_rad", 0.0);
  wave_mult = pin->GetOrAddReal("problem", "wave_mult", 1.0);
#endif
  //////////////////

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Solar wind test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

#if 0
  gamma_param = peos->GetGamma();
  // would need to change this if not cylindrical
  // get minimum x1 grid location get inner radius
  Real x1_min = pin->GetReal("mesh", "x1min");
  inner_radius = pin->GetOrAddReal("problem", "radius", x1_min);
  // r_measure = pin->GetOrAddReal("problem", "r_measure", 1.0);
  // m_p = pin->GetOrAddReal("problem", "m_p", 1.67E-27);
  // Real omega = pin->GetOrAddReal("problem", "omega_sun", 2.87e-6);
  //assume orbital_advection is enabled
  // Real omega_sun = pin->GetReal("orbital_advection", "Omega0");
  // TODO need to make it clear the coordinate systems and inputs
  v1_inner = calc_v_radial_inner(pin);
  v2_inner = 0.0; // create calc for v2 and v3 when orbital advection is off
  v3_inner = 0.0;
  n_inner = calc_n_inner(pin, v1_inner);
  e_inner = calc_energy_inner(pin, n_inner, gamma_param);

  // defaults
  b1 = calc_b_radial_inner(pin);
  b2 = calc_b_azimuthal_inner(pin, b1, v1_inner);
  b3 = 0.0;

  //output for each meshblock
  std::cout << "n_inner = " << n_inner << std::endl;
  std::cout << "v1_inner = " << v1_inner << std::endl;
  std::cout << "e_inner = " << e_inner << std::endl;
  std::cout << "b1 = " << b1 << std::endl;
  std::cout << "b3 = " << b2 << std::endl;
  
  // CME data
  Real t_o = pin->GetReal("problem", "t_o");
  Real vo = pin->GetReal("problem", "vo");

  //CME_start = pin->GetReal("problem", "CME_timestart_hrs")*3600/t_o;
  CME_start = pin->GetOrAddReal("problem", "CME_timestart_hrs",1e20)*3600/t_o;
  CME_duration = pin->GetOrAddReal("problem", "CME_timespan_hrs",0.)*3600/t_o;

  CME_density = pin->GetOrAddReal("problem", "CME_density",1.);
  CME_velocity = pin->GetOrAddReal("problem", "CME_velocity",0.)*1000/vo;
  CME_energy = e_inner*2;

  CME_phi_0 = (PI/180.0)*pin->GetOrAddReal("problem", "CME_phi_0",100.);
  CME_phi_extent = (PI/180.0)*pin->GetOrAddReal("problem", "CME_phi_extent",0.);

  // Real b0; //, angle;
  // if (MAGNETIC_FIELDS_ENABLED) {
  //   b0 = pin->GetReal("problem", "b0");
  //   angle = (PI/180.0)*pin->GetReal("problem", "angle");
  // }

  // get coordinates of center of sun, and convert to Cartesian if necessary
  // should all be zero for now
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x_0 = x1_0;
    y_0 = x2_0;
    z_0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x_0 = x1_0*std::cos(x2_0);
    y_0 = x1_0*std::sin(x2_0);
    z_0 = x3_0;
    b2 = calc_b_azimuthal_inner(pin, b1, v1_inner);
    b3 = 0.0; // z direction is zero
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x_0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y_0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z_0 = x1_0*std::cos(x2_0);
    b2 = 0.0; // polar is zero 
    b3 = calc_b_azimuthal_inner(pin, b1, v1_inner);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in wind.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // calculate normalized gravity potential
  if (sun_gravity) {
      GM_norm = calc_GM_sun(pin);
      // std::cout << GM_norm << " yo" << std::endl;
  }

  if (enable_user_dt) {
    user_dt = pin->GetOrAddReal("problem", "user_dt", 0.01);
  }

  // oscillation of current sheet
  wave_rad = pin->GetOrAddReal("problem", "wave_rad", 0.0);
  wave_mult = pin->GetOrAddReal("problem", "wave_mult", 1.0);
#endif

  // setup uniform ambient medium
  // Modifies density, and energy (non-barotropic eos and relativistic dynamics)
  Real rad, den, energy, x, y, z;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          x = pcoord->x1v(i);
          y = pcoord->x2v(j);
          z = pcoord->x3v(k);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          z = pcoord->x3v(k);
        } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          z = pcoord->x1v(i)*std::cos(pcoord->x2v(j)); 
        }
        rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));

        den = n_inner;
        energy = e_inner;
        if (rad >= inner_radius) {
          den *= SQR((inner_radius / rad));
          energy *= std::pow((inner_radius / rad), 2.0*gamma_param);
        }
        // ambient density and conserved variables
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i) * v1_inner;
        //phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i) * v2_inner;
        //phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i) * v3_inner;

        // currently assumes that inner velocities are given in same coordinate system
        // functions that generate values need to be functions if input
        // values need to be converted to correct coordinate system
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i) * v2_inner;
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i) * v3_inner;
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i) * v2_inner * rad;
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i) * v3_inner;
        } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          // v2_inner is polar and v3_inner is azimuthal
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i) * v2_inner * rad;
          phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i) * v3_inner * rad * std::sin(pcoord->x2v(j));
        }

        if (NON_BAROTROPIC_EOS) {
          // 0.5 * den *(v_r^2 + v_phi^2 + v_theta^2)
          // but move den inside and use momentum
          Real v_contrib;
          v_contrib = (0.5 / phydro->u(IDN,k,j,i)
                       * (SQR(phydro->u(IM1,k,j,i))
                          + SQR(phydro->u(IM2,k,j,i))
                          + SQR(phydro->u(IM3,k,j,i))
                          )
          );
          phydro->u(IEN,k,j,i) = energy + v_contrib;

          if (RELATIVISTIC_DYNAMICS) {
            std::stringstream msg;
            msg << "### FATAL ERROR in wind.cpp ProblemGenerator" << std::endl
            << "can't handle RELATIVISTIC_DYNAMICS=" << RELATIVISTIC_DYNAMICS << std::endl;
            ATHENA_ERROR(msg);
          }
          // if (RELATIVISTIC_DYNAMICS) { // this should only ever be SR with this file
          //   phydro->u(IEN,k,j,i) += den;
          // }
        }
      }
    }
  }

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    Real polar_dependence = 1.0;
    Real sign = 1.0;
    // b.x1f
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x = pcoord->x1f(i);
            y = pcoord->x2v(j);
            z = pcoord->x3v(k);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            x = pcoord->x1f(i)*std::cos(pcoord->x2v(j));
            y = pcoord->x1f(i)*std::sin(pcoord->x2v(j));
            z = pcoord->x3v(k);
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            x = pcoord->x1f(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
            y = pcoord->x1f(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
            z = pcoord->x1f(i)*std::cos(pcoord->x2v(j));
            // sign = std::copysign(1.0, std::cos(pcoord->x2v(j) + wave_rad*std::cos(wave_mult*pcoord->x3v(k))));
            sign = sign_radial_mag_field(pcoord->x2v(j), pcoord->x3v(k));
          }

          rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
          pfield->b.x1f(k,j,i) = b1*sign;
          if (rad >= inner_radius) {
            pfield->b.x1f(k,j,i) *= pow(inner_radius/rad, 2.0);
          }
        }
      }
    }

    //b.x2f
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x = pcoord->x1v(i);
            y = pcoord->x2f(j);
            z = pcoord->x3v(k);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            x = pcoord->x1v(i)*std::cos(pcoord->x2f(j));
            y = pcoord->x1v(i)*std::sin(pcoord->x2f(j));
            z = pcoord->x3v(k);
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            x = pcoord->x1v(i)*std::sin(pcoord->x2f(j))*std::cos(pcoord->x3v(k));
            y = pcoord->x1v(i)*std::sin(pcoord->x2f(j))*std::sin(pcoord->x3v(k));
            z = pcoord->x1v(i)*std::cos(pcoord->x2f(j));
          }
          rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
          pfield->b.x2f(k,j,i) = b2;
          if (rad >= inner_radius) {
            pfield->b.x2f(k,j,i) *= (inner_radius/rad);
          }
        }
      }
    }
    // b.x3f
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x = pcoord->x1v(i);
            y = pcoord->x2v(j);
            z = pcoord->x3f(k);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
            y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
            z = pcoord->x3f(k);
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3f(k));
            y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3f(k));
            z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
            polar_dependence = std::sin(pcoord->x2v(j));
          }
          rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
          pfield->b.x3f(k,j,i) = b3*polar_dependence; //*sign;
          if (rad >= inner_radius) {
            pfield->b.x3f(k,j,i) *= (inner_radius/rad);
          }
        }
      }
    }

    // add magnetic field contribution to total energy
    // using face averaged values
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) += 0.5 * (SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
                                           SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
                                           SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))
                                           );
          }
        }
      }
    }
  }
}

// add gravity potential term
// currently defined spherically for the sun at the center of the simulation
// converts to correct coordinates
void SunGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar) {

  Real rad, rho, rho_0, x, y, z;
  Real grav_force, f1, f2, f3;

  // rho is the radius direction in cylindrical coordinates
  rho_0 = std::sqrt(SQR(x_0) + SQR(y_0));
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          x = pmb->pcoord->x1v(i);
          y = pmb->pcoord->x2v(j);
          z = pmb->pcoord->x3v(k);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          x = pmb->pcoord->x1v(i)*std::cos(pmb->pcoord->x2v(j));
          y = pmb->pcoord->x1v(i)*std::sin(pmb->pcoord->x2v(j));
          z = pmb->pcoord->x3v(k);
          rho = std::sqrt(SQR(x) + SQR(y));
        } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          x = pmb->pcoord->x1v(i)*std::sin(pmb->pcoord->x2v(j))*std::cos(pmb->pcoord->x3v(k));
          y = pmb->pcoord->x1v(i)*std::sin(pmb->pcoord->x2v(j))*std::sin(pmb->pcoord->x3v(k));
          z = pmb->pcoord->x1v(i)*std::cos(pmb->pcoord->x2v(j)); 
        }

        rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
        grav_force = - GM_norm / SQR(rad) * cons(IDN,k,j,i);

        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          f1 = dt*grav_force * (x-x_0) / rad;
          f2 = dt*grav_force * (y-y_0) / rad;
          f3 = dt*grav_force * (z-z_0) / rad;
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // not 100% sure this is correct
          f1 = dt*grav_force * (rho-rho_0) / rad;
          f2 = 0.0;
          f3 = dt*grav_force * (z-z_0) / rad;
        } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          f1 = dt*grav_force;
          f2 = 0.0;
          f3 = 0.0;
        }

        cons(IM1,k,j,i) += f1;
        cons(IM2,k,j,i) += f2;
        cons(IM3,k,j,i) += f3;

        // multiply gravitational potential by smoothing function
        // cons(IM3,k,j,i) -= dt*den*SQR(Omega_0)*x3*fsmooth;
        if (NON_BAROTROPIC_EOS) {
          cons(IEN,k,j,i) += dt*prim(IDN,k,j,i) * (prim(IVX,k,j,i) * f1
                                                   + prim(IVY,k,j,i) * f2
                                                   + prim(IVZ,k,j,i) * f3);
        }
      }
    }
  }
  return;
}

Real SetMinimumTimeStep(MeshBlock *pmb)
{
  return user_dt;
}

//! User-defined boundary Conditions: sets solution in ghost zones to initial values
void CMEInnerX1(MeshBlock *pmb, Coordinates *pcoord,
                 AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // std::cout << "il: "  << il << " iu: " << iu;
  // std::cout << " jl: " << jl << " ju: " << ju;
  // std::cout << " kl: " << kl << " ku: " << ku;
  // std::cout << " ng: " << ngh;
  // std::cout << std::endl;

  Real den, vel, press, rad, x, y, z;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int gi=1; gi<=ngh; ++gi) {
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          x = pcoord->x1v(il-gi);
          y = pcoord->x2v(j);
          z = pcoord->x3v(k);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          x = pcoord->x1v(il-gi)*std::cos(pcoord->x2v(j));
          y = pcoord->x1v(il-gi)*std::sin(pcoord->x2v(j));
          z = pcoord->x3v(k);
        } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          x = pcoord->x1v(il-gi)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          y = pcoord->x1v(il-gi)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          z = pcoord->x1v(il-gi)*std::cos(pcoord->x2v(j)); 
        }
        rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));

        // may need to be looked at again
        den = n_inner;
	//std::cout << "bc: n_inner = " << n_inner << std::endl;
        press = e_inner * (gamma_param-1.0);
        if (rad >= inner_radius) {
          den *= SQR((inner_radius / rad));
          press *= std::pow((inner_radius / rad), 2.0*gamma_param);
        }
	
	vel = v1_inner;
	if (time > CME_start && time < CME_start + CME_duration) {
	  if (CME_phi_extent > 0.) {
	    Real cme_center = CME_phi_0;
	    Real cme_width = CME_phi_extent;
	    Real cme_var = (cme_width/2.3556)*(cme_width/2.3556); //variance of normal distr

	    Real Gauss_prof = exp(-(pcoord->x3v(k) - cme_center)*(pcoord->x3v(k) - cme_center)/(2.*cme_var));
	    Real CME_press = e_inner*(gamma_param-1.0)*(CME_density/n_inner) * 1.33; //*2.
	    Real press_inner = e_inner * (gamma_param-1.0);
	      	    
	    den = n_inner + (CME_density - n_inner)*Gauss_prof;
	    vel = v1_inner + (CME_velocity - v1_inner)*Gauss_prof;
	    //press = e_inner*(gamma_param-1.0)*(CME_density/n_inner) * 2.;
	    press = press_inner + (CME_press - press_inner)*Gauss_prof;
	  }
	  else {
	    den = CME_density;
	    vel = CME_velocity;
	    press = e_inner*(gamma_param-1.0)*(CME_density/n_inner) * 2.;
	  }
	}
	
	//fprintf(stderr, "den = %e \n", den);

        prim(IDN,k,j,il-gi) = den;
        ////prim(IVX,k,j,il-gi) = v1_inner;
	prim(IVX,k,j,il-gi) = vel;



        // prim(IVY,k,j,il-gi) = v2_inner;
        // prim(IVZ,k,j,il-gi) = v3_inner;
        // currently assumes that inner velocities are given in same coordinate system
        // functions that generate values need to be functions if input
        // values need to be converted to correct coordinate system
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          prim(IVY,k,j,il-gi) = v2_inner;
          prim(IVZ,k,j,il-gi) = v3_inner;
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          prim(IVY,k,j,il-gi) = v2_inner * rad;
          prim(IVZ,k,j,il-gi) = v3_inner;
        } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          prim(IVY,k,j,il-gi) = v2_inner * rad;
          prim(IVZ,k,j,il-gi) = v3_inner * rad * std::sin(pcoord->x2v(j));
        }

        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,il-gi) = press;
        }
      }
    }
  }

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real polar_dependence = 1.0;
    Real sign = 1.0;
    // std::cout << "x1f" << std::endl;
    // x1 direction , use volume centered for x2 and x3
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int gi=1; gi<=ngh; ++gi) {
          // mirror BC 
          //b.x1f(k,j,il-gi) = b.x1f(k,j,il);
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x = pcoord->x1f(il-gi);
            y = pcoord->x2v(j);
            z = pcoord->x3v(k);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            x = pcoord->x1f(il-gi)*std::cos(pcoord->x2v(j));
            y = pcoord->x1f(il-gi)*std::sin(pcoord->x2v(j));
            z = pcoord->x3v(k);
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            x = pcoord->x1f(il-gi)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
            y = pcoord->x1f(il-gi)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
            z = pcoord->x1f(il-gi)*std::cos(pcoord->x2v(j));
            // sign = std::copysign(1.0, std::cos(pcoord->x2v(j) + wave_rad*std::cos(wave_mult*pcoord->x3v(k))));
            sign = sign_radial_mag_field(pcoord->x2v(j), pcoord->x3v(k));
          }

          rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
          b.x1f(k,j,il-gi) = b1*sign;
          if (rad >= inner_radius) {
            b.x1f(k,j,il-gi) *= SQR((inner_radius/rad));
          }
        }
      }
    }

    // std::cout << "x2f" << std::endl;
    // x2 direction , use volume centered for x1 and x3
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int gi=1; gi<=ngh; ++gi) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x = pcoord->x1v(il-gi);
            y = pcoord->x2f(j);
            z = pcoord->x3v(k);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            x = pcoord->x1v(il-gi)*std::cos(pcoord->x2f(j));
            y = pcoord->x1v(il-gi)*std::sin(pcoord->x2f(j));
            z = pcoord->x3v(k);
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            x = pcoord->x1v(il-gi)*std::sin(pcoord->x2f(j))*std::cos(pcoord->x3v(k));
            y = pcoord->x1v(il-gi)*std::sin(pcoord->x2f(j))*std::sin(pcoord->x3v(k));
            z = pcoord->x1v(il-gi)*std::cos(pcoord->x2f(j)); 
          }
          // std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
          rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
          b.x2f(k,j,il-gi) = b2;
          if (rad >= inner_radius) {
            b.x2f(k,j,il-gi) *= (inner_radius/rad);
          }
        }
      }
    }

    // std::cout << "x3f" << std::endl;
    // x3 direction , use volume centered for x1 and x2
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int gi=1; gi<=ngh; ++gi) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            x = pcoord->x1v(il-gi);
            y = pcoord->x2v(j);
            z = pcoord->x3f(k);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            x = pcoord->x1v(il-gi)*std::cos(pcoord->x2v(j));
            y = pcoord->x1v(il-gi)*std::sin(pcoord->x2v(j));
            z = pcoord->x3f(k);
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            x = pcoord->x1v(il-gi)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3f(k));
            y = pcoord->x1v(il-gi)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3f(k));
            z = pcoord->x1v(il-gi)*std::cos(pcoord->x2v(j));
            polar_dependence = std::sin(pcoord->x2v(j));
            // sign = -1.0*std::copysign(1.0, std::cos(pcoord->x2v(j) + wave_rad*std::cos(wave_mult*pcoord->x3v(k))));
          }
          // std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
          rad = std::sqrt(SQR(x - x_0) + SQR(y - y_0) + SQR(z - z_0));
          b.x3f(k,j,il-gi) = b3*polar_dependence; //*sign;
          if (rad >= inner_radius) {
            b.x3f(k,j,il-gi) *= (inner_radius/rad);
          }
        }
      }
    }

    // std::stringstream msg;
    // msg << "### poopy" << std::endl;
    // ATHENA_ERROR(msg);
    // add magnetic field contribution to total pressure
    // using face averaged values
    if (NON_BAROTROPIC_EOS) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int gi=1; gi<=ngh; ++gi) {
            prim(IPR,k,j,il-gi) += 0.5 * (SQR(0.5*(b.x1f(k,j,il-gi) + b.x1f(k,j,il-gi+1))) +
                                         SQR(0.5*(b.x2f(k,j,il-gi) + b.x2f(k,j+1,il-gi))) +
                                         SQR(0.5*(b.x3f(k,j,il-gi) + b.x3f(k+1,j,il-gi)))
                                        );
          }
        }
      }
    }
  }
}

//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void CMEOuterX1(MeshBlock *pmb,Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy variables into ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int gi=1; gi<=ngh; ++gi) {
        prim(IDN,k,j,iu+gi) = prim(IDN,k,j,iu);
        prim(IVX,k,j,iu+gi) = prim(IVX,k,j,iu);
        prim(IVY,k,j,iu+gi) = prim(IVY,k,j,iu);
        prim(IVZ,k,j,iu+gi) = prim(IVZ,k,j,iu);

        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,iu+gi) = prim(IPR,k,j,iu);
        }
      }
    }
  }

  // set magnetic field in outlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int gi=1; gi<=ngh; ++gi) {
            b.x1f(k,j,iu+gi) = b.x1f(k,j,iu);
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int gi=1; gi<=ngh; ++gi) {
          b.x2f(k,j,iu+gi) = b.x2f(k,j,iu);
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int gi=1; gi<=ngh; ++gi) {
          b.x3f(k,j,iu+gi) = b.x3f(k,j,iu);
        }
      }
    }

//     if (NON_BAROTROPIC_EOS) {
//       for (int k=kl; k<=ku; ++k) {
//         for (int j=jl; j<=ju; ++j) {
// #pragma omp simd
//           for (int gi=1; gi<=ngh; ++gi) {
//             prim(IPR,k,j,iu+gi) -= 0.5 * (SQR(0.5*(b.x1f(k,j,iu+gi) + b.x1f(k,j,iu+gi+1))) +
//                                          SQR(0.5*(b.x2f(k,j,iu+gi) + b.x2f(k,j+1,iu+gi))) +
//                                          SQR(0.5*(b.x3f(k,j,iu+gi) + b.x3f(k+1,j,iu+gi)))
//                                         );
//           }
//         }
//       }
//     }
  }
}

// refinement condition: check the maximum pressure gradient
int RefinementConditionPressure(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  if (pmb->pmy_mesh->f3) {
    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                               +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                               +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
          maxeps = std::max(maxeps, eps);
        }
      }
    }
  } else if (pmb->pmy_mesh->f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                             + SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
        maxeps = std::max(maxeps, eps);
      }
    }
  } else {
    return 0;
  }

  if (maxeps > press_threshold) return 1;
  if (maxeps < 0.25*press_threshold) return -1;
  return 0;
}

// refinement condition: density jump
int RefinementConditionDensityJump(MeshBlock *pmb) {
  int f2 = pmb->pmy_mesh->f2, f3 = pmb->pmy_mesh->f3;
  AthenaArray<Real> &w = pmb->phydro->w;
  // maximum intercell density ratio
  Real drmax = 1.0;
  for (int k=pmb->ks-f3; k<=pmb->ke+f3; k++) {
    for (int j=pmb->js-f2; j<=pmb->je+f2; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        if (w(IDN,k,j,i-1)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j,i-1)/w(IDN,k,j,i);
        if (w(IDN,k,j,i+1)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j,i+1)/w(IDN,k,j,i);
        if (w(IDN,k,j,i)/w(IDN,k,j,i-1) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j,i-1);
        if (w(IDN,k,j,i)/w(IDN,k,j,i+1) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j,i+1);
        if (f2) {
          if (w(IDN,k,j-1,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j-1,i)/w(IDN,k,j,i);
          if (w(IDN,k,j+1,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j+1,i)/w(IDN,k,j,i);
          if (w(IDN,k,j,i)/w(IDN,k,j-1,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j-1,i);
          if (w(IDN,k,j,i)/w(IDN,k,j+1,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j+1,i);
        }
        if (f3) {
          if (w(IDN,k-1,j,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k-1,j,i)/w(IDN,k,j,i);
          if (w(IDN,k+1,j,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k+1,j,i)/w(IDN,k,j,i);
          if (w(IDN,k,j,i)/w(IDN,k-1,j,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k-1,j,i);
          if (w(IDN,k,j,i)/w(IDN,k+1,j,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k+1,j,i);
        }
      }
    }
  }
  if (drmax > 1.5) return 1;
  else if (drmax < 1.2) return -1;
  return 0;
}

// get normalized gravitational constant
Real calc_GM_sun(ParameterInput *pin) {
  Real AU = pin->GetReal("problem", "AU");
  Real t_o = pin->GetReal("problem", "t_o");
  Real GM_sun = pin->GetOrAddReal("problem", "GM_sun", 1.32712440018e20);
  Real GM_0 = GM_sun * SQR(t_o) / std::pow(AU,3);
  return GM_0;
}

// convert velocity to inner boundary velocity
// empirical formulation
Real calc_v_radial_inner(ParameterInput *pin) {
  Real v_measure = pin->GetReal("problem", "v_measure");
  Real v_measure_extra = pin->GetOrAddReal("problem", "v_measure_extra", 1.0);
  Real vo = pin->GetReal("problem", "vo");
  Real v_in = (std::sqrt(v_measure / 430.7) * 0.8231 *0.6
                  * v_measure_extra * v_measure * 1000.0
                 )
                  / vo;

  //Real t_o = pin->GetReal("problem", "t_o");
  //Real CME_start = pin->GetReal("problem", "CME_timestart_hrs")*3600/t_o;
  //Real CME_duration = pin->GetReal("problem", "CME_timespan_hrs")*3600/t_o;

  //if (time > CME_start && time < CME_start + CME_duration)
  //  v_in = pin->GetReal("problem", "CME_velocity")*1000/vo;

  return v_in;
}



// convert b_radial magnetic field to inner boundary 
Real calc_b_radial_inner(ParameterInput *pin) {
  Real min_radius = pin->GetReal("mesh", "x1min");
  Real r_meas = pin->GetReal("problem", "r_measure");
  Real v_measure = pin->GetReal("problem", "v_measure");
  Real b_measure = pin->GetReal("problem", "b_measure");
  // Real omega = pin->GetReal("orbital_advection", "Omega0");
  Real omega = pin->GetReal("problem", "omega_sun");
  Real AU = pin->GetReal("problem", "AU");
  Real bo = pin->GetReal("problem", "bo");

  Real B_AU = b_measure * 1.0E-9 / bo;
  Real ratio = (omega * AU) / (v_measure * 1000.0);
  Real B1_AU = B_AU / (std::sqrt(1.0 + SQR(ratio)));
  Real b1 = B1_AU * SQR((r_meas / min_radius));

  return b1;
}

// convert b2 magnetic field to inner boundary 
Real calc_b_azimuthal_inner(ParameterInput *pin, Real _b1, Real _v1_inner) {
  Real min_radius = pin->GetReal("mesh", "x1min");
  // Real omega = pin->GetReal("orbital_advection", "Omega0");
  Real omega = pin->GetReal("problem", "omega_sun");
  Real AU = pin->GetReal("problem", "AU");
  Real vo = pin->GetReal("problem", "vo");
  Real b_extra = pin->GetOrAddReal("problem", "b_measure_extra", 1.0);

  Real ratio_inner = omega * AU * min_radius / v1_inner / vo;
  Real _b_az = -1.0* _b1 * ratio_inner * b_extra * 0.73 ;
  return _b_az;
}

// convert density to inner boundary 
Real calc_n_inner(ParameterInput *pin, Real v1_inner) {
  Real min_radius = pin->GetReal("mesh", "x1min");
  Real r_meas = pin->GetReal("problem", "r_measure");
  Real v_measure = pin->GetReal("problem", "v_measure");
  Real n_measure = pin->GetReal("problem", "n_measure");
  Real vo = pin->GetReal("problem", "vo");

  Real n_in = (n_measure * SQR((r_meas / min_radius)) * v_measure * 1000.0) / (v1_inner * vo)*0.73;
  return n_in;
}

// convert density and temperature to energy density 
Real calc_energy_inner(ParameterInput *pin, Real n_inner, Real gamma) {
  Real min_radius = pin->GetReal("mesh", "x1min");
  Real r_meas = pin->GetReal("problem", "r_measure");
  Real T_measure = pin->GetReal("problem", "T_measure");
  // Real gamma = peos->GetGamma();

  // convert temperature to inner boundary temperature
  Real k_inner = T_measure * pow((r_meas/min_radius), (4.0/3.0));
  // convert density and temperature to energy density 
  Real es = (3.0 * n_inner * k_inner) / (gamma - 1.0);
  return es;
}

// return sign of magnetic field
Real sign_radial_mag_field(Real theta, Real phi) {
  Real sign, theta_perturb;
  // theta_perturb = wave_rad*std::cos(phi - omega_sun*(time - radius/v1_inner));
  theta_perturb = wave_rad*std::cos(wave_mult*phi);
  sign = 1.; ////std::copysign(1.0, std::cos(theta + theta_perturb));
  return sign;
}
