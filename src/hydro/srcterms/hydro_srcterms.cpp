//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_srcterms.cpp
//! \brief Class to implement source terms in the hydro equations

// C headers

// C++ headers
#include <cstring>    // strcmp
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../../parameter_input.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//! HydroSourceTerms constructor

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin) {
  pmy_hydro_ = phyd;
  hydro_sourceterms_defined = false;

  // read point mass or constant acceleration parameters from input block

  // set the point source only when the coordinate is spherical or 2D
  // It works even for cylindrical with the orbital advection.
  flag_point_mass_ = false;
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  bool orbital_advection_defined
         = (pin->GetOrAddInteger("orbital_advection","OAorder",0)!=0)?
           true : false;
  if (gm_ != 0.0) {
    if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0
        && std::strcmp(COORDINATE_SYSTEM, "cylindrical") != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point mass gravity works only in the cylindrical and "
          << "spherical polar coordinates." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (orbital_advection_defined) {
      hydro_sourceterms_defined = true;
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
               && phyd->pmy_block->block_size.nx3>1) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point mass gravity deos not work in the 3D cylindrical "
          << "coordinates without orbital advection." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    } else {
      flag_point_mass_ = true;
      hydro_sourceterms_defined = true;
    }
  }
  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  if (g1_ != 0.0) hydro_sourceterms_defined = true;

  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  if (g2_ != 0.0) hydro_sourceterms_defined = true;

  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  if (g3_ != 0.0) hydro_sourceterms_defined = true;

  // read shearing box parameters from input block
  Omega_0_ = pin->GetOrAddReal("orbital_advection","Omega0",0.0);
  qshear_  = pin->GetOrAddReal("orbital_advection","qshear",0.0);
  ShBoxCoord_ = pin->GetOrAddInteger("orbital_advection","shboxcoord",1);

  // check flag for shearing source
  flag_shearing_source_ = 0;
  if(orbital_advection_defined) { // orbital advection source terms
    if(ShBoxCoord_ == 1) {
      flag_shearing_source_ = 1;
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "OrbitalAdvection does NOT work with shboxcoord = 2." << std::endl
          << "Check <orbital_advection> shboxcoord parameter in the input file."
          << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if ((Omega_0_ !=0.0) && (qshear_ != 0.0)
             && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    flag_shearing_source_ = 2; // shearing box source terms
  } else if ((Omega_0_ != 0.0) &&
             (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
              || std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)) {
    flag_shearing_source_ = 3; // rotating system source terms
  }

  if (flag_shearing_source_ != 0)
    hydro_sourceterms_defined = true;

  if (SELF_GRAVITY_ENABLED) hydro_sourceterms_defined = true;

  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
  if (UserSourceTerm != nullptr) hydro_sourceterms_defined = true;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::AddHydroSourceTerms
//! \brief Adds source terms to conserved variables
//! This function is not only for hydro variables but also for passive scalars.

void HydroSourceTerms::AddSourceTerms(const Real time, const Real dt,
                                      const AthenaArray<Real> *flux,
                                      const AthenaArray<Real> &prim,
                                      const AthenaArray<Real> &prim_scalar,
                                      const AthenaArray<Real> &bcc,
                                      AthenaArray<Real> &cons,
                                      AthenaArray<Real> &cons_scalar) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // accleration due to point mass (MUST BE AT ORIGIN)
  if (flag_point_mass_)
    PointMass(dt, flux, prim, cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0)
    ConstantAcceleration(dt, flux, prim, cons);

  // Add new source terms here
  if (SELF_GRAVITY_ENABLED)
    SelfGravity(dt, flux, prim, cons);

  // Sorce terms for orbital advection, shearing box, or rotating system
  if (flag_shearing_source_ == 1)
    OrbitalAdvectionSourceTerms(dt, flux, prim, cons);
  else if (flag_shearing_source_ == 2)
    ShearingBoxSourceTerms(dt, flux, prim, cons);
  else if (flag_shearing_source_ == 3)
    RotatingSystemSourceTerms(dt, flux, prim, cons);

  // MyNewSourceTerms()

  //  user-defined source terms
  if (UserSourceTerm != nullptr) {
    UserSourceTerm(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }

  return;
}
