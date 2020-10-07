//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file planet2D.cpp
//  \brief problem generator of a 2D disk with a planet
//  This works only in Cartesian or Cylindrical coordinates.
//
//======================================================================================

// C++ headers
#include <cmath>      // sqrt()
#include <iomanip>   // setprecision
#include <iostream>   // cout, endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires NOT MHD"
#endif

namespace {
Real CsProfileCyl(const Real rad, const Real phi, const Real z);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real OrbitalVelocityCyl(Real &den, Real &cs,
                        const Real rad, const Real phi, const Real z);

//parameter
static Real gamma_gas;
static Real d0, dslope;
static Real c0, cslope;
static Real dfloor, pfloor;
static Real Omega0, qshear; // shearing box parameter for Cartesian
static Real gm;             // central star gravity for Cylindrical
static Real ap, Mp, epsilon2;
static bool LocallyIsothermal, IndirectTerm;
} // namespace

// B.C.
// fixed boundary
//sets BCs on inner-x1 (left edge) of grid.
void fix_iib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu,
             int jl, int ju, int kl, int ku, int ngh);
//sets BCs on outer-x1 (right edge) of grid.
void fix_oib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu,
             int jl, int ju, int kl, int ku, int ngh);

// Sorce function
void SourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                AthenaArray<Real> &cons_scalar);

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (mesh_size.nx2 == 1 || mesh_size.nx3 > 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in InitUserMeshData" << std::endl
        << "This problem generator works only in 2D."<< std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // initialize global variables
  Real float_min = std::numeric_limits<float>::min();
  d0     = pin->GetOrAddReal("problem","d0",1.0);
  dslope = pin->GetOrAddReal("problem","dslope",0.0);
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));
  pfloor=pin->GetOrAddReal("hydro","pfloor",(1024*(float_min)));

  if (NON_BAROTROPIC_EOS) {
    c0     = pin->GetOrAddReal("problem","c0", 0.05);
    cslope = pin->GetOrAddReal("problem","cslope", 0.0);
    LocallyIsothermal = pin->GetOrAddBoolean("problem","LocallyIsothermal", false);
    gamma_gas      = pin->GetReal("hydro","gamma");
  } else { // isothermal
    c0     = pin->GetReal("hydro","iso_sound_speed");
    cslope = 0.0;
  }

  Mp     = pin->GetOrAddReal("problem","Mp",0.0);
  epsilon2  = SQR(pin->GetOrAddReal("problem","epsilon",0.03));
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    ap     = pin->GetOrAddReal("problem","ap",1.0);
    IndirectTerm = pin->GetOrAddBoolean("problem","IndirectTerm", false);
    if (ap <= 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in InitUserMeshData" << std::endl
          << "Semi-major axis of a planet must be positive."<< std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }

  //Source term
  EnrollUserExplicitSourceFunction(SourceTerm);

  //Boundary Conditions
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, fix_iib);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, fix_oib);
  }

  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief problem generator for a 2D disk with a planet.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Omega0 = porb->Omega0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    qshear  = porb->qshear;
    if (dslope != 0.0 || cslope != 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ProblemGenerator" << std::endl
          << "Initial profiles should be flat in catesian coordinates."<< std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= d0*qshear*Omega0*pcoord->x1v(i);
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = SQR(c0)*d0/(gamma_gas-1.0)
                                   + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                         +SQR(phydro->u(IM2,k,j,i))
                                         +SQR(phydro->u(IM3,k,j,i)))/d0;
          }
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    gm = porb->gm;
    Real den, cs;
    Real x1, x2, x3;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          x2 = pcoord->x2v(j);
          x3 = pcoord->x3v(k);
          Real ovel = OrbitalVelocityCyl(den, cs, x1, x2, x3);
          if(porb->orbital_advection_defined)
            ovel -= porb->OrbitalVelocity(porb, x1, x2, x3);
          phydro->u(IDN,k,j,i) = den;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = den*ovel;
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = SQR(cs)*den/(gamma_gas-1.0)
                                   +0.5*(SQR(phydro->u(IM1,k,j,i))
                                        +SQR(phydro->u(IM2,k,j,i))
                                        +SQR(phydro->u(IM3,k,j,i)))/den;
          }
        }
      }
    }
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//!\f give sound speed profile for cylindrical or spherical_polar coordinates
Real CsProfileCyl(const Real rad, const Real phi, const Real z) {
  return c0*std::pow(rad, cslope);
}

//----------------------------------------------------------------------------------------
//!\f give density profile for cylindrical coordinates
Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den = d0*std::pow(rad,dslope);
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//!\f give orbital velocity profile for cylindrical coordinates
Real OrbitalVelocityCyl(Real &den, Real &cs, const Real rad,
                        const Real phi, const Real z) {
  cs  = CsProfileCyl(rad, phi, z);
  den = DenProfileCyl(rad, phi, z);
  Real temp = (gm+(dslope+2.0*cslope)*SQR(cs))/rad;
  return std::sqrt(temp)-rad*Omega0;
}
} // namespace

//======================================================================================
//  Source term
//======================================================================================
void SourceTerm(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                AthenaArray<Real> &cons_scalar) {
  Real pgv1, pgv2, pgv3; // planet gravity
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    Real x1, x2, x3;
    // put central star at 20H away
    if (Mp>0.0) {
      Real GMp = 8000.0*(c0*c0*c0)/Omega0*Mp;
      for(int k=pmb->ks; k<=pmb->ke; ++k) {
        x3 = pmb->pcoord->x3v(k);
        for(int j=pmb->js; j<=pmb->je; ++j) {
          x2 = pmb->pcoord->x2v(j);
          for(int i=pmb->is; i<=pmb->ie; ++i) {
            x1 = pmb->pcoord->x1v(i);
            Real d2 = SQR(x1)+SQR(x2)+SQR(x3);
            Real temp = dt*prim(IDN,k,j,i)*GMp
                        *(d2+2.5*epsilon2)/(SQR(d2+epsilon2)*std::sqrt(d2+epsilon2));
            pgv1 = -temp*x1;
            pgv2 = -temp*x2;
            pgv3 = -temp*x3;
            cons(IM1,k,j,i) += -temp*x1;
            cons(IM2,k,j,i) += -temp*x2;
            cons(IM3,k,j,i) += -temp*x3;
            if (NON_BAROTROPIC_EOS && !LocallyIsothermal) {
              cons(IEN,k,j,i) += -temp
                                  *(x1*prim(IVX,k,j,i)
                                   +x2*prim(IVY,k,j,i)
                                   +x3*prim(IVZ,k,j,i));
            }
          }
        }
      }
    }
    if (NON_BAROTROPIC_EOS && LocallyIsothermal) {
      for(int k=pmb->ks; k<=pmb->ke; ++k) {
        for(int j=pmb->js; j<=pmb->je; ++j) {
          for(int i=pmb->is; i<=pmb->ie; ++i) {
            Real den = cons(IDN,k,j,i);
            cons(IEN,k,j,i) = SQR(c0)*den/(gamma_gas-1.0)
                                +0.5*(SQR(cons(IM1,k,j,i))
                                     +SQR(cons(IM2,k,j,i))
                                     +SQR(cons(IM3,k,j,i)))/den;
          }
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    Real rad, phi, z;
    if (Mp>0.0) {
      Real planet_omega = std::sqrt(gm/ap)/ap-Omega0;
      Real GMp = gm*Mp;
      Real ppos = std::fmod(planet_omega*time, 2.0*PI);
      for(int k=pmb->ks; k<=pmb->ke; ++k) {
        z = pmb->pcoord->x3v(k);
        for(int j=pmb->js; j<=pmb->je; ++j) {
          phi = pmb->pcoord->x2v(j);
          for(int i=pmb->is; i<=pmb->ie; ++i) {
            rad = pmb->pcoord->x1v(i);
            Real d2 = SQR(rad)+SQR(ap)-2.0*ap*rad*std::cos(phi-ppos)+SQR(z);
            Real temp = dt*prim(IDN,k,j,i)*GMp
                        *(d2+2.5*epsilon2)/(SQR(d2+epsilon2)*std::sqrt(d2+epsilon2));
            pgv1 = -temp*(rad-ap*std::cos(phi-ppos));
            pgv2 = -temp*ap*std::sin(phi-ppos);
            pgv3 = -temp*z;
            if (IndirectTerm) {
              temp = dt*prim(IDN,k,j,i)*GMp/SQR(ap);
              pgv1 += -temp*std::cos(phi-ppos);
              pgv2 +=  temp*std::sin(phi-ppos);
            }
            cons(IM1,k,j,i) += pgv1;
            cons(IM2,k,j,i) += pgv2;
            cons(IM3,k,j,i) += pgv3;
            if (NON_BAROTROPIC_EOS && !LocallyIsothermal) {
              cons(IEN,k,j,i) += pgv1*prim(IVX,k,j,i)+
                                 pgv2*prim(IVY,k,j,i)+
                                 pgv3*prim(IVZ,k,j,i);
            }
          }
        }
      }
    }
    if (NON_BAROTROPIC_EOS && LocallyIsothermal) {
      for(int k=pmb->ks; k<=pmb->ke; ++k) {
        for(int j=pmb->js; j<=pmb->je; ++j) {
          for(int i=pmb->is; i<=pmb->ie; ++i) {
            rad = pmb->pcoord->x1v(i);
            phi = pmb->pcoord->x2v(j);
            z   = pmb->pcoord->x3v(k);
            Real den = cons(IDN,k,j,i);
            Real cs = CsProfileCyl(rad,phi,z);
            cons(IEN,k,j,i) = SQR(cs)*den/(gamma_gas-1.0)
                               +0.5*(SQR(cons(IM1,k,j,i))
                                    +SQR(cons(IM2,k,j,i))
                                    +SQR(cons(IM3,k,j,i)))/den;
          }
        }
      }
    }
  }
  return;
}

//======================================================================================
//  Boundary Condtions, fixed, ix1, ox1
//======================================================================================
void fix_iib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int is, int ie, int js,
             int je, int ks, int ke, int ngh) {
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=1; i<=(NGHOST); i++) {
          prim(IDN,k,j,is-i) = d0;
          prim(IM1,k,j,is-i) = 0.0;
          prim(IM2,k,j,is-i) = 0.0;
          if(!pmb->porb->orbital_advection_defined)
            prim(IM2,k,j,is-i) -= qshear*Omega0*pco->x1v(is-i);
          prim(IM3,k,j,is-i) = 0.0;
          if(NON_BAROTROPIC_EOS)
            prim(IPR,k,j,is-i) = SQR(c0)*d0;
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=1; i<=(NGHOST); i++) {
          Real den, cs;
          Real ovel = OrbitalVelocityCyl(den, cs, pco->x1v(is-i),
                                         pco->x2v(j), pco->x3v(k));
          if(pmb->porb->orbital_advection_defined)
            ovel -= pmb->porb->vKc(k,is-i);
          prim(IDN,k,j,is-i) = den;
          prim(IM1,k,j,is-i) = 0.0;
          prim(IM2,k,j,is-i) = ovel;
          prim(IM3,k,j,is-i) = 0.0;
          if(NON_BAROTROPIC_EOS) prim(IPR,k,j,is-i) = SQR(cs)*den;
        }
      }
    }
  }
  return;
}

void fix_oib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int is, int ie,
             int js, int je, int ks, int ke, int ngh) {
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=1; i<=(NGHOST); i++) {
          prim(IDN,k,j,ie+i) = d0;
          prim(IM1,k,j,ie+i) = 0.0;
          prim(IM2,k,j,ie+i) = 0.0;
          if(!pmb->porb->orbital_advection_defined)
            prim(IM2,k,j,ie+i) -= qshear*Omega0*pco->x1v(ie+i);
          prim(IM3,k,j,ie+i) = 0.0;
          if(NON_BAROTROPIC_EOS)
            prim(IPR,k,j,ie+i) = SQR(c0)*d0;
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    Real den, cs;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=1; i<=(NGHOST); i++) {
          Real ovel = OrbitalVelocityCyl(den, cs, pco->x1v(ie+i),
                                         pco->x2v(j), pco->x3v(k));
          if(pmb->porb->orbital_advection_defined)
            ovel -= pmb->porb->vKc(k,ie+i);
          prim(IDN,k,j,ie+i) = den;
          prim(IM1,k,j,ie+i) = 0.0;
          prim(IM2,k,j,ie+i) = ovel;
          prim(IM3,k,j,ie+i) = 0.0;
          if(NON_BAROTROPIC_EOS) prim(IPR,k,j,ie+i) = SQR(cs)*den;
        }
      }
    }
  }
  return;
}
