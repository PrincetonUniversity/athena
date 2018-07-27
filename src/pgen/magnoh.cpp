//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file magnoh.cpp
//  \brief Magnetized Noh with perturbation in B_phi. Authored by A. Beresnyak
//
// 2D collapse on center, r distance to center
// initial conditions with r in cm:
// rho  = rho0*r^alpha [g/cm^3]
// V    = V0    [cm/s]
// Bphi = Bphi0*r^beta   [gauss] azimuthal
// Bz   = Bz0*  r^beta   [gauss] axial
// pressure = 1.E-6*B^2   actually zero in the exact solution
//
// Can apply sine wave perturbation in a form
//   *(1+perturb*cos(mphi*phi)) to the magnetic potential Az
//
// REFERENCES:
// 1) Velikovich, Giuliani, Zalesak, Gardiner, "Exact self-similar solutions for the
// magnetized Noh Z pinch problem", Phys. of Plasmas, vol.19, p.012707 (2012)
//
// 2) Giuliani, Velikovich, Beresnyak, Zalesak, Gianakon, Rousculp, "Self-similar
// solutions for the magnetized Noh problem with axial and azimuthal field", Phys. of
// Plasmas, in prep (2018)

// C++ headers
#include <algorithm>
#include <cmath>      // sqrt()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

static Real gm1;
static Real alpha, beta, rho0, P0, pcoeff, vr, perturb, mphi;
static Real bphi0, bz;
// static Real nu_iso, eta_ohm;

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  // nu_iso = pin->GetOrAddReal("problem", "nu_iso", 0.0);
  // eta_ohm = pin->GetOrAddReal("problem", "eta_ohm", 0.0);

  alpha  = pin->GetReal("problem", "alpha");
  beta  = pin->GetReal("problem", "beta");
  pcoeff= pin->GetReal("problem", "pcoeff");
  rho0  = pin->GetReal("problem", "d");
  vr =  pin->GetReal("problem", "vr");
  // convert from CGS to Athena Heaviside units:
  bphi0 = pin->GetReal("problem", "bphi")/std::sqrt(4*M_PI);
  bz = pin->GetReal("problem", "bz")/std::sqrt(4*M_PI);
  P0 = 4*M_PI*pcoeff*(bphi0*bphi0+bz*bz);

  perturb = pin->GetOrAddReal("problem","perturb",0.0);
  mphi = pin->GetOrAddReal("problem","mphi",1.0);

  // Add AMR condition based on density gradient
  if (adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);

  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for zpinch problem
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;

  if (COORDINATE_SYSTEM != "cartesian" &&  COORDINATE_SYSTEM != "cylindrical") {
    std::stringstream msg;
    msg << "### FATAL ERROR in magnoh.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl
        << "Only Cartesian and cylindrical are supported for this problem" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // initialize vector potential for inflowing B
  // (only initializing 2D array for vec potential)
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  AthenaArray<Real> az;
  az.NewAthenaArray(nx2,nx1);

  for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      Real rad,phi;
      if (COORDINATE_SYSTEM == "cylindrical") {
        rad = pcoord->x1f(i);
        phi = pcoord->x2f(j);
      } else { // cartesian
        Real x1  = pcoord->x1f(i);
        Real x2  = pcoord->x2f(j);
        rad = std::sqrt(SQR(x1) + SQR(x2));
        phi = atan2(x2, x1);
      }
      // if (rad==0.0) rad=3.0/100000;
      az(j,i) = (bphi0/(beta+1))*pow(rad,beta+1)*(1+perturb*cos(mphi*phi));
    }
  }

  // initialize conserved variables
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        // Volume centered coordinates and quantities
         Real rad,x1,x2;
         if(COORDINATE_SYSTEM == "cylindrical") {
             rad = pcoord->x1v(i);
         } else { // cartesian
             x1=pcoord->x1v(i);
             x2=pcoord->x2v(j);
             rad = std::sqrt(SQR(x1) + SQR(x2));
         }
         Real rho = rho0*pow(rad, alpha);
         Real P   = P0  *pow(rad, 2*beta);

         phydro->u(IDN,k,j,i) = rho;

         if(COORDINATE_SYSTEM == "cylindrical") {
             phydro->u(IM1,k,j,i) = rho*vr;
             phydro->u(IM2,k,j,i) = 0.0;
         } else { // cartesian
             phydro->u(IM1,k,j,i) = rho*vr*x1/rad;
             phydro->u(IM2,k,j,i) = rho*vr*x2/rad;
         }

         phydro->u(IM3,k,j,i) = 0.0;
         phydro->u(IEN,k,j,i) = P/gm1 + 0.5*rho*SQR(vr);
      }
    }
  }

  // initialize face-averaged magnetic fields
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real geom_coeff=1.0;
          if (COORDINATE_SYSTEM == "cylindrical") geom_coeff=1.0/pcoord->x1f(i);
          pfield->b.x1f(k,j,i) = geom_coeff*(az(j+1,i) - az(j,i))/pcoord->dx2f(j);
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real geom_coeff=1.0;
          if (COORDINATE_SYSTEM == "cylindrical") geom_coeff=-1.0; // Left hand system?
          pfield->b.x2f(k,j,i) = geom_coeff*(az(j,i) - az(j,i+1))/pcoord->dx1f(i);
        }
      }
    }
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real rad;
          if (COORDINATE_SYSTEM == "cylindrical") rad = pcoord->x1v(i);
          else rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)));
          pfield->b.x3f(k,j,i) = bz*pow(rad,beta);
        }
      }
    }
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            phydro->u(IEN,k,j,i) +=
                // second-order accurate assumption about volume-averaged field
                0.5*0.25*(SQR(pfield->b.x1f(k,j,i)+pfield->b.x1f(k,j,i+1))
                          + SQR(pfield->b.x2f(k,j,i)+pfield->b.x2f(k,j+1,i))
                          + SQR(pfield->b.x3f(k,j,i)+pfield->b.x3f(k+1,j,i)));
          }
        }
      }
    }
  }
  az.DeleteAthenaArray();
  return;
}

// refinement condition: density jump
int RefinementCondition(MeshBlock *pmb) {
  int f2=0, f3=0;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real drmax=1.0;
  if (pmb->block_size.nx2 > 1) f2 = 1;
  if (pmb->block_size.nx3 > 1) f3 = 1;
  for (int k=pmb->ks-f3; k<=pmb->ke+f3; k++) {
    for (int j=pmb->js-f2; j<=pmb->je+f2; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        if (w(IDN,k,j,i-1)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j,i-1)/w(IDN,k,j,i);
        if (w(IDN,k,j,i+1)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j,i+1)/w(IDN,k,j,i);
        if (w(IDN,k,j,i)/w(IDN,k,j,i-1) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j,i-1);
        if (w(IDN,k,j,i)/w(IDN,k,j,i+1) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j,i+1);
        if (f2==1) {
          if (w(IDN,k,j-1,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j-1,i)/w(IDN,k,j,i);
          if (w(IDN,k,j+1,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j+1,i)/w(IDN,k,j,i);
          if (w(IDN,k,j,i)/w(IDN,k,j-1,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j-1,i);
          if (w(IDN,k,j,i)/w(IDN,k,j+1,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j+1,i);
        }
        if (f3==1) {
          if (w(IDN,k-1,j,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k-1,j,i)/w(IDN,k,j,i);
          if (w(IDN,k+1,j,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k+1,j,i)/w(IDN,k,j,i);
          if (w(IDN,k,j,i)/w(IDN,k-1,j,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k-1,j,i);
          if (w(IDN,k,j,i)/w(IDN,k+1,j,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k+1,j,i);
        }
      }
    }
  }
  if (drmax > 1.5) return 1;
  else if (drmax < 1.2) return 1;
  return 0;
}
