//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file magnoh.cpp
//  \brief Magnetized Noh with perturbation in B_phi. Authored by A. Beresnyak
//========================================================================================

// C++ headers
#include <cmath>      // sqrt()

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
static Real nu_iso, eta_ohm;

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  nu_iso = pin->GetOrAddReal("problem", "nu_iso", 0.0);
  eta_ohm = pin->GetOrAddReal("problem", "eta_ohm", 0.0);

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

  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for zpinch problem
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;

  // initialize conserved variables
  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie+1; i++) {
        // Volume centered coordinates
        Real x1=pcoord->x1v(i);
        Real x2=pcoord->x2v(j);
        Real rad = std::sqrt(SQR(x1)+SQR(x2));
        if (rad==0.0) rad=3.0/100000;
        // Real phi=atan2(x2,x1); // unused variable

        // Face centered coordinates
        Real  x1f=pcoord->x1f(i);
        Real  x2f=pcoord->x2f(j);
        Real   rad_f1 = std::sqrt(SQR(x1f)+SQR(x2 ));
        Real   rad_f2 = std::sqrt(SQR(x1 )+SQR(x2f));
        if(rad_f1==0.0) rad_f1=3.0/100000;
        if(rad_f2==0.0) rad_f2=3.0/100000;
        Real phi_f1=atan2(x2 ,x1f);
        Real phi_f2=atan2(x2f,x1 );
        // Volume-centered variables
        Real rho=rho0*pow(rad,alpha);
        Real P  = P0*pow(rad,2*beta);
        Real V1 = vr*x1/rad;
        Real V2 = vr*x2/rad;
        Real V3 = 0.0;
        Real B3f = bz*pow(rad,beta);

        Real B1f = bphi0*pow(rad_f1,beta-1)*(x2+x1f*perturb*cos(mphi*phi_f1));
        Real B2f = bphi0*pow(rad_f2,beta-1)*(-x1+x2f*perturb*cos(mphi*phi_f2));
        Real EN= P/gm1 + 0.5*rho*(SQR(V1)+SQR(V2)+SQR(V3));

        if (i!=ie+1 && j!=je+1 && k!=ke+1) {
          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = rho*V1;
          phydro->u(IM2,k,j,i) = rho*V2;
          phydro->u(IM3,k,j,i) = rho*V3;
          phydro->u(IEN,k,j,i) = EN;
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          if (j!=je+1 && k!=ke+1)
            pfield->b.x1f(k,j,i) = B1f;
          if (i!=ie+1 && k!=ke+1)
            pfield->b.x2f(k,j,i) = B2f;
          if (i!=ie+1 && j!=je+1)
            pfield->b.x3f(k,j,i) = B3f;
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) +=  0.5*0.25*(
              SQR(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))
              + SQR(pfield->b.x2f(k,j,i)+pfield->b.x2f(k,j+1,i))
              + SQR(pfield->b.x3f(k,j,i)+pfield->b.x3f(k+1,j,i)));
        }
      }
    }
  }
  return;
}
