//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file collapse.cpp
//! \brief Problem generator for collapse of a Bonnor-Ebert like sphere with AMR or SMR

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"

#if SELF_GRAVITY_ENABLED != 2
#error "This problem generator requires Multigrid gravity solver."
#endif

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


namespace {


// dimension-less constants
constexpr Real four_pi_G = 1.0;
constexpr Real rc = 6.45; // the BE radius in the normalized unit system
constexpr Real rcsq = 26.0 / 3.0;      // the parameter of the BE profile
constexpr Real bemass = 197.561;       // the total mass of the critical BE sphere

// dimensional constants
constexpr Real pi   = M_PI;
constexpr Real cs10 = 1.9e4;        // sound speed at 10K, cm / s
constexpr Real msun = 1.9891e33;    // solar mass, g
constexpr Real pc   = 3.0857000e18; // parsec, cm
constexpr Real au   = 1.4959787e13; // astronomical unit, cm
constexpr Real yr   = 3.15569e7;    // year, s
constexpr Real G    = 6.67259e-8;   // gravitational constant, dyn cm^2 g^-2

// units in cgs
Real m0, v0, t0, l0, rho0, gauss;

// parameters and derivatives
Real mass, temp, f, rhocrit, omega, bz, mu, amp;

// AMR parameter
Real njeans; // Real is used intentionally

Real totalm;
} // namespace

// Mask the density outside the initial sphere
void SourceMask(AthenaArray<Real> &src, int is, int ie, int js, int je,
                int ks, int ke, const MGCoordinates &coord) {
  const Real rc2 = rc*rc;
  for (int k=ks; k<=ke; ++k) {
    Real z = coord.x3v(k);
    for (int j=js; j<=je; ++j) {
      Real y = coord.x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x = coord.x1v(i);
        Real r2 = x*x + y*y + z*z;
        if (r2 > rc2)
          src(k, j, i) = 0.0;
      }
    }
  }
  return;
}


// AMR refinement condition
int JeansCondition(MeshBlock *pmb) {
  Real njmin = 1e300;
  const Real dx = pmb->pcoord->dx1f(0); // assuming uniform cubic cells
  if (MAGNETIC_FIELDS_ENABLED) {
    if (NON_BAROTROPIC_EOS) {
      const Real gamma = pmb->peos->GetGamma();
      const Real fac = 2.0 * pi / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real v = std::sqrt(gamma * pmb->phydro->w(IPR,k,j,i)
                              / pmb->phydro->w(IDN,k,j,i))+
                   + std::sqrt((SQR(pmb->pfield->bcc(IB1,k,j,i))
                              + SQR(pmb->pfield->bcc(IB2,k,j,i))
                              + SQR(pmb->pfield->bcc(IB3,k,j,i)))
                              / pmb->phydro->w(IDN,k,j,i));
            Real nj = v / std::sqrt(pmb->phydro->w(IDN,k,j,i));
            njmin = std::min(njmin, nj);
          }
        }
      }
      njmin *= fac;
    } else {
      const Real cs = pmb->peos->GetIsoSoundSpeed();
      const Real fac = 2.0 * pi / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real v = cs + std::sqrt((SQR(pmb->pfield->bcc(IB1,k,j,i))
                                   + SQR(pmb->pfield->bcc(IB2,k,j,i))
                                   + SQR(pmb->pfield->bcc(IB3,k,j,i)))
                                   / pmb->phydro->w(IDN,k,j,i));
            Real nj = v / std::sqrt(pmb->phydro->w(IDN,k,j,i));
            njmin = std::min(njmin, nj);
          }
        }
      }
      njmin *= fac;
    }
  } else {
    if (NON_BAROTROPIC_EOS) {
      const Real gamma = pmb->peos->GetGamma();
      const Real fac = 2.0 * pi * std::sqrt(gamma) / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real nj = std::sqrt(pmb->phydro->w(IPR,k,j,i))
                              / pmb->phydro->w(IDN,k,j,i);
            njmin = std::min(njmin, nj);
          }
        }
      }
      njmin *= fac;
    } else {
      const Real cs = pmb->peos->GetIsoSoundSpeed();
      const Real fac = 2.0 * pi / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real nj = cs / std::sqrt(pmb->phydro->w(IDN,k,j,i));
            njmin = std::min(njmin, nj);
          }
        }
      }
      njmin *= fac;
    }
  }
  if (njmin < njeans)
    return 1;
  if (njmin > njeans * 2.5)
    return -1;
  return 0;
}

// Approximated Bonnor-Ebert profile
// Tomida 2011, PhD Thesis
Real BEProfile(Real r) {
  return std::pow(1.0+r*r/rcsq, -1.5);
}


void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  const Real gm1 = pmb->peos->GetGamma() - 1.0;
  const Real igm1 = 1.0 / gm1;

  // Fixed boundary condition outside the initial core
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real z = pmb->pcoord->x3v(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real y = pmb->pcoord->x2v(j);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x = pmb->pcoord->x1v(i);
        Real r = std::sqrt(SQR(x) + SQR(y) + SQR(z));
        if (r > rc) {
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
        }
      }
    }
  }

  // Set the internal energy to follow the barotropic relation
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real ke = 0.5 / cons(IDN,k,j,i)
                  * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
          Real me = 0.5*(SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)));
          Real te = igm1 * cons(IDN,k,j,i)
                  * std::sqrt(1.0+std::pow(cons(IDN,k,j,i)/rhocrit, 2.0*gm1));
          cons(IEN,k,j,i) = te + ke + me;
        }
      }
    }
  } else {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real ke = 0.5 / cons(IDN,k,j,i)
                  * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
          Real te = igm1 * cons(IDN,k,j,i)
                  * std::sqrt(1.0+std::pow(cons(IDN,k,j,i)/rhocrit, 2.0*gm1));
          cons(IEN,k,j,i) = te + ke;
        }
      }
    }
  }

  return;
}




//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  SetFourPiG(four_pi_G); // 4piG = 1.0
  mass = pin->GetReal("problem", "mass"); // solar mass
  temp = pin->GetReal("problem", "temperature");
  f = pin->GetReal("problem", "f"); // Density enhancement factor; f = 1 is critical
  amp = pin->GetOrAddReal("problem", "amp", 0.0); // perturbation amplitude
  mu = pin->GetOrAddReal("problem", "mu", 0.0); // micro gauss
  m0 = mass * msun / (bemass*f); // total mass = 1.0
  v0 = cs10 * std::sqrt(temp/10.0); // cs at 10K = 1.0
  rho0 = (v0*v0*v0*v0*v0*v0) / (m0*m0) /(64.0*pi*pi*pi*G*G*G);
  t0 = 1.0/std::sqrt(4.0*pi*G*rho0); // sqrt(1/4pi G rho0) = 1.0
  l0 = v0 * t0;
  gauss = std::sqrt(rho0*v0*v0*4.0*pi);
  rhocrit = pin->GetReal("problem", "rhocrit") / rho0;
  totalm = 0.0;
  Real tff = std::sqrt(3.0/8.0/f)*pi;
  Real omegatff = pin->GetOrAddReal("problem", "omegatff", 0.0);
  omega = omegatff/tff;
  Real mucrit1 = 0.53/(3.0*pi)*std::sqrt(5.0/G);
  Real mucrit2 = 1.0/(2.0*pi*std::sqrt(G));
  bz = mass*msun/mucrit1/mu/pi/SQR(rc*l0)/gauss;
  Real bzug = bz*gauss*1e6;
  if (Globals::my_rank == 0 && ncycle == 0) {
    std::cout << std::endl
      << "---  Dimensional parameters of the simulation  ---" << std::endl
      << "Total mass          : " << mass      << " \t\t[Msun]" << std::endl
      << "Initial temperature : " << temp      << " \t\t[K]" << std::endl
      << "Sound speed         : " << v0        << " \t\t[cm s^-1]" << std::endl
      << "Central density     : " << rho0*f    << " \t[g cm^-3]" << std::endl
      << "Cloud radius        : " << rc*l0/au  << " \t\t[au]" << std::endl
      << "Free fall time      : " << tff*t0/yr << " \t\t[yr]" << std::endl
      << "Angular velocity    : " << omega/t0  << " \t[s^-1]" << std::endl
      << "Angular velocity    : " << omega/t0*yr << " \t[yr^-1]" << std::endl
      << "Magnetic field      : " << bzug      << " \t\t[uGauss]" << std::endl
      << "Density Enhancement : " << f         << std::endl << std::endl
      << "---   Normalization Units of the simulation    ---" << std::endl
      << "Mass                : " << m0        << " \t[g]" << std::endl
      << "Mass                : " << m0/msun   << " \t[Msun]" << std::endl
      << "Length              : " << l0        << " \t[cm]" << std::endl
      << "Length              : " << l0/au     << " \t\t[au]" << std::endl
      << "Length              : " << l0/pc     << " \t[pc]" << std::endl
      << "Time                : " << t0        << " \t[s]" << std::endl
      << "Time                : " << t0/yr     << " \t\t[yr]" << std::endl
      << "Velocity            : " << v0        << " \t\t[cm s^-1]" << std::endl
      << "Density             : " << rho0      << " \t[g cm^-3]" << std::endl
      << "Magnetic field      : " << gauss     << " \t[Gauss]" << std::endl << std::endl
      << "--- Dimensionless parameters of the simulation ---" << std::endl
      << "Total mass          : " << bemass*f  << std::endl
      << "Sound speed at " << temp << " K : "  << 1.0 << std::endl
      << "Central density     : " << 1.0       << std::endl
      << "Cloud radius        : " << rc        << std::endl
      << "Free fall time      : " << tff       << std::endl
      << "m=2 perturbation    : " << amp       << std::endl
      << "Omega * tff         : " << omegatff << std::endl
      << "Mass-to-flux ratio  : " << mu        << std::endl << std::endl;
  }

  EnrollUserMGGravitySourceMaskFunction(SourceMask);

  if (NON_BAROTROPIC_EOS)
    EnrollUserExplicitSourceFunction(Cooling);

  if (adaptive) {
    njeans = pin->GetReal("problem","njeans");
    EnrollUserRefinementCondition(JeansCondition);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real igm1 = 1.0 / (peos->GetGamma() - 1.0);
  Real dx = pcoord->dx1f(is);

  for (int k=ks; k<=ke; ++k) {
    Real z = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real y = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x = pcoord->x1v(i);
        Real r = std::sqrt(SQR(x) + SQR(y) + SQR(z));
        r = std::min(r, rc); // pressure confinement - constant beyond the cloud radius
        phydro->u(IDN,k,j,i) = f * BEProfile(r)
                             * (1.0+amp*SQR(r)/SQR(rc)*cos(2.0*std::atan2(y,x)));
        if (r < rc) {
          phydro->u(IM1,k,j,i) =  phydro->u(IDN,k,j,i)*omega*y;
          phydro->u(IM2,k,j,i) = -phydro->u(IDN,k,j,i)*omega*x;
          phydro->u(IM3,k,j,i) = 0.0;
        } else {
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
        }
        if (NON_BAROTROPIC_EOS)
          phydro->u(IEN,k,j,i) = igm1 * phydro->u(IDN,k,j,i) // c_s = 1
                               + 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                               + SQR(phydro->u(IM3,k,j,i)))/SQR(phydro->u(IDN,k,j,i));
      }
    }
  }
  // initialize interface B, uniform
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = bz;
        }
      }
    }
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) += 0.5*SQR(bz);
          }
        }
      }
    }
  }
}

