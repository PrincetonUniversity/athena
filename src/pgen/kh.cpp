//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file kh.cpp
//! \brief Problem generator for KH instability.
//!
//! Sets up several different problems:
//!   - iprob=1: slip surface with random perturbations
//!   - iprob=2: tanh profile, with single-mode perturbation (Frank et al. 1996)
//!   - iprob=3: tanh profiles for v and d, SR test problem in Beckwith & Stone (2011)
//!   - iprob=4: tanh profiles for v and d, "Lecoanet" test
//!   - iprob=5: two resolved slip-surfaces with m=2 perturbation for the AMR test

// C headers

// C++ headers
#include <algorithm>  // min, max
#include <cmath>      // log
#include <cstring>    // strcmp()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../defs.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp"

namespace {
Real vflow;
int iprob;
Real PassiveDyeEntropy(MeshBlock *pmb, int iout);
} // namespace

Real threshold;
int RefinementCondition(MeshBlock *pmb);

//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  vflow = pin->GetReal("problem","vflow");
  iprob = pin->GetInteger("problem","iprob");

  if (adaptive) {
    threshold = pin->GetReal("problem", "thr");
    EnrollUserRefinementCondition(RefinementCondition);
  }
  if (iprob == 4 && NSCALARS > 0) {
    AllocateUserHistoryOutput(1);
    EnrollUserHistoryOutput(0, PassiveDyeEntropy, "tot-S");
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Kelvin-Helmholtz test

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::int64_t iseed = -1 - gid;
  Real gm1 = peos->GetGamma() - 1.0;

  //--- iprob=1.  Uniform stream with density ratio "drat" located in region -1/4<y<1/4
  // moving at (-vflow) seperated by two slip-surfaces from background medium with d=1
  // moving at (+vflow), random perturbations.  This is the classic, unresolved K-H test.

  if (iprob == 1) {
    // Read problem parameters
    Real drat = pin->GetReal("problem","drat");
    Real amp = pin->GetReal("problem","amp");
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 1.0;
          phydro->u(IM1,k,j,i) = vflow + amp*(ran2(&iseed) - 0.5);
          phydro->u(IM2,k,j,i) = amp*(ran2(&iseed) - 0.5);
          phydro->u(IM3,k,j,i) = 0.0;
          if (std::abs(pcoord->x2v(j)) < 0.25) {
            phydro->u(IDN,k,j,i) = drat;
            phydro->u(IM1,k,j,i) = -drat*(vflow + amp*(ran2(&iseed) - 0.5));
            phydro->u(IM2,k,j,i) = drat*amp*(ran2(&iseed) - 0.5);
          }
          // Pressure scaled to give a sound speed of 1 with gamma=1.4
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                2.5/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=2. Uniform density medium moving at +/-vflow seperated by a single shear
  // layer with tanh() profile at y=0 with a single mode perturbation, reflecting BCs at
  // top/bottom.  Based on Frank et al., ApJ 460, 777, 1996.

  if (iprob == 2) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem", "amp");
    Real a = 0.02;
    Real sigma = 0.2;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 1.0;
          phydro->u(IM1,k,j,i) = vflow*std::tanh((pcoord->x2v(j))/a);
          phydro->u(IM2,k,j,i) = amp*std::cos(TWO_PI*pcoord->x1v(i))
                                 *std::exp(-(SQR(pcoord->x2v(j)))/SQR(sigma));
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=3.  Test in SR paper (Beckwith & Stone, ApJS 193, 6, 2011).  Gives two
  // resolved shear layers with tanh() profiles for velocity and density located at
  // y = +/- 0.5, density one in middle and 0.01 elsewhere, single mode perturbation.

  if (iprob == 3) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem","amp");
    Real a = 0.01;
    Real sigma = 0.1;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = 0.505 + 0.495
                                 *std::tanh((std::abs(pcoord->x2v(j))-0.5)/a);
          phydro->u(IM1,k,j,i) = vflow*std::tanh((std::abs(pcoord->x2v(j))-0.5)/a);
          phydro->u(IM2,k,j,i) =
              amp*vflow*std::sin(TWO_PI*pcoord->x1v(i))
              *std::exp(-((std::abs(pcoord->x2v(j))-0.5)
                          *(std::abs(pcoord->x2v(j))-0.5))/(sigma*sigma));
          if (pcoord->x2v(j) < 0.0) phydro->u(IM2,k,j,i) *= -1.0;
          phydro->u(IM1,k,j,i) *= phydro->u(IDN,k,j,i);
          phydro->u(IM2,k,j,i) *= phydro->u(IDN,k,j,i);
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                1.0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                               SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  //--- iprob=4.  "Lecoanet" test, resolved shear layers with tanh() profiles for velocity
  // and density located at z1=0.5, z2=1.5 two-mode perturbation for fully periodic BCs

  // To promote symmetry of FP errors about midplanes, rescale z' = z - 1. ; x' = x - 0.5
  // so that domain x1 = [-0.5, 0.5] and x2 = [-1.0, 1.0] is centered about origin
  if (iprob == 4) {
    // Read/set problem parameters
    Real amp = pin->GetReal("problem","amp");
    // unstratified problem is the default
    Real drho_rho0 = pin->GetOrAddReal("problem", "drho_rho0", 0.0);
    // set background vx to nonzero to evolve the KHI in a moving frame
    Real vboost = pin->GetOrAddReal("problem", "vboost", 0.0);
    Real P0 = 10.0;
    Real a = 0.05;
    Real sigma = 0.2;
    // Initial condition's reflect-and-shift symmetry, x1-> x1 + 1/2, x2-> -x2
    // is preserved in new coordinates; hence, the same flow is solved twice in this prob.
    Real z1 = -0.5;  // z1' = z1 - 1.0
    Real z2 = 0.5;   // z2' = z2 - 1.0

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          // Lecoanet (2015) equation 8a)
          Real dens = 1.0 + 0.5*drho_rho0*(std::tanh((pcoord->x2v(j) - z1)/a) -
                                           std::tanh((pcoord->x2v(j) - z2)/a));
          phydro->u(IDN,k,j,i) = dens;

          Real v1 = vflow*(std::tanh((pcoord->x2v(j) - z1)/a)
                           - std::tanh((pcoord->x2v(j) - z2)/a) - 1.0) // 8b)
                    + vboost;
          // Currently, the midpoint approx. is applied in the momenta and energy calc
          phydro->u(IM1,k,j,i) = v1*dens;

          // NOTE ON FLOATING-POINT SHIFT SYMMETRY IN X1:
          // There is no scaling + translation coordinate transformation that would
          // naturally preserve this symmetry when calculating x1 coordinates in
          // floating-point representation. Need to correct for the asymmetry of FP error
          // by modifying the operands.  Centering the domain on x1=0.0 ensures reflective
          // symmetry, x1' -> -x1 NOT shift symmetry, x1' -> x1 + 0.5 (harder guarantee)

          // For example, consider a cell in the right half of the domain with x1v > 0.0,
          // so that shift symmetry should hold w/ another cell's x1v'= -0.5 + x1v < 0.0

          // ISSUE: sin(2*pi*(-0.5+x1v)) != -sin(2*pi*x1v) in floating-point calculations
          // The primary FP issues are twofold: 1) different rounding errors in x1v, x1v'
          // and 2) IEEE-754 merely "recommends" that sin(), cos(), etc. functions are
          // correctly rounded. Note, glibc library doesn't provide correctly-rounded fns

          // 1) Since x1min = -0.5 can be perfectly represented in binary as -2^{-1}:
          // double(x1v')= double(double(-0.5) + double(x1v)) = double(-0.5 + double(x1v))
          // Even if x1v is also a dyadic rational -> has exact finite FP representation:
          // x1v'= double(-0.5 + double(x1v)) = double(-0.5 + x1v) ?= (-0.5 + x1v) exactly

          // Sterbenz's Lemma does not hold for any nx1>4, so cannot guarantee exactness.
          // However, for most nx1 = power of two, the calculation of ALL cell center
          // positions x1v will be exact. For nx1 != 2^n, differences are easily observed.

          // 2) Even if the rounding error of x1v (and hence x1v') is zero, the exact
          // periodicity of trigonometric functions (even after range reduction of input
          // to [-pi/4, pi/4], e.g.) is NOT guaranteed:
          // sin(2*pi*(-0.5+x1v)) = sin(-pi + 2*pi*x1v) != -sin(2*pi*x1v)

          // WORKAROUND: Average inexact sin() with -sin() sample on opposite x1-half of
          // domain The assumption of periodic domain with x1min=-0.5 and x1max=0.5 is
          // hardcoded here (v2 is the only quantity in the IC with x1 dependence)

          Real ave_sine = std::sin(TWO_PI*pcoord->x1v(i));
          if (pcoord->x1v(i) > 0.0) {
            ave_sine -= std::sin(TWO_PI*(-0.5 + pcoord->x1v(i)));
          } else {
            ave_sine -= std::sin(TWO_PI*(0.5 + pcoord->x1v(i)));
          }
          ave_sine /= 2.0;

          // translated x1= x - 1/2 relative to Lecoanet (2015) shifts sine function by pi
          // (half-period) and introduces U_z sign change:
          Real v2 = -amp*ave_sine
                    *(std::exp(-(SQR(pcoord->x2v(j) - z1))/(sigma*sigma)) +
                      std::exp(-(SQR(pcoord->x2v(j) - z2))/(sigma*sigma))); // 8c), mod.
          phydro->u(IM2,k,j,i) = v2*dens;

          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = P0/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                                 + SQR(phydro->u(IM2,k,j,i))
                                                 + SQR(phydro->u(IM3,k,j,i)) )
                                   /phydro->u(IDN,k,j,i);
          }
          // color concentration of passive scalar
          if (NSCALARS > 0) {
            Real concentration = 0.5*(std::tanh((pcoord->x2v(j) - z2)/a)  // 8e)
                                      - std::tanh((pcoord->x2v(j) - z1)/a) + 2.0);
            // uniformly fill all scalar species to have equal concentration
            constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0;
            for (int n=0; n<NSCALARS; ++n) {
              pscalars->s(n,k,j,i) = 1.0/scalar_norm*concentration*phydro->u(IDN,k,j,i);
            }
          }
        }
      }
    }
    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem", "b0");
      b0 = b0/std::sqrt(4.0*(PI));
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0*std::tanh((std::abs(pcoord->x2v(j)) - 0.5)/a);
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*SQR(pfield->b.x1f(k,j,i));
            }
          }
        }
      }
    }
  }

  //--- iprob=5. Uniform stream with density ratio "drat" located in region -1/4<y<1/4
  // moving at (-vflow) seperated by two resolved slip-surfaces from background medium
  // with d=1 moving at (+vflow), with m=2 perturbation, for the AMR test.

  if (iprob == 5) {
    // Read problem parameters
    Real a = pin->GetReal("problem","a");
    Real sigma = pin->GetReal("problem","sigma");
    Real drat = pin->GetReal("problem","drat");
    Real amp = pin->GetReal("problem","amp");
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real w=(std::tanh((std::abs(pcoord->x2v(j))-0.25)/a)+1.0)*0.5;
          phydro->u(IDN,k,j,i) = w+(1.0-w)*drat;
          phydro->u(IM1,k,j,i) = w*vflow-(1.0-w)*vflow*drat;
          phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*amp
                                 * std::sin(2.0*TWO_PI*pcoord->x1v(i))
                                 * std::exp(-SQR(std::abs(pcoord->x2v(j))-0.25)
                                            /(sigma*sigma));
          phydro->u(IM3,k,j,i) = 0.0;
          // Pressure scaled to give a sound speed of 1 with gamma=1.4
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) =
                2.5/gm1 + 0.25*(SQR(phydro->u(IM1,k,j,i)) +
                                SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
          }
        }
      }
    }

    // initialize uniform interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }
  }

  return;
}


// refinement condition: velocity gradient

int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vgmax = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real vgy = std::abs(w(IVY,k,j,i+1) - w(IVY,k,j,i-1))*0.5;
        Real vgx = std::abs(w(IVX,k,j+1,i) - w(IVX,k,j-1,i))*0.5;
        Real vg  = std::sqrt(vgx*vgx + vgy*vgy);
        if (vg > vgmax) vgmax = vg;
      }
    }
  }
  if (vgmax > threshold) return 1;
  if (vgmax < 0.5*threshold) return -1;
  return 0;
}

namespace {
Real PassiveDyeEntropy(MeshBlock *pmb, int iout) {
  Real total_entropy = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &r = pmb->pscalars->r;
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=is; i<=ie; i++) {
        // no loop over NSCALARS; hardcode assumption that NSCALARS=1
        Real specific_entropy = -r(0,k,j,i)*std::log(r(0,k,j,i));
        total_entropy += volume(i)*w(IDN,k,j,i)*specific_entropy;  // Lecoanet (2016) eq 5
      }
    }
  }
  return total_entropy;
}
} // namespace
