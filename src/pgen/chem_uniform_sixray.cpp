//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_uniform_sixray.cpp
//! \brief problem generator, uniform chemistry and radiation with six-ray
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // std::find()
#include <cstdio>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

// User defined boundary conditions
void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// Radiation boundary
namespace {
  AthenaArray<Real> G0_iang;
  Real G0, cr_rate;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SixRayBoundaryInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SixRayBoundaryOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x2, SixRayBoundaryInnerX2);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x2, SixRayBoundaryOuterX2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, SixRayBoundaryInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, SixRayBoundaryOuterX3);
  G0 = pin->GetOrAddReal("chem_radiation", "G0", 0.);
  G0_iang.NewAthenaArray(6);
  G0_iang(BoundaryFace::inner_x1) = pin->GetOrAddReal("chem_radiation","G0_inner_x1",G0);
  G0_iang(BoundaryFace::inner_x2) = pin->GetOrAddReal("chem_radiation","G0_inner_x2",G0);
  G0_iang(BoundaryFace::inner_x3) = pin->GetOrAddReal("chem_radiation","G0_inner_x3",G0);
  G0_iang(BoundaryFace::outer_x1) = pin->GetOrAddReal("chem_radiation","G0_outer_x1",G0);
  G0_iang(BoundaryFace::outer_x2) = pin->GetOrAddReal("chem_radiation","G0_outer_x2",G0);
  G0_iang(BoundaryFace::outer_x3) = pin->GetOrAddReal("chem_radiation","G0_outer_x3",G0);
  cr_rate = pin->GetOrAddReal("chem_radiation", "CR", 2e-16);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem of uniform chemistry and radiation
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // read density and radiation field strength
  const Real nH = pin->GetReal("problem", "nH");
  const Real vx = pin->GetOrAddReal("problem", "vx_kms", 0);
  const Real r_init = pin->GetOrAddReal("problem", "r_init", 0.);
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = nH*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // density
        phydro->u(IDN, k, j, i) = nH;
        // velocity, x direction
        phydro->u(IM1, k, j, i) = nH*vx;
        // energy
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*nH*SQR(vx);
        }
      }
    }
  }

  // intialize radiation field
  if (CHEMRADIATION_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifreq=0; ifreq < pchemrad->nfreq; ++ifreq) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              pchemrad->ir(k, j, i, ifreq * pchemrad->nang + iang) = G0_iang(iang);
            }
          }
          if (CHEMISTRY_ENABLED) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              // cr rate
              pchemrad->ir(k, j, i,
                  pscalars->chemnet.index_cr_ * pchemrad->nang + iang) = cr_rate;
            }
          }
        }
      }
    }
    // calculate the average radiation field for output of the initial condition
    pchemrad->pchemradintegrator->CopyToOutput();
  }

  // intialize chemical species
  if (NSPECIES > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            pscalars->s(ispec, k, j, i) = r_init*nH;
            if (CHEMISTRY_ENABLED) {
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec*nH;
              }
            }
          }
        }
      }
    }
  }
  return;
}


void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          pmb->pscalars->s(n,k,j,il-i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,il-i) = 0;
        }
      }
    }
  }
  return;
}


void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,k,jl-j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,jl-j,i) = 0;
        }
      }
    }
  }
  return;
}


void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,kl-k,j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,kl-k,j,i) = 0;
        }
      }
    }
  }
  return;
}


void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          pmb->pscalars->s(n,k,j,iu+i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,iu+i) = 0;
        }
      }
    }
  }
  return;
}


void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,k,ju+j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,ju+j,i) = 0;
        }
      }
    }
  }
  return;
}


void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,ku+k,j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,ku+k,j,i) = 0;
        }
      }
    }
  }
  return;
}
