//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena++ headers
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

//---------------------------------------------------------------------------------------
// Calculate isotropic thermal conduction

void HydroDiffusion::ThermalFlux_iso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                     AthenaArray<Real> *flx)
{
  return;
}


//---------------------------------------------------------------------------------------
// Calculate anisotropic thermal conduction

void HydroDiffusion::ThermalFlux_aniso(const AthenaArray<Real> &p,const AthenaArray<Real> &c,
                     AthenaArray<Real> *flx)
{
  return;
}




//----------------------------------------------------------------------------------------
// constant viscosity

void ConstConduction(HydroDiffusion *phdif, const AthenaArray<Real> &prim,
    const AthenaArray<Real> &bcc, int is, int ie, int js, int je, int ks, int ke)
{
  if (phdif->coeff_kiso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ISO,k,j,i) = phdif->coeff_kiso;
      }
    }
  }
  if (phdif->coeff_kani > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ANI,k,j,i) = phdif->coeff_kani;
      }
    }
  }
  return;
}

