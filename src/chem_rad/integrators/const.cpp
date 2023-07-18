//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file const.cpp
//! \brief implementation of radiation integrators: constant radiation

// this class header
#include "rad_integrators.hpp"

// Athena++ headers
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../chem_rad.hpp"

//----------------------------------------------------------------------------------------
//! constructor, for constant radiation integrator
ChemRadIntegrator::ChemRadIntegrator(ChemRadiation *pchemrad, ParameterInput *pin)
    : col(0, 0, 0, 0),
    col_bvar(pchemrad->pmy_block, &col) {
  pmy_mb = pchemrad->pmy_block;
  pmy_rad = pchemrad;
}

//----------------------------------------------------------------------------------------
//! destructor
ChemRadIntegrator::~ChemRadIntegrator() {}

//----------------------------------------------------------------------------------------
//! \fn void ChemRadIntegrator::CopyToOutput()
//! \brief average radiation field over all angles and copy values to output
void ChemRadIntegrator::CopyToOutput() {
  int is = pmy_mb->is;
  int ie = pmy_mb->ie;
  int js = pmy_mb->js;
  int je = pmy_mb->je;
  int ks = pmy_mb->ks;
  int ke = pmy_mb->ke;
  int jl, ju, kl, ku;
  if (js == 0 && je == 0) {
    jl = ju = 0;
  } else {
    jl = js-NGHOST;
    ju = je+NGHOST;
  }
  if (ks == 0 && ke == 0) {
    kl = ku = 0;
  } else {
    kl = ks-NGHOST;
    ku = ke+NGHOST;
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        for (int ifreq=0; ifreq < pmy_rad->nfreq; ++ifreq) {
          pmy_rad->ir_avg(ifreq, k, j, i) =
            pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemRadIntegrator::UpdateRadiation()
//! \brief update radiation field
void ChemRadIntegrator::UpdateRadiation() {}

void ChemRadIntegrator::GetColMB(BoundaryFace direction) {}
void ChemRadIntegrator::UpdateCol(BoundaryFace direction) {}
