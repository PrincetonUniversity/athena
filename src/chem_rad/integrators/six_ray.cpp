//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file six_ray.cpp
//! \brief implementation of radiation integrators: six_ray

//c++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream

// Athena++ headers
#include "../../bvals/bvals.hpp"
#include "../../chemistry/utils/shielding.hpp"
#include "../../chemistry/utils/thermo.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../../units/units.hpp"
#include "../chem_rad.hpp"

// this class header
#include "rad_integrators.hpp"

namespace {
  AthenaArray<Real> G0_iang; // diffuse radiation field strength in Draine 1987 unit
  Real G0, cr_rate; // cosmic ray ionisation rate
  Real f_cell, f_prev; // fraction of the column in the cell that goes to shielding
  Real lunit; // length unit in cm
}

//----------------------------------------------------------------------------------------
//! constructor, for six_ray radiation integrator

ChemRadIntegrator::ChemRadIntegrator(ChemRadiation *pchemrad, ParameterInput *pin) :
    col(6, pchemrad->pmy_block->ncells3, pchemrad->pmy_block->ncells2,
        pchemrad->pmy_block->ncells1,  pchemrad->pmy_block->pscalars->chemnet.n_cols_),
    col_bvar(pchemrad->pmy_block, &col) {
  pmy_mb = pchemrad->pmy_block;
  pmy_rad = pchemrad;
  G0 = pin->GetOrAddReal("chem_radiation", "G0", 0.);
  G0_iang.NewAthenaArray(6);
  G0_iang(BoundaryFace::inner_x1) = pin->GetOrAddReal("chem_radiation","G0_inner_x1",G0);
  G0_iang(BoundaryFace::inner_x2) = pin->GetOrAddReal("chem_radiation","G0_inner_x2",G0);
  G0_iang(BoundaryFace::inner_x3) = pin->GetOrAddReal("chem_radiation","G0_inner_x3",G0);
  G0_iang(BoundaryFace::outer_x1) = pin->GetOrAddReal("chem_radiation","G0_outer_x1",G0);
  G0_iang(BoundaryFace::outer_x2) = pin->GetOrAddReal("chem_radiation","G0_outer_x2",G0);
  G0_iang(BoundaryFace::outer_x3) = pin->GetOrAddReal("chem_radiation","G0_outer_x3",G0);
  cr_rate = pin->GetOrAddReal("chem_radiation", "CR", 2e-16);
  f_cell = pin->GetOrAddReal("chem_radiation", "shielding_fraction_cell", 0.5);
  f_prev = 1. - f_cell;
  std::stringstream msg; // error message
  if (pmy_rad->nang != 6) {
    msg << "### FATAL ERROR in ChemRadIntegrator constructor [ChemRadIntegrator]"
      << std::endl
      << "Six-ray scheme with nang != 6." << std::endl;
    ATHENA_ERROR(msg);
  }
  if (CHEMISTRY_ENABLED) {
    pmy_chemnet = &pmy_mb->pscalars->chemnet;
    lunit = pmy_mb->pmy_mesh->punit->code_length_cgs;
    ncol = pmy_chemnet->n_cols_;
    // allocate array for column density
    // enroll SixRayBoundaryVariable object
    // to enable functions such as yCopyVariableBufferSameProcess()
    col_bvar.bvar_index = pmy_mb->pbval->bvars.size();
    pmy_mb->pbval->bvars.push_back(&col_bvar);
    if (DEBUG) {
      col_avg.NewAthenaArray(ncol, pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1);
      col_Htot.NewAthenaArray(6, pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1);
      col_H2.NewAthenaArray(6, pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1);
      col_CO.NewAthenaArray(6, pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1);
      col_C.NewAthenaArray(6, pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1);
    }
  }
}


//----------------------------------------------------------------------------------------
//! destructor
ChemRadIntegrator::~ChemRadIntegrator() {}

//----------------------------------------------------------------------------------------
//! \fn void ChemRadIntegrator::CopyToOutput()
//! \brief average radiation field over all angles and copy values to output

void ChemRadIntegrator::CopyToOutput() {
  const int is = pmy_mb->is;
  const int js = pmy_mb->js;
  const int ks = pmy_mb->ks;
  const int ie = pmy_mb->ie;
  const int je = pmy_mb->je;
  const int ke = pmy_mb->ke;
  int iang_arr[6] =
  {BoundaryFace::inner_x1, BoundaryFace::inner_x2, BoundaryFace::inner_x3,
   BoundaryFace::outer_x1, BoundaryFace::outer_x2, BoundaryFace::outer_x3};
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
        for (int iang=0; iang < 6; iang++) {
          // column densities
          if (CHEMISTRY_ENABLED) {
            if (DEBUG) {
            col_Htot(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNHtot_);
            col_H2(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNH2_);
            col_CO(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNCO_);
            col_C(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNC_);
            //  angle-averaged column densities
            for (int icol=0; icol<ncol; icol++) {
              if (iang == 0) {
                col_avg(icol, k, j, i) = 0;
              }
              col_avg(icol, k, j, i) += col(iang, k, j, i, icol)/6.;
            }
            }
          }
          // radiation
          for (int ifreq=0; ifreq < pmy_rad->nfreq; ++ifreq) {
            if (iang == 0) {
              pmy_rad->ir_avg(ifreq, k, j, i) = 0;
            }
            pmy_rad->ir_avg(ifreq, k, j, i) +=
              pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang+iang)/6.;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemRadIntegrator::UpdateRadiation()
//! \brief calcuate total column and update radiation

void ChemRadIntegrator::UpdateRadiation() {
  const Real Zd = pmy_chemnet->zdg_;
  const Real bH2 = 3.0e5; // H2 velocity dispersion
  const int iph_H2 = ChemNetwork::iph_H2_;
  const int iph_CO = ChemNetwork::iph_CO_;
  const int iph_C = ChemNetwork::iph_C_;
  const int iPE = pmy_chemnet->index_gpe_;
  const int iCR = pmy_chemnet->index_cr_;
  const Real sigmaPE = Thermo::sigmaPE_;
  const Real NH0_CR = 9.35e20;
  Real NH, AV, NCO, NC, NH2, fs_CO, fs_H2, fs_CI;
  for (int direction=0; direction < 6; direction++) {
    int iang = direction;
    for (int k=pmy_mb->ks; k<=pmy_mb->ke; ++k) {
      for (int j=pmy_mb->js; j<=pmy_mb->je; ++j) {
        for (int i=pmy_mb->is; i<=pmy_mb->ie; ++i) {
          NH = col(direction, k, j, i, pmy_chemnet->iNHtot_);
          NH2 = col(direction, k, j, i,  pmy_chemnet->iNH2_);
          NCO = col(direction, k, j, i,  pmy_chemnet->iNCO_);
          NC = col(direction, k, j, i,  pmy_chemnet->iNC_);
          AV = NH * Zd / 1.87e21;
          // photo-reactions
          for (int ifreq=0; ifreq < pmy_rad->nfreq-2; ++ifreq) {
            pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang+iang) = G0_iang(iang)
              * std::exp( -ChemNetwork::kph_avfac_[ifreq] * AV );
          }
          // H2, CO and CI self-shielding
          fs_CO = Shielding::fShield_CO_V09(NCO, NH2);
          fs_H2 = Shielding::fShield_H2(NH2, bH2);
          fs_CI = Shielding::fShield_C(NC, NH2);
          pmy_rad->ir(k, j, i, iph_H2 * pmy_rad->nang+iang) *=
            fs_H2;
          pmy_rad->ir(k, j, i, iph_CO * pmy_rad->nang+iang) *=
            fs_CO;
          pmy_rad->ir(k, j, i, iph_C * pmy_rad->nang+iang) *=
            fs_CI;
          // GPE
          pmy_rad->ir(k, j, i, iPE * pmy_rad->nang+iang) = G0_iang(iang)
            *  std::exp(-NH * sigmaPE * Zd);
          // CR
          if (pmy_chemnet->is_cr_shielding_) {
            if (NH <= NH0_CR) {
              pmy_rad->ir(k, j, i, iCR * pmy_rad->nang+iang)
                = cr_rate;
            } else {
              pmy_rad->ir(k, j, i, iCR * pmy_rad->nang+iang)
                = cr_rate * (NH0_CR/NH);
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GetColMB(BoundaryFace direction)
//! \brief calculate column densities within the meshblock
//!
//! direction: 0:+x, 1:-x, 2:+y, 3:-y, 4:+y, 5:-z (bvals/bvals_interfaces.hpp)

void ChemRadIntegrator::GetColMB(BoundaryFace direction) {
  const int iH2 = pmy_chemnet->iH2_;
  const int iCO = pmy_chemnet->iCO_;
  const int iCplus = pmy_chemnet->iCplus_;
  const int iHCOplus = pmy_chemnet->iHCOplus_;
  const int iCHx = pmy_chemnet->iCHx_;
  const Real xCtot =  pmy_chemnet->xC_;
  const int is = pmy_mb->is;
  const int js = pmy_mb->js;
  const int ks = pmy_mb->ks;
  const int ie = pmy_mb->ie;
  const int je = pmy_mb->je;
  const int ke = pmy_mb->ke;
  Real NHtot_cell, NHtot_cell_prev;
  Real xCI, xCI_prev;
  std::stringstream msg; // error message
  if (direction == BoundaryFace::inner_x1) {
    // +x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx1f(i)
            * lunit;
          if (i == is) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell;
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI;
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j, i-1)
              * pmy_mb->pcoord->dx1f(i-1) * lunit;
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell
              +  NHtot_cell_prev*f_prev
              + col(direction, k, j, i-1, pmy_chemnet->iNHtot_);
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iH2, k, j, i-1)
              + col(direction, k, j, i-1, pmy_chemnet->iNH2_);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iCO, k, j, i-1)
              + col(direction, k, j, i-1, pmy_chemnet->iNCO_);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j, i-1)
                  - pmy_mb->pscalars->r(iCplus, k, j, i-1)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i-1)
                  - pmy_mb->pscalars->r(iCHx, k, j, i-1);
            if (xCI_prev < 0.0) {
              xCI_prev = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI
              + NHtot_cell_prev*f_prev * xCI_prev
              + col(direction, k, j, i-1, pmy_chemnet->iNC_);
          }
        }
      }
    }
  } else if (direction == BoundaryFace::outer_x1) {
    // -x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=ie; i>=is; --i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx1f(i)
            * lunit;
          if (i == ie) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell;
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI;
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j, i+1)
              * pmy_mb->pcoord->dx1f(i+1) * lunit;
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell
              + NHtot_cell_prev*f_prev
              + col(direction, k, j, i+1, pmy_chemnet->iNHtot_);
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iH2, k, j, i+1)
              + col(direction, k, j, i+1, pmy_chemnet->iNH2_);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iCO, k, j, i+1)
              + col(direction, k, j, i+1, pmy_chemnet->iNCO_);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j, i+1)
                  - pmy_mb->pscalars->r(iCplus, k, j, i+1)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i+1)
                  - pmy_mb->pscalars->r(iCHx, k, j, i+1);
            if (xCI_prev < 0.0) {
              xCI_prev = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI
              + NHtot_cell_prev*f_prev * xCI_prev
              + col(direction, k, j, i+1, pmy_chemnet->iNC_);
          }
        }
      }
    }
  } else if (direction == BoundaryFace::inner_x2) {
    // +y
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx2f(j)
            * lunit;
          if (j == js) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell;
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI;
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j-1, i)
              * pmy_mb->pcoord->dx2f(j-1) * lunit;
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell
              + NHtot_cell_prev*f_prev
              + col(direction, k, j-1, i, pmy_chemnet->iNHtot_);
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iH2, k, j-1, i)
              + col(direction, k, j-1, i, pmy_chemnet->iNH2_);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iCO, k, j-1, i)
              + col(direction, k, j-1, i, pmy_chemnet->iNCO_);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j-1, i)
                  - pmy_mb->pscalars->r(iCplus, k, j-1, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j-1, i)
                  - pmy_mb->pscalars->r(iCHx, k, j-1, i);
            if (xCI_prev < 0.0) {
              xCI_prev = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI
              + NHtot_cell_prev*f_prev * xCI_prev
              + col(direction, k, j-1, i, pmy_chemnet->iNC_);
          }
        }
      }
    }
  } else if (direction == BoundaryFace::outer_x2) {
    // -y
    for (int k=ks; k<=ke; ++k) {
      for (int j=je; j>=js; --j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx2f(j)
            * lunit;
          if (j == je) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell;
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI;
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j+1, i)
              * pmy_mb->pcoord->dx2f(j+1) * lunit;
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell
              + NHtot_cell_prev*f_prev
              + col(direction, k, j+1, i, pmy_chemnet->iNHtot_);
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iH2, k, j+1, i)
              + col(direction, k, j+1, i, pmy_chemnet->iNH2_);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iCO, k, j+1, i)
              + col(direction, k, j+1, i, pmy_chemnet->iNCO_);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j+1, i)
                  - pmy_mb->pscalars->r(iCplus, k, j+1, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j+1, i)
                  - pmy_mb->pscalars->r(iCHx, k, j+1, i);
            if (xCI_prev < 0.0) {
              xCI_prev = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI
              + NHtot_cell_prev*f_prev * xCI_prev
              + col(direction, k, j+1, i, pmy_chemnet->iNC_);
          }
        }
      }
    }
  } else if (direction ==  BoundaryFace::inner_x3) {
    // +z
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx3f(k)
            * lunit;
          if (k == ks) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell;
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI;
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k-1, j, i)
              * pmy_mb->pcoord->dx3f(k-1) * lunit;
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell
              + NHtot_cell_prev*f_prev
              + col(direction, k-1, j, i, pmy_chemnet->iNHtot_);
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iH2, k-1, j, i)
              + col(direction, k-1, j, i, pmy_chemnet->iNH2_);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iCO, k-1, j, i)
              + col(direction, k-1, j, i, pmy_chemnet->iNCO_);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k-1, j, i)
                  - pmy_mb->pscalars->r(iCplus, k-1, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k-1, j, i)
                  - pmy_mb->pscalars->r(iCHx, k-1, j, i);
            if (xCI_prev < 0.0) {
              xCI_prev = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI
              + NHtot_cell_prev*f_prev * xCI_prev
              + col(direction, k-1, j, i, pmy_chemnet->iNC_);
          }
        }
      }
    }
  } else if (direction == BoundaryFace::outer_x3) {
    // -z
    for (int k=ke; k>=ks; --k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx3f(k)
            * lunit;
          if (k == ke) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell;
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI;
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k+1, j, i)
              * pmy_mb->pcoord->dx3f(k+1) * lunit;
            col(direction, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell*f_cell
              + NHtot_cell_prev*f_prev
              + col(direction, k+1, j, i, pmy_chemnet->iNHtot_);
            col(direction, k, j, i, pmy_chemnet->iNH2_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iH2, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iH2, k+1, j, i)
              + col(direction, k+1, j, i, pmy_chemnet->iNH2_);
            col(direction, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell*f_cell * pmy_mb->pscalars->r(iCO, k, j, i)
              + NHtot_cell_prev*f_prev * pmy_mb->pscalars->r(iCO, k+1, j, i)
              + col(direction, k+1, j, i, pmy_chemnet->iNCO_);
            xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, i)
                  - pmy_mb->pscalars->r(iCplus, k, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i)
                  - pmy_mb->pscalars->r(iCHx, k, j, i);
            if (xCI < 0.0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k+1, j, i)
                  - pmy_mb->pscalars->r(iCplus, k+1, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k+1, j, i)
                  - pmy_mb->pscalars->r(iCHx, k+1, j, i);
            if (xCI_prev < 0.0) {
              xCI_prev = 0.;
            }
            col(direction, k, j, i, pmy_chemnet->iNC_) =
              NHtot_cell*f_cell * xCI
              + NHtot_cell_prev*f_prev * xCI_prev
              + col(direction, k+1, j, i, pmy_chemnet->iNC_);
          }
        }
      }
    }
  } else {
    msg << "### FATAL ERROR in ChemRadIntegrator six_ray [GetColMB]" << std::endl
      << "direction {0,1,2,3,4,5}:" << direction << " unknown." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemRadIntegrator::UpdateCol(BoundaryFace direction)
//! \brief update the columns after receiving boundary

void ChemRadIntegrator::UpdateCol(BoundaryFace direction) {
  const int iH2 = pmy_chemnet->iH2_;
  const int iCO = pmy_chemnet->iCO_;
  const int iCplus = pmy_chemnet->iCplus_;
  const int iHCOplus = pmy_chemnet->iHCOplus_;
  const int iCHx = pmy_chemnet->iCHx_;
  const Real xCtot =  pmy_chemnet->xC_;
  const int is = pmy_mb->is;
  const int js = pmy_mb->js;
  const int ks = pmy_mb->ks;
  const int ie = pmy_mb->ie;
  const int je = pmy_mb->je;
  const int ke = pmy_mb->ke;
  Real xCI;
  Real NH_ghostzone;
  Real NH_boundary, NH2_boundary, NCO_boundary, NC_boundary;
  if (direction == BoundaryFace::inner_x1) {
    // +x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        NH_ghostzone = lunit *
          pmy_mb->phydro->w(IDN, k, j, is-1) * pmy_mb->pcoord->dx1f(is-1) * f_prev;
        NH_boundary = col(direction, k, j, is-1, pmy_chemnet->iNHtot_) + NH_ghostzone;
        NH2_boundary = col(direction, k, j, is-1, pmy_chemnet->iNH2_)
          + pmy_mb->pscalars->r(iH2, k, j, is-1) * NH_ghostzone;
        NCO_boundary = col(direction, k, j, is-1, pmy_chemnet->iNCO_)
          + pmy_mb->pscalars->r(iCO, k, j, is-1) * NH_ghostzone;
        xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, is-1)
          - pmy_mb->pscalars->r(iCplus, k, j, is-1)
          - pmy_mb->pscalars->r(iHCOplus, k, j, is-1)
          - pmy_mb->pscalars->r(iCHx, k, j, is-1);
        if (xCI < 0.0) {
          xCI = 0.;
        }
        NC_boundary = col(direction, k, j, is-1, pmy_chemnet->iNC_)
          + xCI * NH_ghostzone;
        for (int i=is; i<=ie; ++i) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
          col(direction, k, j, i, pmy_chemnet->iNC_) += NC_boundary;
        }
      }
    }
  } else if (direction == BoundaryFace::outer_x1) {
    // -x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        NH_ghostzone = lunit *
          pmy_mb->phydro->w(IDN, k, j, ie+1) * pmy_mb->pcoord->dx1f(ie+1) * f_prev;
        NH_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNHtot_) + NH_ghostzone;
        NH2_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNH2_)
          + pmy_mb->pscalars->r(iH2, k, j, ie+1) * NH_ghostzone;
        NCO_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNCO_)
          + pmy_mb->pscalars->r(iCO, k, j, ie+1) * NH_ghostzone;
        xCI = xCtot - pmy_mb->pscalars->r(iCO, k, j, ie+1)
          - pmy_mb->pscalars->r(iCplus, k, j, ie+1)
          - pmy_mb->pscalars->r(iHCOplus, k, j, ie+1)
          - pmy_mb->pscalars->r(iCHx, k, j, ie+1);
        if (xCI < 0.0) {
          xCI = 0.;
        }
        NC_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNC_)
          + xCI * NH_ghostzone;
        for (int i=ie; i>=is; --i) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
          col(direction, k, j, i, pmy_chemnet->iNC_) += NC_boundary;
        }
      }
    }
  } else if (direction == BoundaryFace::inner_x2) {
    // +y
    if (js != 0) { // y dimension included
      for (int k=ks; k<=ke; ++k) {
        for (int i=is; i<=ie; ++i) {
          NH_ghostzone = lunit *
            pmy_mb->phydro->w(IDN, k, js-1, i) * pmy_mb->pcoord->dx2f(js-1) * f_prev;
          NH_boundary = col(direction, k, js-1, i, pmy_chemnet->iNHtot_) + NH_ghostzone;
          NH2_boundary = col(direction, k, js-1, i, pmy_chemnet->iNH2_)
            + pmy_mb->pscalars->r(iH2, k, js-1, i) * NH_ghostzone;
          NCO_boundary = col(direction, k, js-1, i, pmy_chemnet->iNCO_)
            + pmy_mb->pscalars->r(iCO, k, js-1, i) * NH_ghostzone;
          xCI = xCtot - pmy_mb->pscalars->r(iCO, k, js-1, i)
            - pmy_mb->pscalars->r(iCplus, k, js-1, i)
            - pmy_mb->pscalars->r(iHCOplus, k, js-1, i)
            - pmy_mb->pscalars->r(iCHx, k, js-1, i);
          if (xCI < 0.0) {
            xCI = 0.;
          }
          NC_boundary = col(direction, k, js-1, i, pmy_chemnet->iNC_)
            + xCI * NH_ghostzone;
          for (int j=js; j<=je; ++j) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
            col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
            col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
            col(direction, k, j, i, pmy_chemnet->iNC_) += NC_boundary;
          }
        }
      }
    }
  } else if (direction == BoundaryFace::outer_x2) {
    // -y
    if (je != 0) { // y dimension included
      for (int k=ks; k<=ke; ++k) {
        for (int i=is; i<=ie; ++i) {
          NH_ghostzone = lunit *
            pmy_mb->phydro->w(IDN, k, je+1, i) * pmy_mb->pcoord->dx2f(je+1) * f_prev;
          NH_boundary = col(direction, k, je+1, i, pmy_chemnet->iNHtot_) + NH_ghostzone;
          NH2_boundary = col(direction, k, je+1, i, pmy_chemnet->iNH2_)
            + pmy_mb->pscalars->r(iH2, k, je+1, i) * NH_ghostzone;
          NCO_boundary = col(direction, k, je+1, i, pmy_chemnet->iNCO_)
            + pmy_mb->pscalars->r(iCO, k, je+1, i) * NH_ghostzone;
          xCI = xCtot - pmy_mb->pscalars->r(iCO, k, je+1, i)
            - pmy_mb->pscalars->r(iCplus, k, je+1, i)
            - pmy_mb->pscalars->r(iHCOplus, k, je+1, i)
            - pmy_mb->pscalars->r(iCHx, k, je+1, i);
          if (xCI < 0.0) {
            xCI = 0.;
          }
          NC_boundary = col(direction, k, je+1, i, pmy_chemnet->iNC_)
            + xCI * NH_ghostzone;
          for (int j=je; j>=js; --j) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
            col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
            col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
            col(direction, k, j, i, pmy_chemnet->iNC_) += NC_boundary;
          }
        }
      }
    }
  } else if (direction == BoundaryFace::inner_x3) {
    // +z
    if (ks != 0) { // z dimension included
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NH_ghostzone = lunit *
            pmy_mb->phydro->w(IDN, ks-1, j, i) * pmy_mb->pcoord->dx3f(ks-1) * f_prev;
          NH_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNHtot_) + NH_ghostzone;
          NH2_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNH2_)
            + pmy_mb->pscalars->r(iH2, ks-1, j, i) * NH_ghostzone;
          NCO_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNCO_)
            + pmy_mb->pscalars->r(iCO, ks-1, j, i) * NH_ghostzone;
          xCI = xCtot - pmy_mb->pscalars->r(iCO, ks-1, j, i)
            - pmy_mb->pscalars->r(iCplus, ks-1, j, i)
            - pmy_mb->pscalars->r(iHCOplus, ks-1, j, i)
            - pmy_mb->pscalars->r(iCHx, ks-1, j, i);
          if (xCI < 0.0) {
            xCI = 0.;
          }
          NC_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNC_)
            + xCI * NH_ghostzone;
          for (int k=ks; k<=ke; ++k) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
            col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
            col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
            col(direction, k, j, i, pmy_chemnet->iNC_) += NC_boundary;
          }
        }
      }
    }
  } else if (direction == BoundaryFace::outer_x3) {
    // -z
    if (ke != 0) { // z dimension included
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NH_ghostzone = lunit *
            pmy_mb->phydro->w(IDN, ke+1, j, i) * pmy_mb->pcoord->dx3f(ke+1) * f_prev;
          NH_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNHtot_) + NH_ghostzone;
          NH2_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNH2_)
            + pmy_mb->pscalars->r(iH2, ke+1, j, i) * NH_ghostzone;
          NCO_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNCO_)
            + pmy_mb->pscalars->r(iCO, ke+1, j, i) * NH_ghostzone;
          xCI = xCtot - pmy_mb->pscalars->r(iCO, ke+1, j, i)
            - pmy_mb->pscalars->r(iCplus, ke+1, j, i)
            - pmy_mb->pscalars->r(iHCOplus, ke+1, j, i)
            - pmy_mb->pscalars->r(iCHx, ke+1, j, i);
          if (xCI < 0.0) {
            xCI = 0.;
          }
          NC_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNC_)
            + xCI * NH_ghostzone;
          for (int k=ke; k>=ks; --k) {
            col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
            col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
            col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
            col(direction, k, j, i, pmy_chemnet->iNC_) += NC_boundary;
          }
        }
      }
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ChemRadIntegrator six_ray [UpdateCol]" << std::endl
      << "direction:" << direction << " unknown." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}
