//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file six_ray.cpp
//! \brief implementation of radiation integrators: six_ray

// this class header
#include "rad_integrators.hpp"

//c++ headers
#include <sstream>    // stringstream
#include <iostream>   // endl

// Athena++ headers
#include "../../coordinates/coordinates.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../radiation.hpp"
#include "../../scalars/scalars.hpp"
#include "../../bvals/bvals.hpp"
#include "../../utils/units.hpp"

namespace {
  Real rad_G0; //diffuse radiation field strength in Draine 1987 unit
  Real f_cell, f_prev; //fraction of the column in the cell that goes to shielding
#ifdef INCLUDE_CHEMISTRY
  Real lunit; //length unit in cm
#endif //INCLUDE_CHEMISTRY
}

//----------------------------------------------------------------------------------------
//! constructor, for six_ray radiation integrator
RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin) {
  pmy_mb = prad->pmy_block;
  pmy_rad = prad;
  rad_G0 = pin->GetReal("radiation", "G0");
  f_cell = pin->GetOrAddReal("radiation", "shielding_fraction_cell", 0.5);
  std::stringstream msg; //error message
  if (pmy_rad->nang != 6) {
    msg << "### FATAL ERROR in RadIntegrator constructor [RadIntegrator]" << std::endl
      << "Six-ray scheme with nang != 6." << std::endl;
    ATHENA_ERROR(msg);
  }
#ifdef INCLUDE_CHEMISTRY
  pmy_chemnet = &pmy_mb->pscalars->chemnet;
  lunit = pmy_chemnet->punit->Length;
  ncol = pmy_chemnet->n_cols_;
  //allocate array for column density
  int ncells1 = pmy_mb->block_size.nx1 + 2; //one ghost cell on each side
  int ncells2 = 1, ncells3 = 1;
  if (pmy_mb->block_size.nx2 > 1) ncells2 = pmy_mb->block_size.nx2 + 2;
  if (pmy_mb->block_size.nx3 > 1) ncells3 = pmy_mb->block_size.nx3 + 2;
  col.NewAthenaArray(6, ncells3, ncells2, ncells1, ncol);
#ifdef DEBUG
  col_avg.NewAthenaArray(ncol, ncells3, ncells2, ncells1);
  col_Htot.NewAthenaArray(6, ncells3, ncells2, ncells1);
  col_H2.NewAthenaArray(6, ncells3, ncells2, ncells1);
  col_CO.NewAthenaArray(6, ncells3, ncells2, ncells1);
  col_C.NewAthenaArray(6, ncells3, ncells2, ncells1);
#endif //DEBUG
#endif //INCLUDE_CHEMISTRY
}

//----------------------------------------------------------------------------------------
//! destructor
RadIntegrator::~RadIntegrator() {
#ifdef INCLUDE_CHEMISTRY
  col.DeleteAthenaArray();
#ifdef DEBUG
  col_avg.DeleteAthenaArray();
  col_Htot.DeleteAthenaArray();
  col_H2.DeleteAthenaArray();
  col_CO.DeleteAthenaArray();
  col_C.DeleteAthenaArray();
#endif //DEBUG
#endif //INCLUDE_CHEMISTRY
}

//----------------------------------------------------------------------------------------
//! \fn void RadIntegrator::CopyToOutput()
//! \brief average radiation field over all angles and copy values to output
void RadIntegrator::CopyToOutput() {
  int is = pmy_mb->is;
  int js = pmy_mb->js;
  int ks = pmy_mb->ks;
  int ie = pmy_mb->ie;
  int je = pmy_mb->je;
  int ke = pmy_mb->ke;
  int iang_arr[6] =
  {BoundaryFace::inner_x1, BoundaryFace::inner_x2, BoundaryFace::inner_x3,
   BoundaryFace::outer_x1, BoundaryFace::outer_x2, BoundaryFace::outer_x3};
  for (int k=ks-NGHOST; k<=ke+NGHOST; ++k) {
    for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        for (int iang=0; iang < 6; iang++) {
          //column densities
#ifdef DEBUG
          col_Htot(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNHtot_);
          col_H2(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNH2_);
          col_CO(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNCO_);
          col_C(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNC_);
          //angel averaged column densities
          for (int icol=0; icol<ncol; icol++) {
            if (iang == 0) {
              col_avg(icol, k, j, i) = 0;
            }
            col_avg(icol, k, j, i) += col(iang, k, j, i, icol)/6.;
          }
#endif //DEBUG
          //radiation
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
//! \fn void RadIntegrator::UpdateRadiation(int direction)
//! \brief calcuate total column and update radiation
void RadIntegrator::UpdateRadiation(int direction) {}

#ifdef INCLUDE_CHEMISTRY
//----------------------------------------------------------------------------------------
//! \fn void void GetColMB(int direction)
//! \brief calculate column densities within the meshblock
//!
//! direction: 0:+x, 1:+y, 2:+z, 3:-x, 4:-y, 5:-z
void RadIntegrator::GetColMB(int direction) {
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
  std::stringstream msg; //error message
  if (direction == BoundaryFace::inner_x1) {
    //+x
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
            if (xCI < 0) {
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
            if (xCI < 0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j, i-1)
                  - pmy_mb->pscalars->r(iCplus, k, j, i-1)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i-1)
                  - pmy_mb->pscalars->r(iCHx, k, j, i-1);
            if (xCI_prev < 0) {
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
    //-x
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
            if (xCI < 0) {
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
            if (xCI < 0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j, i+1)
                  - pmy_mb->pscalars->r(iCplus, k, j, i+1)
                  - pmy_mb->pscalars->r(iHCOplus, k, j, i+1)
                  - pmy_mb->pscalars->r(iCHx, k, j, i+1);
            if (xCI_prev < 0) {
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
    //+y
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
            if (xCI < 0) {
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
            if (xCI < 0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j-1, i)
                  - pmy_mb->pscalars->r(iCplus, k, j-1, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j-1, i)
                  - pmy_mb->pscalars->r(iCHx, k, j-1, i);
            if (xCI_prev < 0) {
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
    //-y
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
            if (xCI < 0) {
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
            if (xCI < 0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k, j+1, i)
                  - pmy_mb->pscalars->r(iCplus, k, j+1, i)
                  - pmy_mb->pscalars->r(iHCOplus, k, j+1, i)
                  - pmy_mb->pscalars->r(iCHx, k, j+1, i);
            if (xCI_prev < 0) {
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
    //+z
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
            if (xCI < 0) {
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
            if (xCI < 0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k-1, j, i)
                  - pmy_mb->pscalars->r(iCplus, k-1, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k-1, j, i)
                  - pmy_mb->pscalars->r(iCHx, k-1, j, i);
            if (xCI_prev < 0) {
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
    //-z
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
            if (xCI < 0) {
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
            if (xCI < 0) {
              xCI = 0.;
            }
            xCI_prev = xCtot - pmy_mb->pscalars->r(iCO, k+1, j, i)
                  - pmy_mb->pscalars->r(iCplus, k+1, j, i)
                  - pmy_mb->pscalars->r(iHCOplus, k+1, j, i)
                  - pmy_mb->pscalars->r(iCHx, k+1, j, i);
            if (xCI_prev < 0) {
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
    msg << "### FATAL ERROR in RadIntegrator six_ray [GetColMB]" << std::endl
      << "direction {0,1,2,3,4,5}:" << direction << " unknown." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}
#endif //INCLUDE_CHEMISTRY
