//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of reconstruction-related functions in class Radiation

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates

//----------------------------------------------------------------------------------------
// Donor cell (1st-order) reconstruction in x1-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationDonorCellX1(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie+1; ++i) {
        ii_l_(lm,i) = intensity(lm,k,j,i-1);
        ii_r_(lm,i) = intensity(lm,k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Donor cell (1st-order) reconstruction in x2-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationDonorCellX2(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie; ++i) {
        ii_l_(lm,i) = intensity(lm,k,j-1,i);
        ii_r_(lm,i) = intensity(lm,k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Donor cell (1st-order) reconstruction in x3-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationDonorCellX3(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie; ++i) {
        ii_l_(lm,i) = intensity(lm,k-1,j,i);
        ii_r_(lm,i) = intensity(lm,k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Donor cell (1st-order) reconstruction in zeta-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationDonorCellA1(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze+1; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm_l = AngleInd(l - 1, m, false, false);
      int lm_c = AngleInd(l, m, true, false);
      int lm_r = AngleInd(l, m, false, false);
      for (int i = is; i <= ie; ++i) {
        if (na1_n_0_(l,m,k,j,i) != 0.0) {
          Real n_0_face = na1_n_0_(l,m,k,j,i) / na1_(l,m,k,j,i);
          Real n_0_cell = n0_n_mu_(0,l-1,m,k,j,i) / nmu_(0,l-1,m,k,j,i);
          ii_l_(lm_c,i) = intensity(lm_l,k,j,i) * n_0_cell / n_0_face;
          n_0_cell = n0_n_mu_(0,l,m,k,j,i) / nmu_(0,l,m,k,j,i);
          ii_r_(lm_c,i) = intensity(lm_r,k,j,i) * n_0_cell / n_0_face;
        } else {
          ii_l_(lm_c,i) = 0.0;
          ii_r_(lm_c,i) = 0.0;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Donor cell (1st-order) reconstruction in psi-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationDonorCellA2(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      int lm_l = AngleInd(l, m - 1, false, false);
      int lm_c = AngleInd(l, m, false, true);
      int lm_r = AngleInd(l, m, false, false);
      for (int i = is; i <= ie; ++i) {
        if (na2_n_0_(l,m,k,j,i) != 0.0) {
          Real n_0_face = na2_n_0_(l,m,k,j,i) / na2_(l,m,k,j,i);
          Real n_0_cell = n0_n_mu_(0,l,m-1,k,j,i) / nmu_(0,l,m-1,k,j,i);
          ii_l_(lm_c,i) = intensity(lm_l,k,j,i) * n_0_cell / n_0_face;
          n_0_cell = n0_n_mu_(0,l,m,k,j,i) / nmu_(0,l,m,k,j,i);
          ii_r_(lm_c,i) = intensity(lm_r,k,j,i) * n_0_cell / n_0_face;
        } else {
          ii_l_(lm_c,i) = 0.0;
          ii_r_(lm_c,i) = 0.0;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Piecewise linear (2nd-order) reconstruction in x1-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearX1(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie+1; ++i) {
        Real x_lll = pmy_block->pcoord->x1v(i-2);
        Real x_ll = pmy_block->pcoord->x1f(i-1);
        Real x_l = pmy_block->pcoord->x1v(i-1);
        Real x_c = pmy_block->pcoord->x1f(i);
        Real x_r = pmy_block->pcoord->x1v(i);
        Real x_rr = pmy_block->pcoord->x1f(i+1);
        Real x_rrr = pmy_block->pcoord->x1v(i+1);
        Real c_l = (x_c - x_l) / (x_c - x_ll);
        Real c_r = (x_r - x_c) / (x_rr - x_c);
        Real cb_l = (x_l - x_lll) / (x_l - x_ll);
        Real cf_l = (x_r - x_l) / (x_c - x_l);
        Real cb_r = (x_r - x_l) / (x_r - x_c);
        Real cf_r = (x_rrr - x_r) / (x_rr - x_r);
        Real q_ll = intensity(lm,k,j,i-2);
        Real q_l = intensity(lm,k,j,i-1);
        Real q_r = intensity(lm,k,j,i);
        Real q_rr = intensity(lm,k,j,i+1);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq2_l = dq_l * dq_c;
        Real dq2_r = dq_c * dq_r;
        Real dqm_l = (dq2_l > 0.0) ? dq2_l * (cf_l * dq_l + cb_l * dq_c)
            / (SQR(dq_l) + SQR(dq_c) + (cf_l + cb_l - 2.0) * dq2_l) : 0.0;
        Real dqm_r = (dq2_r > 0.0) ? dq2_r * (cf_r * dq_c + cb_r * dq_r)
            / (SQR(dq_c) + SQR(dq_r) + (cf_r + cb_r - 2.0) * dq2_r) : 0.0;
        ii_l_(lm,i) = q_l + c_l * dqm_l;
        ii_r_(lm,i) = q_r - c_r * dqm_r;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Piecewise linear (2nd-order) reconstruction in x2-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearX2(const AthenaArray<Real> &intensity, int k,
    int j) {
  Real x_lll = pmy_block->pcoord->x2v(j-2);
  Real x_ll = pmy_block->pcoord->x2f(j-1);
  Real x_l = pmy_block->pcoord->x2v(j-1);
  Real x_c = pmy_block->pcoord->x2f(j);
  Real x_r = pmy_block->pcoord->x2v(j);
  Real x_rr = pmy_block->pcoord->x2f(j+1);
  Real x_rrr = pmy_block->pcoord->x2v(j+1);
  Real c_l = (x_c - x_l) / (x_c - x_ll);
  Real c_r = (x_r - x_c) / (x_rr - x_c);
  Real cb_l = (x_l - x_lll) / (x_l - x_ll);
  Real cf_l = (x_r - x_l) / (x_c - x_l);
  Real cb_r = (x_r - x_l) / (x_r - x_c);
  Real cf_r = (x_rrr - x_r) / (x_rr - x_r);
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie; ++i) {
        Real q_ll = intensity(lm,k,j-2,i);
        Real q_l = intensity(lm,k,j-1,i);
        Real q_r = intensity(lm,k,j,i);
        Real q_rr = intensity(lm,k,j+1,i);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq2_l = dq_l * dq_c;
        Real dq2_r = dq_c * dq_r;
        Real dqm_l = (dq2_l > 0.0) ? dq2_l * (cf_l * dq_l + cb_l * dq_c)
            / (SQR(dq_l) + SQR(dq_c) + (cf_l + cb_l - 2.0) * dq2_l) : 0.0;
        Real dqm_r = (dq2_r > 0.0) ? dq2_r * (cf_r * dq_c + cb_r * dq_r)
            / (SQR(dq_c) + SQR(dq_r) + (cf_r + cb_r - 2.0) * dq2_r) : 0.0;
        ii_l_(lm,i) = q_l + c_l * dqm_l;
        ii_r_(lm,i) = q_r - c_r * dqm_r;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Piecewise linear (2nd-order) reconstruction in x3-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearX3(const AthenaArray<Real> &intensity, int k,
    int j) {
  Real x_lll = pmy_block->pcoord->x3v(k-2);
  Real x_ll = pmy_block->pcoord->x3f(k-1);
  Real x_l = pmy_block->pcoord->x3v(k-1);
  Real x_c = pmy_block->pcoord->x3f(k);
  Real x_r = pmy_block->pcoord->x3v(k);
  Real x_rr = pmy_block->pcoord->x3f(k+1);
  Real x_rrr = pmy_block->pcoord->x3v(k+1);
  Real c_l = (x_c - x_l) / (x_c - x_ll);
  Real c_r = (x_r - x_c) / (x_rr - x_c);
  Real cb_l = (x_l - x_lll) / (x_l - x_ll);
  Real cf_l = (x_r - x_l) / (x_c - x_l);
  Real cb_r = (x_r - x_l) / (x_r - x_c);
  Real cf_r = (x_rrr - x_r) / (x_rr - x_r);
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie; ++i) {
        Real q_ll = intensity(lm,k-2,j,i);
        Real q_l = intensity(lm,k-1,j,i);
        Real q_r = intensity(lm,k,j,i);
        Real q_rr = intensity(lm,k+1,j,i);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq2_l = dq_l * dq_c;
        Real dq2_r = dq_c * dq_r;
        Real dqm_l = (dq2_l > 0.0) ? dq2_l * (cf_l * dq_l + cb_l * dq_c)
            / (SQR(dq_l) + SQR(dq_c) + (cf_l + cb_l - 2.0) * dq2_l) : 0.0;
        Real dqm_r = (dq2_r > 0.0) ? dq2_r * (cf_r * dq_c + cb_r * dq_r)
            / (SQR(dq_c) + SQR(dq_r) + (cf_r + cb_r - 2.0) * dq2_r) : 0.0;
        ii_l_(lm,i) = q_l + c_l * dqm_l;
        ii_r_(lm,i) = q_r - c_r * dqm_r;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Piecewise linear (2nd-order) reconstruction in zeta-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearA1(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze+1; ++l) {
    Real x_lll = zetav(l-2);
    Real x_ll = zetaf(l-1);
    Real x_l = zetav(l-1);
    Real x_c = zetaf(l);
    Real x_r = zetav(l);
    Real x_rr = zetaf(l+1);
    Real x_rrr = zetav(l+1);
    Real c_l = (x_c - x_l) / (x_c - x_ll);
    Real c_r = (x_r - x_c) / (x_rr - x_c);
    Real cb_l = (x_l - x_lll) / (x_l - x_ll);
    Real cf_l = (x_r - x_l) / (x_c - x_l);
    Real cb_r = (x_r - x_l) / (x_r - x_c);
    Real cf_r = (x_rrr - x_r) / (x_rr - x_r);
    for (int m = ps; m <= pe; ++m) {
      int lm_ll = AngleInd(l - 2, m, false, false);
      int lm_l = AngleInd(l - 1, m, false, false);
      int lm_c = AngleInd(l, m, true, false);
      int lm_r = AngleInd(l, m, false, false);
      int lm_rr = AngleInd(l + 1, m, false, false);
      for (int i = is; i <= ie; ++i) {
        Real n_0_cell = n0_n_mu_(0,l-2,m,k,j,i) / nmu_(0,l-2,m,k,j,i);
        Real q_ll = intensity(lm_ll,k,j,i) * n_0_cell;
        n_0_cell = n0_n_mu_(0,l-1,m,k,j,i) / nmu_(0,l-1,m,k,j,i);
        Real q_l = intensity(lm_l,k,j,i) * n_0_cell;
        n_0_cell = n0_n_mu_(0,l,m,k,j,i) / nmu_(0,l,m,k,j,i);
        Real q_r = intensity(lm_r,k,j,i) * n_0_cell;
        n_0_cell = n0_n_mu_(0,l+1,m,k,j,i) / nmu_(0,l+1,m,k,j,i);
        Real q_rr = intensity(lm_rr,k,j,i) * n_0_cell;
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq2_l = dq_l * dq_c;
        Real dq2_r = dq_c * dq_r;
        Real dqm_l = (dq2_l > 0.0) ? dq2_l * (cf_l * dq_l + cb_l * dq_c)
            / (SQR(dq_l) + SQR(dq_c) + (cf_l + cb_l - 2.0) * dq2_l) : 0.0;
        Real dqm_r = (dq2_r > 0.0) ? dq2_r * (cf_r * dq_c + cb_r * dq_r)
            / (SQR(dq_c) + SQR(dq_r) + (cf_r + cb_r - 2.0) * dq2_r) : 0.0;
        if (na1_n_0_(l,m,k,j,i) != 0.0) {
          Real n_0_face = na1_n_0_(l,m,k,j,i) / na1_(l,m,k,j,i);
          ii_l_(lm_c,i) = (q_l + c_l * dqm_l) / n_0_face;
          ii_r_(lm_c,i) = (q_r - c_r * dqm_r) / n_0_face;
        } else {
          ii_l_(lm_c,i) = 0.0;
          ii_r_(lm_c,i) = 0.0;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Piecewise linear (2nd-order) reconstruction in psi-direction
// Inputs:
//   intensity: radiation intensity
//   k, j: x3- and x2-indices
// Outputs:
//   this->ii_l_, this->ii_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearA2(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      int lm_ll = AngleInd(l, m - 2, false, false);
      int lm_l = AngleInd(l, m - 1, false, false);
      int lm_c = AngleInd(l, m, false, true);
      int lm_r = AngleInd(l, m, false, false);
      int lm_rr = AngleInd(l, m + 1, false, false);
      Real x_lll = psiv(m-2);
      Real x_ll = psif(m-1);
      Real x_l = psiv(m-1);
      Real x_c = psif(m);
      Real x_r = psiv(m);
      Real x_rr = psif(m+1);
      Real x_rrr = psiv(m+1);
      Real c_l = (x_c - x_l) / (x_c - x_ll);
      Real c_r = (x_r - x_c) / (x_rr - x_c);
      Real cb_l = (x_l - x_lll) / (x_l - x_ll);
      Real cf_l = (x_r - x_l) / (x_c - x_l);
      Real cb_r = (x_r - x_l) / (x_r - x_c);
      Real cf_r = (x_rrr - x_r) / (x_rr - x_r);
      for (int i = is; i <= ie; ++i) {
        Real n_0_cell = n0_n_mu_(0,l,m-2,k,j,i) / nmu_(0,l,m-2,k,j,i);
        Real q_ll = intensity(lm_ll,k,j,i) * n_0_cell;
        n_0_cell = n0_n_mu_(0,l,m-1,k,j,i) / nmu_(0,l,m-1,k,j,i);
        Real q_l = intensity(lm_l,k,j,i) * n_0_cell;
        n_0_cell = n0_n_mu_(0,l,m,k,j,i) / nmu_(0,l,m,k,j,i);
        Real q_r = intensity(lm_r,k,j,i) * n_0_cell;
        n_0_cell = n0_n_mu_(0,l,m+1,k,j,i) / nmu_(0,l,m+1,k,j,i);
        Real q_rr = intensity(lm_rr,k,j,i) * n_0_cell;
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq2_l = dq_l * dq_c;
        Real dq2_r = dq_c * dq_r;
        Real dqm_l = (dq2_l > 0.0) ? dq2_l * (cf_l * dq_l + cb_l * dq_c)
            / (SQR(dq_l) + SQR(dq_c) + (cf_l + cb_l - 2.0) * dq2_l) : 0.0;
        Real dqm_r = (dq2_r > 0.0) ? dq2_r * (cf_r * dq_c + cb_r * dq_r)
            / (SQR(dq_c) + SQR(dq_r) + (cf_r + cb_r - 2.0) * dq2_r) : 0.0;
        if (na2_n_0_(l,m,k,j,i) != 0.0) {
          Real n_0_face = na2_n_0_(l,m,k,j,i) / na2_(l,m,k,j,i);
          ii_l_(lm_c,i) = (q_l + c_l * dqm_l) / n_0_face;
          ii_r_(lm_c,i) = (q_r - c_r * dqm_r) / n_0_face;
        } else {
          ii_l_(lm_c,i) = 0.0;
          ii_r_(lm_c,i) = 0.0;
        }
      }
    }
  }
  return;
}
