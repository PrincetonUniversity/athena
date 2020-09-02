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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationDonorCellX1(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie+1; ++i) {
        rad_l_(lm,i) = intensity(lm,k,j,i-1);
        rad_r_(lm,i) = intensity(lm,k,j,i);
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationDonorCellX2(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie; ++i) {
        rad_l_(lm,i) = intensity(lm,k,j-1,i);
        rad_r_(lm,i) = intensity(lm,k,j,i);
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationDonorCellX3(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie; ++i) {
        rad_l_(lm,i) = intensity(lm,k-1,j,i);
        rad_r_(lm,i) = intensity(lm,k,j,i);
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationDonorCellA1(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze+1; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm_l = AngleInd(l - 1, m, false, false);
      int lm_c = AngleInd(l, m, true, false);
      int lm_r = AngleInd(l, m, false, false);
      for (int i = is; i <= ie; ++i) {
        rad_l_(lm_c,i) = intensity(lm_l,k,j,i);
        rad_r_(lm_c,i) = intensity(lm_r,k,j,i);
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationDonorCellA2(const AthenaArray<Real> &intensity, int k, int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      int lm_l = AngleInd(l, m - 1, false, false);
      int lm_c = AngleInd(l, m, false, true);
      int lm_r = AngleInd(l, m, false, false);
      for (int i = is; i <= ie; ++i) {
        rad_l_(lm_c,i) = intensity(lm_l,k,j,i);
        rad_r_(lm_c,i) = intensity(lm_r,k,j,i);
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearX1(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int i = is; i <= ie+1; ++i) {
        Real x_l = pmy_block->pcoord->x1v(i-1);
        Real x_c = pmy_block->pcoord->x1f(i);
        Real x_r = pmy_block->pcoord->x1v(i);
        Real dx_l = x_c - x_l;
        Real dx_r = x_r - x_c;
        Real q_ll = intensity(lm,k,j,i-2);
        Real q_l = intensity(lm,k,j,i-1);
        Real q_r = intensity(lm,k,j,i);
        Real q_rr = intensity(lm,k,j,i+1);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq_2_l = dq_l * dq_c;
        Real dq_2_r = dq_c * dq_r;
        Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
        Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
        rad_l_(lm,i) = q_l + dx_l * dq_m_l;
        rad_r_(lm,i) = q_r - dx_r * dq_m_r;
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearX2(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      Real x_l = pmy_block->pcoord->x2v(j-1);
      Real x_c = pmy_block->pcoord->x2f(j);
      Real x_r = pmy_block->pcoord->x2v(j);
      Real dx_l = x_c - x_l;
      Real dx_r = x_r - x_c;
      for (int i = is; i <= ie; ++i) {
        Real q_ll = intensity(lm,k,j-2,i);
        Real q_l = intensity(lm,k,j-1,i);
        Real q_r = intensity(lm,k,j,i);
        Real q_rr = intensity(lm,k,j+1,i);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq_2_l = dq_l * dq_c;
        Real dq_2_r = dq_c * dq_r;
        Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
        Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
        rad_l_(lm,i) = q_l + dx_l * dq_m_l;
        rad_r_(lm,i) = q_r - dx_r * dq_m_r;
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearX3(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      Real x_l = pmy_block->pcoord->x3v(k-1);
      Real x_c = pmy_block->pcoord->x3f(k);
      Real x_r = pmy_block->pcoord->x3v(k);
      Real dx_l = x_c - x_l;
      Real dx_r = x_r - x_c;
      for (int i = is; i <= ie; ++i) {
        Real q_ll = intensity(lm,k-2,j,i);
        Real q_l = intensity(lm,k-1,j,i);
        Real q_r = intensity(lm,k,j,i);
        Real q_rr = intensity(lm,k+1,j,i);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq_2_l = dq_l * dq_c;
        Real dq_2_r = dq_c * dq_r;
        Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
        Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
        rad_l_(lm,i) = q_l + dx_l * dq_m_l;
        rad_r_(lm,i) = q_r - dx_r * dq_m_r;
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearA1(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze+1; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm_ll = AngleInd(l - 2, m, false, false);
      int lm_l = AngleInd(l - 1, m, false, false);
      int lm_c = AngleInd(l, m, true, false);
      int lm_r = AngleInd(l, m, false, false);
      int lm_rr = AngleInd(l + 1, m, false, false);
      Real x_l = zetav(l-1);
      Real x_c = zetaf(l);
      Real x_r = zetav(l);
      Real dx_l = x_c - x_l;
      Real dx_r = x_r - x_c;
      for (int i = is; i <= ie; ++i) {
        Real q_ll = intensity(lm_ll,k,j,i);
        Real q_l = intensity(lm_l,k,j,i);
        Real q_r = intensity(lm_r,k,j,i);
        Real q_rr = intensity(lm_rr,k,j,i);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq_2_l = dq_l * dq_c;
        Real dq_2_r = dq_c * dq_r;
        Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
        Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
        rad_l_(lm_c,i) = q_l + dx_l * dq_m_l;
        rad_r_(lm_c,i) = q_r - dx_r * dq_m_r;
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
//   this->rad_l_, this->rad_r_: reconstructed intensities set

void Radiation::RadiationPiecewiseLinearA2(const AthenaArray<Real> &intensity, int k,
    int j) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      int lm_ll = AngleInd(l, m - 2, false, false);
      int lm_l = AngleInd(l, m - 1, false, false);
      int lm_c = AngleInd(l, m, false, true);
      int lm_r = AngleInd(l, m, false, false);
      int lm_rr = AngleInd(l, m + 1, false, false);
      Real x_l = psiv(m-1);
      Real x_c = psif(m);
      Real x_r = psiv(m);
      Real dx_l = x_c - x_l;
      Real dx_r = x_r - x_c;
      for (int i = is; i <= ie; ++i) {
        Real q_ll = intensity(lm_ll,k,j,i);
        Real q_l = intensity(lm_l,k,j,i);
        Real q_r = intensity(lm_r,k,j,i);
        Real q_rr = intensity(lm_rr,k,j,i);
        Real dq_l = q_l - q_ll;
        Real dq_c = q_r - q_l;
        Real dq_r = q_rr - q_r;
        Real dq_2_l = dq_l * dq_c;
        Real dq_2_r = dq_c * dq_r;
        Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
        Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
        rad_l_(lm_c,i) = q_l + dx_l * dq_m_l;
        rad_r_(lm_c,i) = q_r - dx_r * dq_m_r;
      }
    }
  }
  return;
}
