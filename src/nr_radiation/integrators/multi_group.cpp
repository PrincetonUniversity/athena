//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file multi_group.cpp
//  \brief  Add multi_group  source terms
//======================================================================================

// C headers

// C++ headers
#include <sstream>  // msg
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
#include "../radiation.hpp"

// this class header
#include "./rad_integrators.hpp"


//--------------------------------------------------------------------------------------
//The function calculate the amount of shift for each I we need
// The default frequency grid is nu_grid with nfreq frequency bins
// This covers from 0 to infty
// the original frequency grid is nu_grid(0,...,nfreq-2)
// the shifted frequency grid is cm_nu(0,...,nfreq-2)
// for each angle, the shift is the same \gamma (1-nv/c)
// The procedure is: Interpolate the intensity to the frequency grid nu_grid from cm_nu
// solve the source terms and update co-moving frame intensity
// interpolate the intensity back to the frequency grid cm_nu from nu_grid for each angle
// Then transform the intensity back to the lab-frame

//Giving the frequency integrated intensity for each group in the co-moving frame,
// get the monochromatic intensity at the face and center of the frequency grid
//monochromatic intensity is always 0 at nu=0, the values at frequency grid faces should
// satisfy 0.5*(ir_face(f)+ir_face(f+1))*delta_nu = ir_f
// as we assume piecewise linear spectrum shape, and continuous intensity across frequency
void RadIntegrator::MapLabToCmFrequency(Real &tran_coef,
                   AthenaArray<Real> &split_ratio,
                   AthenaArray<int> &map_start, AthenaArray<int> &map_end,
                   AthenaArray<Real> &ir_cm, AthenaArray<Real> &ir_shift) {
  const int& nfreq=pmy_rad->nfreq;

  // prepare the frequency bin width
  for (int ifr=0; ifr<nfreq-1; ++ifr) {
    delta_nu_n_(ifr) = pmy_rad->delta_nu(ifr) * tran_coef;
  }

  GetCmMCIntensity(ir_cm, delta_nu_n_, ir_face_);
  // calculate the shift ratio
  ForwardSplitting(tran_coef, ir_cm, ir_face_, split_ratio, map_start, map_end);
  MapIrcmFrequency(split_ratio, map_start, map_end, ir_cm,ir_shift);
  return;
}

void RadIntegrator::GetCmMCIntensity(AthenaArray<Real> &ir_cm,
                          AthenaArray<Real> &delta_nu_n, AthenaArray<Real> &ir_face) {
  const int& nfreq = pmy_rad->nfreq;
  if (nfreq > 1) {
    ir_face(0) = 0.0;
    for (int ifr=1; ifr<nfreq; ++ifr) {
      ir_face(ifr) = std::max(2.0*ir_cm(ifr-1)/delta_nu_n(ifr-1) -
                                             ir_face(ifr-1), 0.0);
    }
  }
}

// general function to split any array in the frequency bin [\Gamma nu_f]
// to the frequency bin [nu_f]
// ir_last_bin is used to determine the shift in the last frequency bin
// assuming BlackBody spectrum in the last frequency bin
// In other bins, we assume piecewise constant, split each bin according
// to frequency overlap
void RadIntegrator::ForwardSplitting(
    Real &tran_coef, AthenaArray<Real> &ir_cm, AthenaArray<Real> &ir_face,
    AthenaArray<Real> &split_ratio, AthenaArray<int> &map_start,
    AthenaArray<int> &map_end) {
  const int& nfreq = pmy_rad->nfreq;
  Real *nu_lab = &(pmy_rad->nu_grid(0));
  // check to make sure nfreq > 2
  if (nfreq < 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ForwardSplitting]"
        << std::endl << "nfreq '" << nfreq <<
          "' is smaller than 2! ";
    ATHENA_ERROR(msg);
  }
  // map intensity to the desired bin
  // This is a generic function to shift any array
  if (tran_coef >= 1) {
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      Real nu_l = nu_lab[ifr] * tran_coef;
      Real nu_r = nu_lab[ifr+1] * tran_coef;
      int l_bd = ifr;
      int r_bd = ifr;

      while((nu_l > nu_lab[l_bd+1]) && (l_bd < nfreq-1))   l_bd++;
      r_bd = l_bd; // r_bd always > l_bd
      while((nu_r > nu_lab[r_bd+1]) && (r_bd < nfreq-1))   r_bd++;

      if (r_bd-l_bd+1 > nmax_map_) {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [ForwardSplitting]"
          << std::endl << "Frequency shift '" << r_bd-l_bd+1 <<
          "' larger than maximum allowed " << nmax_map_;
        ATHENA_ERROR(msg);
      }
      map_start(ifr) = l_bd;
      map_end(ifr) = r_bd;
      SplitFrequencyBinLinear(l_bd, r_bd, nu_lab, nu_l, nu_r, ir_face(ifr),
                                    ir_face(ifr+1), &(split_ratio(ifr,0)));
    }
    // the last frequency bin
    map_start(nfreq-1) = nfreq-1;
    map_end(nfreq-1) = nfreq-1;
    split_ratio(nfreq-1,0) = 1.0;
  } else {
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      Real nu_l = nu_lab[ifr] * tran_coef;
      Real nu_r = nu_lab[ifr+1] * tran_coef;
      int l_bd = ifr;
      int r_bd = ifr;

      while((nu_r < nu_lab[r_bd]) && (r_bd > 0))   r_bd--;
      l_bd = r_bd; // r_bd always > l_bd
      while((nu_l < nu_lab[l_bd]) && (l_bd > 0))   l_bd--;

      if (r_bd-l_bd+1 > nmax_map_) {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [ForwardSplitting]"
            << std::endl << "Frequency shift '" << r_bd-l_bd+1 <<
          "' larger than maximum allowed " << nmax_map_;
        ATHENA_ERROR(msg);
      }
      map_start(ifr) = l_bd;
      map_end(ifr) = r_bd;
      SplitFrequencyBinLinear(l_bd, r_bd, nu_lab, nu_l, nu_r, ir_face(ifr),
                                    ir_face(ifr+1), &(split_ratio(ifr,0)));
    }
    Real nu_l = nu_lab[nfreq-1] * tran_coef;
    int r_bd = nfreq-1;
    int l_bd = nfreq-2;// it will always be <= current bin
    while((nu_l < nu_lab[l_bd]) && (l_bd > 0))   l_bd--;
        // This frequency bin now maps to l_bd to r_bd
    map_start(nfreq-1) = l_bd;
    map_end(nfreq-1) = r_bd;
      // nu_l/kt
    Real nu_tr = pmy_rad->EffectiveBlackBody(ir_cm(nfreq-1), nu_l);

      // FitIntPlanckFunc is integral _0 to nu_tr
      // the integral we need is 1 - ori_norm
    Real ori_norm = pmy_rad->FitIntPlanckFunc(nu_tr);
    Real div_ori = 0.0;
    if (1.0 - ori_norm > TINY_NUMBER)
      div_ori = 1.0/(1.0 - ori_norm);


      // the first bin
      // the effective temperature 1/T = nu_tr/nu_l
    Real ratio = pmy_rad->FitIntPlanckFunc(nu_tr*nu_lab[l_bd+1]/nu_l);

      // the difference is (1 - ori_norm) - (1 - ratio)
    split_ratio(nfreq-1,0) = (ratio - ori_norm) * div_ori;
    Real sum = split_ratio(nfreq-1,0);

    for (int m=l_bd+1; m<r_bd; ++m) {
      Real ratio_r = pmy_rad->FitIntPlanckFunc(nu_tr*nu_lab[m+1]/nu_l);
      Real ratio_l = pmy_rad->FitIntPlanckFunc(nu_tr*nu_lab[m]/nu_l);
      split_ratio(nfreq-1,m-l_bd) = (ratio_r - ratio_l) * div_ori;
      sum += split_ratio(nfreq-1,m-l_bd);
    }

    split_ratio(nfreq-1,r_bd-l_bd) = 1.0 - sum;
    // in the case that the blackbody tail is too small, we assume
    // power law spectrum shape to determine the ratio
    // Ir=ir_face(\nu/nu_0)^alpha
    // integral is -ir_face*nu_0/(alpha+1) and alpha < -1
    // ratio is (nu_2^alpha+1- nu_1^alpha+1)/nu_0^alpha+1
  }
  return;
}


// interpolate co-moving frame specific intensity over frequency grid
// every bin is
void RadIntegrator::MapIrcmFrequency( AthenaArray<Real> &split_ratio,
      AthenaArray<int> &map_start, AthenaArray<int> &map_end,
      AthenaArray<Real> &input_array, AthenaArray<Real> &shift_array) {
  const int& nfreq = pmy_rad->nfreq;

  // initialize ir_shift to be 0
  shift_array.ZeroClear();

  // check to make sure nfreq > 2
  if (nfreq < 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [MapIrcmFrequency]"
        << std::endl << "nfreq '" << nfreq <<
          "' is smaller than 2! ";
    ATHENA_ERROR(msg);
  }

  // map intensity to the desired bin
  for (int ifr=0; ifr<nfreq; ++ifr) {
    // map shifted intensity to the nu_grid
    const int& fre_start=map_start(ifr);
    const int& fre_end = map_end(ifr);
    for (int m=fre_start; m<=fre_end; ++m) {
      shift_array(m) += input_array(ifr) * split_ratio(ifr,m-fre_start);
    }
  }
  return;
}

// constrcutre the matrix when map between different frequency grids
// Matrix x intensity in shifted grid = intensity in default grid
// get the matrix inversion if possible
bool RadIntegrator::FreMapMatrix(AthenaArray<Real> &split_ratio,
          Real &tran_coef, AthenaArray<int> &map_bin_start,
          AthenaArray<int> &map_bin_end, AthenaArray<int> &map_count,
                                      AthenaArray<Real> &map_matrix) {
  const int& nfreq = pmy_rad->nfreq;

  map_matrix.ZeroClear();
  map_count.ZeroClear();

  bool invertible = true;
  // first, check the diagonal elements are non-zero
  for (int ifr=0; ifr<nfreq; ++ifr) {
    if (split_ratio(ifr,ifr-map_bin_start(ifr)) < TINY_NUMBER) {
      invertible = false;
      return invertible;
    }
  }

  // now construct the frequency map matrix for each angle
  if (tran_coef >= 1) {
    // fre_map_matrix is lower triangle
    // make sure fre_map_matrix(ifr,0) is always the diagonal
    for (int ifr=0; ifr<nfreq; ++ifr) {
      const int& fre_start=map_bin_start(ifr);
      const int& fre_end = map_bin_end(ifr);
      // m is always larger than ifr
      // lower triangle, when store the matrix, we always start from diagonal
      for (int m=fre_start; m<=fre_end; ++m) {
        map_matrix(m,m-ifr) = split_ratio(ifr,m-fre_start);
        map_count(m)++;
      }
    }
  } else {
    // fre_map_matrix is upper triangle
    // make sure fre_map_matrix(ifr,0) is always the diagonal
    for (int ifr=0; ifr<nfreq; ++ifr) {
      const int& fre_start=map_bin_start(ifr);
      const int& fre_end = map_bin_end(ifr);
      // m is always smaller than ifr
      // upper triangle, when store the matrix, we always start from diagonal
      for (int m=fre_start; m<=fre_end; ++m) {
        map_matrix(m,ifr-m) = split_ratio(ifr,m-fre_start);
        map_count(m)++;
      }
    }
  }
  return invertible;
}

// We have the equation map_matrix * shift_array = input_array
// need to calculate inverse_map_matrix * input_array

// if this method cause negative intensity, we will give up
bool RadIntegrator::InverseMapFrequency(
    Real &tran_coef,
    AthenaArray<int> &map_count, AthenaArray<Real> &map_matrix,
    AthenaArray<Real> &input_array, AthenaArray<Real> &shift_array) {
  const int& nfreq = pmy_rad->nfreq;
  // clear zero
  shift_array.ZeroClear();
  // Now invert the matrix split_ratio
  // we need to do this for each angle, all frequency bins
  if (tran_coef >= 1) {
    // map_matrix is a lower triangle matrix
    // first, ifr=0
    shift_array(0) = input_array(0)/map_matrix(0,0);
    for (int ifr=1; ifr<nfreq; ++ifr) {
      shift_array(ifr) = input_array(ifr);
      for (int m=1; m<map_count(ifr); ++m) {
        shift_array(ifr) -= map_matrix(ifr,m) * shift_array(ifr-m);
      }
      shift_array(ifr) /= map_matrix(ifr,0);
      if (shift_array(ifr) < 0.0)
        return false;
    }
  } else {
    // map_matrix is a upper triangle,
    // we need to start from ifr=nfreq-1

    // when flag for fail, set intensity to 0
    shift_array(nfreq-1) = input_array(nfreq-1)/map_matrix(nfreq-1,0);

    for (int ifr=nfreq-2; ifr>=0; --ifr) {
      shift_array(ifr) = input_array(ifr);
      for (int m=1; m<map_count(ifr); ++m) {
        shift_array(ifr) -= map_matrix(ifr,m) * shift_array(ifr+m);
      }
      shift_array(ifr) /= map_matrix(ifr,0);
      if (shift_array(ifr) < 0.0)
        return false;
    }
  }
  return true;
}

// fit a linear line between nu_l and nu_r for ir_l and ir_r
void RadIntegrator::SplitFrequencyBinLinear(
    int &l_bd, int &r_bd, Real *nu_lab, Real &nu_l, Real &nu_r,
    Real &ir_l, Real &ir_r, Real *split_ratio) {
  Real ir_sum = ir_l + ir_r;
  if ((r_bd == l_bd) || (ir_sum < TINY_NUMBER)) {
    if (r_bd == l_bd) {
      split_ratio[0] = 1.0;
    } else {
    // determine the ratio based on the frequency grid only
      Real delta_nu = 1.0/(nu_r - nu_l);
      Real sum = 0.0;
      split_ratio[0] = (nu_lab[l_bd+1]-nu_l)*delta_nu;
      sum=split_ratio[0];
      for (int m=l_bd+1; m<r_bd; ++m) {
        split_ratio[m-l_bd] = (nu_lab[m+1]-nu_lab[m])*delta_nu;
        sum += split_ratio[m-l_bd];
      }
      split_ratio[r_bd-l_bd] = 1.0 - sum;
    }
  } else {
    Real sum = 0.0;
    Real delta_nu = 1.0/(nu_r - nu_l);
    Real nu_ratio = (nu_lab[l_bd+1] - nu_l)*delta_nu;
    Real ir_l_bd1 = ir_l + (ir_r-ir_l)*nu_ratio;
    // between nu_l and nu_lab[l_bd+1]
    split_ratio[0] = ((ir_l_bd1+ir_l)/ir_sum)*nu_ratio;
    sum = split_ratio[0];
    for (int m=l_bd+1; m< r_bd; ++m) {
      nu_ratio = (nu_lab[m+1] - nu_lab[m])*delta_nu;
      Real nu_ratio2 = (nu_lab[m+1] - nu_l) * delta_nu;
      Real nu_ratio1 = (nu_lab[m] - nu_l) * delta_nu;
      Real ir_l_m =  ir_l + (ir_r-ir_l)*nu_ratio1;
      Real ir_l_m1 =  ir_l + (ir_r-ir_l)*nu_ratio2;
      split_ratio[m-l_bd] = ((ir_l_m+ir_l_m1)/ir_sum)*nu_ratio;
      sum += split_ratio[m-l_bd];
    }
    // make sure the sum is always 1
    split_ratio[r_bd-l_bd] = 1.0 - sum;
  }
}

// Ir_cm and ir_shift are both co-moving frame intensities
// Ir_shift is defined in the same frequency grid for all angles
// Ir_cm is defined in the corresponding frequency grid as transoformed
// from lab frame
void RadIntegrator::MapCmToLabFrequency(Real &tran_coef,
                    AthenaArray<Real> &split_ratio,
                    AthenaArray<int> &map_start, AthenaArray<int> &map_end,
                    AthenaArray<Real> &ir_shift, AthenaArray<Real> &ir_cm) {
  //Real *nu_fixed = &(pmy_rad->nu_grid(0));
  AthenaArray<Real> &delta_nu = pmy_rad->delta_nu;
  // now call the function to get value at frequency face
  GetCmMCIntensity(ir_shift, delta_nu, ir_face_);

  BackwardSplitting(tran_coef, ir_shift, ir_face_, split_ratio,
                                     map_start,map_end);
  MapIrcmFrequency(split_ratio, map_start, map_end, ir_shift, ir_cm);
  return;
}

// general function to split any array in the frequency bin [nu_f]
// to the frequency bin [Gamma nu_f]
// ir_last_bin is used to determine the shift in the last frequency bin
// assuming BlackBody spectrum in the last frequency bin
// In other bins, we assume piecewise constant, split each bin according
// to frequency overlap
void RadIntegrator::BackwardSplitting(Real &tran_coef,
                      AthenaArray<Real> &ir_cm, AthenaArray<Real> &ir_face,
                      AthenaArray<Real> &split_ratio,
                      AthenaArray<int> &map_start,AthenaArray<int> &map_end) {
  const int& nfreq = pmy_rad->nfreq;
  Real *nu_lab = &(pmy_rad->nu_grid(0));
  Real *nu_shift = &(nu_shift_(0));
  // check to make sure nfreq > 2
  if (nfreq < 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ForwardSplitting]"
        << std::endl << "nfreq '" << nfreq <<
          "' is smaller than 2! ";
    ATHENA_ERROR(msg);
  }

  for (int ifr=0; ifr<nfreq; ++ifr)
    nu_shift_(ifr) = nu_lab[ifr] * tran_coef;

  // map intensity to the desired bin
  // This is a generic function to shift any array
  if (tran_coef < 1) {
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      Real nu_l = nu_lab[ifr];
      Real nu_r = nu_lab[ifr+1];
      int l_bd = ifr;
      int r_bd = ifr;

      while((nu_l > nu_shift[l_bd+1]) && (l_bd < nfreq-1))   l_bd++;
      r_bd = l_bd; // r_bd always > l_bd
      while((nu_r > nu_shift[r_bd+1]) && (r_bd < nfreq-1))   r_bd++;

      if (r_bd-l_bd+1 > nmax_map_) {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [ForwardSplitting]"
          << std::endl << "Frequency shift '" << r_bd-l_bd+1 <<
          "' larger than maximum allowed " << nmax_map_;
        ATHENA_ERROR(msg);
      }
      map_start(ifr) = l_bd;
      map_end(ifr) = r_bd;
      SplitFrequencyBinLinear(l_bd, r_bd, nu_shift, nu_l, nu_r, ir_face(ifr),
                                    ir_face(ifr+1), &(split_ratio(ifr,0)));
    }
    // the last frequency bin
    map_start(nfreq-1) = nfreq-1;
    map_end(nfreq-1) = nfreq-1;
    split_ratio(nfreq-1,0) = 1.0;
  } else {
    for (int ifr=0; ifr<nfreq-1; ++ifr) {
      Real nu_l = nu_lab[ifr];
      Real nu_r = nu_lab[ifr+1];
      int l_bd = ifr;
      int r_bd = ifr;

      while((nu_r < nu_shift[r_bd]) && (r_bd > 0))   r_bd--;
      l_bd = r_bd; // r_bd always > l_bd
      while((nu_l < nu_shift[l_bd]) && (l_bd > 0))   l_bd--;

      if (r_bd-l_bd+1 > nmax_map_) {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [ForwardSplitting]"
            << std::endl << "Frequency shift '" << r_bd-l_bd+1 <<
          "' larger than maximum allowed " << nmax_map_;
        ATHENA_ERROR(msg);
      }
      map_start(ifr) = l_bd;
      map_end(ifr) = r_bd;

      SplitFrequencyBinLinear(l_bd, r_bd, nu_shift, nu_l, nu_r, ir_face(ifr),
                                    ir_face(ifr+1), &(split_ratio(ifr,0)));
    }

    //-------------------------------------

    Real nu_l = nu_lab[nfreq-1];
    int r_bd = nfreq-1;
    int l_bd = nfreq-2;// it will always be <= current bin
    while((nu_l < nu_shift[l_bd]) && (l_bd > 0))   l_bd--;
    // This frequency bin now maps to l_bd to r_bd
    map_start(nfreq-1) = l_bd;
    map_end(nfreq-1) = r_bd;
    // nu_l/kt
    Real nu_tr = pmy_rad->EffectiveBlackBody(ir_cm(nfreq-1), nu_l);

    // FitIntPlanckFunc is integral _0 to nu_tr
    // the integral we need is 1 - ori_norm
    Real ori_norm = pmy_rad->FitIntPlanckFunc(nu_tr);
    Real div_ori = 0.0;
    if (1.0 - ori_norm > TINY_NUMBER) {
      div_ori = 1.0/(1.0 - ori_norm);
    }
    // the first bin
    // the effective temperature 1/T = nu_tr/nu_l
    Real ratio = pmy_rad->FitIntPlanckFunc(nu_tr*nu_shift[l_bd+1]/nu_l);

      // the difference is (1 - ori_norm) - (1 - ratio)
    split_ratio(nfreq-1,0) = (ratio - ori_norm) * div_ori;
    Real sum = split_ratio(nfreq-1,0);

    for (int m=l_bd+1; m<r_bd; ++m) {
      Real ratio_r = pmy_rad->FitIntPlanckFunc(nu_tr*nu_shift[m+1]/nu_l);
      Real ratio_l = pmy_rad->FitIntPlanckFunc(nu_tr*nu_shift[m]/nu_l);
      split_ratio(nfreq-1,m-l_bd) = (ratio_r - ratio_l) * div_ori;
      sum += split_ratio(nfreq-1,m-l_bd);
    }
    split_ratio(nfreq-1,r_bd-l_bd) = 1.0 - sum;
  }
  return;
}
