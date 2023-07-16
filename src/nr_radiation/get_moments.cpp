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
//! \file get_moments.cpp
//  \brief calculate the moments of the radiation field
//======================================================================================

// C headers

// C++ headers

// Athena++ headers

#include "../defs.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "./integrators/rad_integrators.hpp"
#include "./radiation.hpp"

//--------------------------------------------------------------------------------------
// \!fn void CalculateMoment()
// \brief function to create the radiation moments
// calculate the frequency integrated moments of the radiation field
// including the ghost zones
void NRRadiation::CalculateMoment(AthenaArray<Real> &ir_in) {
  Real er, frx, fry, frz, prxx, pryy, przz, prxy, prxz, pryz;
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = pmy_block->block_size.nx2;
  int n3z = pmy_block->block_size.nx3;
  if (n2z > 1) n2z += (2*(NGHOST));
  if (n3z > 1) n3z += (2*(NGHOST));
  Real *weight = &(wmu(0));

  // reset the moment arrays to be zero
  // There are 13 3D arrays
  for (int n=0; n<13; ++n)
    for (int k=0; k<n3z; ++k)
      for (int j=0; j<n2z; ++j) {
          Real *i_mom = &(rad_mom(n,k,j,0));
        for (int i=0; i<n1z; ++i) {
          i_mom[i] = 0.0;
        }
      }

  for (int k=0; k<n3z; ++k) {
    for (int j=0; j<n2z; ++j) {
      for (int i=0; i<n1z; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          er=0.0; frx=0.0; fry=0.0; frz=0.0;
          prxx=0.0; pryy=0.0; przz=0.0; prxy=0.0;
          prxz=0.0; pryz=0.0;
          Real *intensity = &(ir_in(k,j,i,ifr*nang));
          Real *cosx = &(mu(0,k,j,i,0));
          Real *cosy = &(mu(1,k,j,i,0));
          Real *cosz = &(mu(2,k,j,i,0));
          for (int n=0; n<nang; ++n) {
            Real irweight = weight[n] * intensity[n];
            er   += irweight;
            frx  += irweight * cosx[n];
            fry  += irweight * cosy[n];
            frz  += irweight * cosz[n];
            prxx += irweight * cosx[n] * cosx[n];
            pryy += irweight * cosy[n] * cosy[n];
            przz += irweight * cosz[n] * cosz[n];
            prxy += irweight * cosx[n] * cosy[n];
            prxz += irweight * cosx[n] * cosz[n];
            pryz += irweight * cosy[n] * cosz[n];
          }
          // assign the moments for each frequency group
          if (nfreq > 1) {
            rad_mom_nu(ifr*13+IER,k,j,i) = er;
            rad_mom_nu(ifr*13+IFR1,k,j,i) = frx;
            rad_mom_nu(ifr*13+IFR2,k,j,i) = fry;
            rad_mom_nu(ifr*13+IFR3,k,j,i) = frz;
            rad_mom_nu(ifr*13+IPR11,k,j,i) = prxx;
            rad_mom_nu(ifr*13+IPR22,k,j,i) = pryy;
            rad_mom_nu(ifr*13+IPR33,k,j,i) = przz;
            rad_mom_nu(ifr*13+IPR12,k,j,i) = prxy;
            rad_mom_nu(ifr*13+IPR13,k,j,i) = prxz;
            rad_mom_nu(ifr*13+IPR23,k,j,i) = pryz;
            rad_mom_nu(ifr*13+IPR21,k,j,i) = prxy;
            rad_mom_nu(ifr*13+IPR31,k,j,i) = prxz;
            rad_mom_nu(ifr*13+IPR32,k,j,i) = pryz;
          }

          // assign the moments
          rad_mom(IER,k,j,i) += er;
          rad_mom(IFR1,k,j,i) += frx;
          rad_mom(IFR2,k,j,i) += fry;
          rad_mom(IFR3,k,j,i) += frz;
          rad_mom(IPR11,k,j,i) += prxx;
          rad_mom(IPR22,k,j,i) += pryy;
          rad_mom(IPR33,k,j,i) += przz;
          rad_mom(IPR12,k,j,i) += prxy;
          rad_mom(IPR13,k,j,i) += prxz;
          rad_mom(IPR23,k,j,i) += pryz;
          rad_mom(IPR21,k,j,i) += prxy;
          rad_mom(IPR31,k,j,i) += prxz;
          rad_mom(IPR32,k,j,i) += pryz;
        }
      }
    }
  }
}


//--------------------------------------------------------------------------------------
// \!fn void CalculateComMoment()

// \brief Calculate the radiation moments in the co-moving frame
// Also load specific intensity for dump

void NRRadiation::CalculateComMoment() {
  Hydro *phydro=pmy_block->phydro;
  Real invcrat = 1.0/crat;
  Real er, frx, fry, frz;
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = pmy_block->block_size.nx2;
  int n3z = pmy_block->block_size.nx3;
  if (n2z > 1) n2z += (2*(NGHOST));
  if (n3z > 1) n3z += (2*(NGHOST));

  AthenaArray<Real> &i_mom = rad_mom_cm;
  Real *weight = &(wmu(0));
  Real *ir_output;
  // Get the temporary array
  AthenaArray<Real> &wmu_cm = pradintegrator->wmu_cm_;
  AthenaArray<Real> &tran_coef = pradintegrator->tran_coef_;
  AthenaArray<Real> &ir_cm = pradintegrator->ir_cm_;
  AthenaArray<Real> &cm_to_lab = pradintegrator->cm_to_lab_;
  AthenaArray<Real> &cosx_cm = cosx_cm_;
  AthenaArray<Real> &cosy_cm = cosy_cm_;
  AthenaArray<Real> &cosz_cm = cosz_cm_;

  AthenaArray<Real> split_ratio;
  AthenaArray<int> map_start, map_end;

  // reset the moment arrays to be zero
  // There are 4 3D arrays
  for (int n=0; n<4; ++n) {
    for (int k=0; k<n3z; ++k) {
      for (int j=0; j<n2z; ++j) {
//#pragma omp simd
        for (int i=0; i<n1z; ++i) {
          i_mom(n,k,j,i) = 0.0;
        }
      }
    }
  }

  for (int k=0; k<n3z; ++k) {
    for (int j=0; j<n2z; ++j) {
      for (int i=0; i<n1z; ++i) {
        Real *cosx = &(mu(0,k,j,i,0));
        Real *cosy = &(mu(1,k,j,i,0));
        Real *cosz = &(mu(2,k,j,i,0));

        Real vx = phydro->u(IM1,k,j,i)/phydro->u(IDN,k,j,i);
        Real vy = phydro->u(IM2,k,j,i)/phydro->u(IDN,k,j,i);
        Real vz = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);
        Real vsq = vx * vx + vy * vy + vz * vz;

        Real ratio = std::sqrt(vsq) * invcrat;
         // Limit the velocity to be smaller than the speed of light
        if (ratio > vmax) {
          Real factor = vmax/ratio;
          vx *= factor;
          vy *= factor;
          vz *= factor;

          vsq *= (factor*factor);
        }

            // square of Lorentz factor
        Real lorzsq = 1.0/(1.0 - vsq * invcrat * invcrat);
        Real lorz = std::sqrt(lorzsq);

        // first, get co-moving frame ir_cm
        Real numsum = 0.0;
        for (int n=0; n<nang; ++n) {
          Real vdotn = vx * cosx[n] + vy * cosy[n] + vz * cosz[n];
          Real vnc = 1.0 - vdotn * invcrat;
          tran_coef(n) = std::sqrt(lorzsq) * vnc;
          wmu_cm(n) = weight[n]/(lorzsq * vnc * vnc);
          numsum += wmu_cm(n);
          cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
          Real angcoef = lorz * invcrat * (1.0
                         - lorz * vdotn * invcrat/(1.0 + lorz));
          Real incoef = 1.0/(lorz * vnc);
          cosx_cm(n) = (cosx[n] - vx * angcoef) * incoef;
          cosy_cm(n) = (cosy[n] - vy * angcoef) * incoef;
          cosz_cm(n) = (cosz[n] - vz * angcoef) * incoef;
        }
        numsum = 1.0/numsum;
//#pragma omp simd
        for (int n=0; n<nang; ++n) {
           wmu_cm(n) *= numsum;
        }

        for (int ifr=0; ifr<nfreq; ++ifr) {
          for (int n=0; n<nang; ++n) {
            ir_cm(n+ifr*nang) = ir(k,j,i,ifr*nang+n) * cm_to_lab(n);
          }
        }
        // in the case of multi-groups
        if (nfreq > 1) {
          // shift intensity from shifted frequency bins

          // map frequency grid
          for (int n=0; n<nang; ++n) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              pradintegrator->ir_ori_(ifr) = ir_cm(ifr*nang+n);
            }
            split_ratio.InitWithShallowSlice(pradintegrator->split_ratio_,
                                                                    3,n,1);
            map_start.InitWithShallowSlice(pradintegrator->map_bin_start_,
                                                                    2,n,1);
            map_end.InitWithShallowSlice(pradintegrator->map_bin_end_,
                                                                    2,n,1);


            pradintegrator->MapLabToCmFrequency(tran_coef(n), split_ratio,
                      map_start, map_end,
                      pradintegrator->ir_ori_, pradintegrator->ir_done_);

            for (int ifr=0; ifr<nfreq; ++ifr)
              ir_cm(ifr*nang+n) = pradintegrator->ir_done_(ifr);
          }
        }

        Real *cm_weight = &(wmu_cm(0));
        for (int ifr=0; ifr<nfreq; ++ifr) {
          ir_output = &(ir_cm(ifr*nang));
          er=0.0; frx=0.0; fry=0.0; frz=0.0;
          for (int n=0; n<nang; ++n) {
            Real ir_weight = ir_output[n] * cm_weight[n];
            er   += ir_weight;
            frx  += ir_weight * cosx_cm(n);
            fry  += ir_weight * cosy_cm(n);
            frz  += ir_weight * cosz_cm(n);
          }

          // assign the moments for each frequency group
          if (nfreq > 1) {
            rad_mom_cm_nu(ifr*4+IER,k,j,i) = er;
            rad_mom_cm_nu(ifr*4+IFR1,k,j,i) = frx;
            rad_mom_cm_nu(ifr*4+IFR2,k,j,i) = fry;
            rad_mom_cm_nu(ifr*4+IFR3,k,j,i) = frz;
          }
          i_mom(IER,k,j,i) += er;
          i_mom(IFR1,k,j,i) += frx;
          i_mom(IFR2,k,j,i) += fry;
          i_mom(IFR3,k,j,i) += frz;
        }

        // prepare the opacity array for output
        for (int ifr=0; ifr<nfreq; ++ifr) {
          output_sigma(3*ifr+OPAS,k,j,i) = sigma_s(k,j,i,ifr);
          output_sigma(3*ifr+OPAA,k,j,i) = sigma_a(k,j,i,ifr);
          output_sigma(3*ifr+OPAP,k,j,i) = sigma_p(k,j,i,ifr);
        }
      }
    }
  }
  return;
}
