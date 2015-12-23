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
//! \file  roe.cpp
//  \brief Roe's linearized Riemann solver.
//
// Computes 1D fluxes using Roe's linearization.  When Roe's method fails because of
// negative density or pressure in the intermediate states, LLF fluxes are used instead.
//
// REFERENCES:
// - P. Roe, "Approximate Riemann solvers, parameter vectors, and difference schemes",
//   JCP, 43, 357 (1981).
//======================================================================================

// C/C++ headers
#include <algorithm>  // max()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../../../athena.hpp"
#include "../../../../athena_arrays.hpp"
#include "../../../hydro.hpp"
#include "../../../eos/eos.hpp"

// this class header
#include "../../hydro_integrator.hpp"

// function to compute eigenvalues and eigenvectors of Roe's matrix A
inline void RoeEigensystem(const Real wroe[], Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVE)], Real left_eigenmatrix[][(NWAVE)]);

// (gamma-1) and isothermal sound speed made global so can be shared with eigensystem
static Real gm1, iso_cs;

void HydroIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[NWAVE],wri[NWAVE],wroe[NWAVE],fl[NWAVE],fr[NWAVE],flxi[NWAVE];
  gm1 = pmy_hydro->peos->GetGamma() - 1.0;
  iso_cs = pmy_hydro->peos->GetIsoSoundSpeed();

  Real coeff[NWAVE];
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real du[NWAVE],a[NWAVE],u[NWAVE];

#pragma simd
  for (int i=il; i<=iu; ++i){

//--- Step 1.  Load L/R states into local variables

    wli[IDN]=wl(IDN,i);
    wli[IVX]=wl(ivx,i);
    wli[IVY]=wl(ivy,i);
    wli[IVZ]=wl(ivz,i);
    if (NON_BAROTROPIC_EOS) wli[IEN]=wl(IEN,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    if (NON_BAROTROPIC_EOS) wri[IEN]=wr(IEN,i);

//--- Step 2.  Compute Roe-averaged data from left- and right-states

    Real sqrtdl = sqrt(wli[IDN]);
    Real sqrtdr = sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN]  = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;

    // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
    // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    Real el,er;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IEN]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
      er = wri[IEN]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
      wroe[IEN] = ((el + wli[IEN])/sqrtdl + (er + wri[IEN])/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute eigenvalues and eigenmatrices using Roe-averaged values

    RoeEigensystem(wroe,ev,rem,lem);

//--- Step 4.  Compute L/R fluxes 

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];

    fl[IDN] = mxl;
    fr[IDN] = mxr;

    fl[IVX] = mxl*wli[IVX];
    fr[IVX] = mxr*wri[IVX];

    fl[IVY] = mxl*wli[IVY];
    fr[IVY] = mxr*wri[IVY];

    fl[IVZ] = mxl*wli[IVZ];
    fr[IVZ] = mxr*wri[IVZ];

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IEN];
      fr[IVX] += wri[IEN];
      fl[IEN] = (el + wli[IEN])*wli[IVX];
      fr[IEN] = (er + wri[IEN])*wri[IVX];
    } else {
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

//--- Step 5.  Compute projection of dU onto L eigenvectors ("vector A")

    du[IDN] = wri[IDN]          - wli[IDN];
    du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
    du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
    du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) du[IEN] = er - el;

    a[IDN]  = lem[IDN][IDN]*du[IDN];
    a[IDN] += lem[IDN][IVX]*du[IVX];
    a[IDN] += lem[IDN][IVY]*du[IVY];
    a[IDN] += lem[IDN][IVZ]*du[IVZ];

    a[IVX]  = lem[IVX][IDN]*du[IDN];
    a[IVX] += lem[IVX][IVX]*du[IVX];
    a[IVX] += lem[IVX][IVY]*du[IVY];
    a[IVX] += lem[IVX][IVZ]*du[IVZ];

    a[IVY]  = lem[IVY][IDN]*du[IDN];
    a[IVY] += lem[IVY][IVX]*du[IVX];
    a[IVY] += lem[IVY][IVY]*du[IVY];
    a[IVY] += lem[IVY][IVZ]*du[IVZ];

    a[IVZ]  = lem[IVZ][IDN]*du[IDN];
    a[IVZ] += lem[IVZ][IVX]*du[IVX];
    a[IVZ] += lem[IVZ][IVY]*du[IVY];
    a[IVZ] += lem[IVZ][IVZ]*du[IVZ];

    if (NON_BAROTROPIC_EOS) {
      a[IDN] += lem[IDN][IEN]*du[IEN];
      a[IVX] += lem[IVX][IEN]*du[IEN];
      a[IVY] += lem[IVY][IEN]*du[IEN];
      a[IVZ] += lem[IVZ][IEN]*du[IEN];

      a[IEN]  = lem[IEN][IDN]*du[IDN];
      a[IEN] += lem[IEN][IVX]*du[IVX];
      a[IEN] += lem[IEN][IVY]*du[IVY];
      a[IEN] += lem[IEN][IVZ]*du[IVZ];
      a[IEN] += lem[IEN][IEN]*du[IEN];
    }

//--- Step 6.  Check that the density and pressure in the intermediate states are
// positive.  If not, set a flag that will be checked below.

    int llf_flag = 0;
    u[IDN] = wli[IDN];
    u[IVX] = wli[IDN]*wli[IVX];
    u[IVY] = wli[IDN]*wli[IVY];
    u[IVZ] = wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) u[IEN] = el;

    // jump across wave[0]
    u[IDN] += a[0]*rem[IDN][0];
    if (u[IDN] < 0.0) llf_flag=1;
    if (NON_BAROTROPIC_EOS) {
      u[IVX] += a[0]*rem[IVX][0];
      u[IVY] += a[0]*rem[IVY][0];
      u[IVZ] += a[0]*rem[IVZ][0];
      u[IEN] += a[0]*rem[IEN][0];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
      if (p < 0.0) llf_flag=2;
    }

    // jump across wave[1]
    u[IDN] += a[1]*rem[IDN][1];
    if (u[IDN] < 0.0) llf_flag=1;
    if (NON_BAROTROPIC_EOS) {
      u[IVX] += a[1]*rem[IVX][1];
      u[IVY] += a[1]*rem[IVY][1];
      u[IVZ] += a[1]*rem[IVZ][1];
      u[IEN] += a[1]*rem[IEN][1];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
      if (p < 0.0) llf_flag=2;
    }

    // jump across wave[2]
    u[IDN] += a[2]*rem[IDN][2];
    if (u[IDN] < 0.0) llf_flag=1;
    if (NON_BAROTROPIC_EOS) {
      u[IVX] += a[2]*rem[IVX][2];
      u[IVY] += a[2]*rem[IVY][2];
      u[IVZ] += a[2]*rem[IVZ][2];
      u[IEN] += a[2]*rem[IEN][2];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
      if (p < 0.0) llf_flag=2;
    }

    if (NON_BAROTROPIC_EOS) {
      // jump across wave[3]
      u[IDN] += a[3]*rem[IDN][3];
      if (u[IDN] < 0.0) llf_flag=1;
      u[IVX] += a[3]*rem[IVX][3];
      u[IVY] += a[3]*rem[IVY][3];
      u[IVZ] += a[3]*rem[IVZ][3];
      u[IEN] += a[3]*rem[IEN][3];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
      if (p < 0.0) llf_flag=2;
    }

//--- Step 7.  Compute Roe flux

    coeff[IDN] = 0.5*fabs(ev[IDN])*a[IDN];
    coeff[IVX] = 0.5*fabs(ev[IVX])*a[IVX];
    coeff[IVY] = 0.5*fabs(ev[IVY])*a[IVY];
    coeff[IVZ] = 0.5*fabs(ev[IVZ])*a[IVZ];

    flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]);
    flxi[IDN] -= coeff[IDN]*rem[IDN][IDN];
    flxi[IDN] -= coeff[IVX]*rem[IDN][IVX];
    flxi[IDN] -= coeff[IVY]*rem[IDN][IVY];
    flxi[IDN] -= coeff[IVZ]*rem[IDN][IVZ];

    flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]);
    flxi[IVX] -= coeff[IDN]*rem[IVX][IDN];
    flxi[IVX] -= coeff[IVX]*rem[IVX][IVX];
    flxi[IVX] -= coeff[IVY]*rem[IVX][IVY];
    flxi[IVX] -= coeff[IVZ]*rem[IVX][IVZ];

    flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]);
    flxi[IVY] -= coeff[IDN]*rem[IVY][IDN];
    flxi[IVY] -= coeff[IVX]*rem[IVY][IVX];
    flxi[IVY] -= coeff[IVY]*rem[IVY][IVY];
    flxi[IVY] -= coeff[IVZ]*rem[IVY][IVZ];

    flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]);
    flxi[IVZ] -= coeff[IDN]*rem[IVZ][IDN];
    flxi[IVZ] -= coeff[IVX]*rem[IVZ][IVX];
    flxi[IVZ] -= coeff[IVY]*rem[IVZ][IVY];
    flxi[IVZ] -= coeff[IVZ]*rem[IVZ][IVZ];

    if (NON_BAROTROPIC_EOS) {
      coeff[IEN] = 0.5*fabs(ev[IEN])*a[IEN];
      flxi[IDN] -= coeff[IEN]*rem[IDN][IEN];
      flxi[IVX] -= coeff[IEN]*rem[IVX][IEN];
      flxi[IVY] -= coeff[IEN]*rem[IVY][IEN];
      flxi[IVZ] -= coeff[IEN]*rem[IVZ][IEN];

      flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]);
      flxi[IEN] -= coeff[IDN]*rem[IEN][IDN];
      flxi[IEN] -= coeff[IVX]*rem[IEN][IVX];
      flxi[IEN] -= coeff[IVY]*rem[IEN][IVY];
      flxi[IEN] -= coeff[IVZ]*rem[IEN][IVZ];
      flxi[IEN] -= coeff[IEN]*rem[IEN][IEN];
    }

//--- Step 8.  Overwrite with upwind flux if flow is supersonic

    if(ev[0] >= 0.0){
      flxi[IDN] = fl[IDN];
      flxi[IVX] = fl[IVX];
      flxi[IVY] = fl[IVY];
      flxi[IVZ] = fl[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fl[IEN];
    }
    if(ev[NWAVE-1] <= 0.0){
      flxi[IDN] = fr[IDN];
      flxi[IVX] = fr[IVX];
      flxi[IVY] = fr[IVY];
      flxi[IVZ] = fr[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fr[IEN];
    }

//--- Step 9.  Overwrite with LLF flux if any of intermediate states are negative

    if (llf_flag != 0) {
      Real a = std::max(fabs(ev[0]), fabs(ev[NWAVE-1]));

      flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]) - a*du[IDN];
      flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]) - a*du[IVX];
      flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]) - a*du[IVY];
      flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]) - a*du[IVZ];
      if (NON_BAROTROPIC_EOS) {
        flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]) - a*du[IEN];
      }
    }

    flx(IDN,i) = flxi[IDN];
    flx(ivx,i) = flxi[IVX];
    flx(ivy,i) = flxi[IVY];
    flx(ivz,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,i) = flxi[IEN];
  }

  return;
}

//--------------------------------------------------------------------------------------
// \!fn RoeEigensystem()
// \brief computes eigenvalues and eigenvectors for hydrodynamics
//
// PURPOSE: Functions to evaluate the eigenvalues, and left- and right-eigenvectors of
// "Roe's matrix A" for the linearized system in the CONSERVED variables, i.e.
// U,t = AU,x, where U=(d,d*vx,d*vy,d*vz,[E],[By,Bz]). The eigenvalues are returned
// through the argument list as a vector of length NWAVE.  The eigenvectors are returned
// as matrices of size (NWAVE)x(NWAVE), with right-eigenvectors stored as COLUMNS
// (so R_i = right_eigenmatrix[*][i]), and left-eigenvectors stored as ROWS
// (so L_i = left_eigenmatrix[i][*]).
//     - Input: v1,v2,v3,h = Roe averaged velocities and enthalpy
//     - Output: eigenvalues[], right_eigenmatrix[][], left_eigenmatrix[][];
//
// REFERENCES:
// - P. Cargo & G. Gallice, "Roe matrices for ideal MHD and systematic construction of
//   Roe matrices for systems of conservation laws", JCP, 136, 446 (1997)
//
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix B  Equation numbers refer to this paper.
//--------------------------------------------------------------------------------------

inline void RoeEigensystem(const Real wroe[], Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVE)], Real left_eigenmatrix[][(NWAVE)])
{
  Real d  = wroe[IDN];
  Real v1 = wroe[IVX];
  Real v2 = wroe[IVY];
  Real v3 = wroe[IVZ];

// Adiabatic hydrodynamics

  if (NON_BAROTROPIC_EOS) {
    Real h = wroe[IEN];
    Real vsq = v1*v1 + v2*v2 + v3*v3;
    Real q = h - 0.5*vsq;
    Real asq = (q < 0.0) ? 0.0 : gm1*q;
    Real a = sqrt(asq);

    // Compute eigenvalues (eq. B2)
    eigenvalues[0] = v1 - a;
    eigenvalues[1] = v1;
    eigenvalues[2] = v1;
    eigenvalues[3] = v1;
    eigenvalues[4] = v1 + a;

    // Right-eigenvectors, stored as COLUMNS (eq. B3)
    right_eigenmatrix[0][0] = 1.0;
    right_eigenmatrix[1][0] = v1 - a;
    right_eigenmatrix[2][0] = v2;
    right_eigenmatrix[3][0] = v3;
    right_eigenmatrix[4][0] = h - v1*a;

    right_eigenmatrix[0][1] = 0.0;
    right_eigenmatrix[1][1] = 0.0;
    right_eigenmatrix[2][1] = 1.0;
    right_eigenmatrix[3][1] = 0.0;
    right_eigenmatrix[4][1] = v2;

    right_eigenmatrix[0][2] = 0.0;
    right_eigenmatrix[1][2] = 0.0;
    right_eigenmatrix[2][2] = 0.0; 
    right_eigenmatrix[3][2] = 1.0;
    right_eigenmatrix[4][2] = v3;

    right_eigenmatrix[0][3] = 1.0;
    right_eigenmatrix[1][3] = v1;
    right_eigenmatrix[2][3] = v2;
    right_eigenmatrix[3][3] = v3;
    right_eigenmatrix[4][3] = 0.5*vsq;

    right_eigenmatrix[0][4] = 1.0;
    right_eigenmatrix[1][4] = v1 + a;
    right_eigenmatrix[2][4] = v2;
    right_eigenmatrix[3][4] = v3;
    right_eigenmatrix[4][4] = h + v1*a;

    // Left-eigenvectors, stored as ROWS (eq. B4)
    Real na = 0.5/asq;
    left_eigenmatrix[0][0] = na*(0.5*gm1*vsq + v1*a);
    left_eigenmatrix[0][1] = -na*(gm1*v1 + a);
    left_eigenmatrix[0][2] = -na*gm1*v2;
    left_eigenmatrix[0][3] = -na*gm1*v3;
    left_eigenmatrix[0][4] = na*gm1;
  
    left_eigenmatrix[1][0] = -v2;
    left_eigenmatrix[1][1] = 0.0;
    left_eigenmatrix[1][2] = 1.0;
    left_eigenmatrix[1][3] = 0.0;
    left_eigenmatrix[1][4] = 0.0;

    left_eigenmatrix[2][0] = -v3;
    left_eigenmatrix[2][1] = 0.0;
    left_eigenmatrix[2][2] = 0.0;
    left_eigenmatrix[2][3] = 1.0;
    left_eigenmatrix[2][4] = 0.0; 

    Real qa = gm1/asq;
    left_eigenmatrix[3][0] = 1.0 - na*gm1*vsq;
    left_eigenmatrix[3][1] = qa*v1;
    left_eigenmatrix[3][2] = qa*v2;
    left_eigenmatrix[3][3] = qa*v3;
    left_eigenmatrix[3][4] = -qa;

    left_eigenmatrix[4][0] = na*(0.5*gm1*vsq - v1*a);
    left_eigenmatrix[4][1] = -na*(gm1*v1 - a);
    left_eigenmatrix[4][2] = left_eigenmatrix[0][2];
    left_eigenmatrix[4][3] = left_eigenmatrix[0][3];
    left_eigenmatrix[4][4] = left_eigenmatrix[0][4];

// Isothermal hydrodynamics

  } else {
    // Compute eigenvalues (eq. B6)
    eigenvalues[0] = v1 - iso_cs;
    eigenvalues[1] = v1;
    eigenvalues[2] = v1;
    eigenvalues[3] = v1 + iso_cs;

    // Right-eigenvectors, stored as COLUMNS (eq. B3)
    right_eigenmatrix[0][0] = 1.0;
    right_eigenmatrix[1][0] = v1 - iso_cs;
    right_eigenmatrix[2][0] = v2;
    right_eigenmatrix[3][0] = v3;

    right_eigenmatrix[0][1] = 0.0;
    right_eigenmatrix[1][1] = 0.0;
    right_eigenmatrix[2][1] = 1.0;
    right_eigenmatrix[3][1] = 0.0;

    right_eigenmatrix[0][2] = 0.0;
    right_eigenmatrix[1][2] = 0.0;
    right_eigenmatrix[2][2] = 0.0;
    right_eigenmatrix[3][2] = 1.0;

    right_eigenmatrix[0][3] = 1.0;
    right_eigenmatrix[1][3] = v1 + iso_cs;
    right_eigenmatrix[2][3] = v2;
    right_eigenmatrix[3][3] = v3;

    // Left-eigenvectors, stored as ROWS (eq. B7)
    left_eigenmatrix[0][0] = 0.5*(1.0 + v1/iso_cs);
    left_eigenmatrix[0][1] = -0.5/iso_cs;
    left_eigenmatrix[0][2] = 0.0;
    left_eigenmatrix[0][3] = 0.0;

    left_eigenmatrix[1][0] = -v2;
    left_eigenmatrix[1][1] = 0.0;
    left_eigenmatrix[1][2] = 1.0;
    left_eigenmatrix[1][3] = 0.0;

    left_eigenmatrix[2][0] = -v3;
    left_eigenmatrix[2][1] = 0.0;
    left_eigenmatrix[2][2] = 0.0;
    left_eigenmatrix[2][3] = 1.0;

    left_eigenmatrix[3][0] = 0.5*(1.0 - v1/iso_cs);
    left_eigenmatrix[3][1] = 0.5/iso_cs;
    left_eigenmatrix[3][2] = 0.0;
    left_eigenmatrix[3][3] = 0.0;
  }
}
