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

// Primary header
#include "../../fluid_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../../athena.hpp"         // enums, macros, Real
#include "../../../../athena_arrays.hpp"  // AthenaArray
#include "../../../fluid.hpp"             // Fluid
#include "../../../eos/eos.hpp"           // GetGamma

// function to compute eigenvalues and eigenvectors of Roe's matrix A
void RoeEigensystem(const Real v1, const Real v2, const Real v3, const Real h,
  Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVE)], Real left_eigenmatrix[][(NWAVE)]);

// (gamma-1) and isothermal sound speed made global so can be shared with eigensystem
static Real gm1, iso_cs;

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

void FluidIntegrator::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[NFLUID],wri[NFLUID],flxi[NFLUID],fl[NFLUID],fr[NFLUID];
  gm1 = pmy_fluid->pf_eos->GetGamma() - 1.0;
  iso_cs = pmy_fluid->pf_eos->GetIsoSoundSpeed();

  Real coeff[NWAVE];
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real du[NWAVE],a[NWAVE],u[NWAVE];

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

    Real droe  = sqrtdl*sqrtdr;
    Real v1roe = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    Real v2roe = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    Real v3roe = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;

    // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
    // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    Real el,er,hroe;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IEN]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
      er = wri[IEN]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
      hroe = ((el + wli[IEN])/sqrtdl + (er + wri[IEN])/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute eigenvalues and eigenmatrices using Roe-averaged values

    RoeEigensystem(v1roe,v2roe,v3roe,hroe,ev,rem,lem);

//--- Step 4.  Compute L/R fluxes 

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];

    fl[IDN]  = mxl;
    fr[IDN]  = mxr;

    fl[IVX] = mxl*wli[IVX];
    fr[IVX] = mxr*wri[IVX];

    fl[IVY] = mxl*wli[IVY];
    fr[IVY] = mxr*wri[IVY];

    fl[IVZ] = mxl*wli[IVZ];
    fr[IVZ] = mxr*wri[IVZ];

    if (NON_BAROTROPIC_EOS) {
      fl[IEN] = (el + wli[IEN])*wli[IVX];
      fr[IEN] = (er + wri[IEN])*wri[IVX];
      fl[IVX] += wli[IEN];
      fr[IVX] += wri[IEN];
    } else {
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

//--- Step 5.  Return upwind flux if flow is supersonic

    if(ev[0] >= 0.0){
      flxi[IDN] = fl[IDN];
      flxi[IVX] = fl[IVX];
      flxi[IVY] = fl[IVY];
      flxi[IVZ] = fl[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fl[IEN];
    } else if(ev[NWAVE-1] <= 0.0){
      flxi[IDN] = fr[IDN];
      flxi[IVX] = fr[IVX];
      flxi[IVY] = fr[IVY];
      flxi[IVZ] = fr[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fr[IEN];
    } else {

//--- Step 6.  Compute projection of dU onto L eigenvectors ("vector A")

      du[IDN] = wri[IDN]          - wli[IDN];
      du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
      du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
      du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
      if (NON_BAROTROPIC_EOS) du[IEN] = er - el;

      a[IDN] = 0.0;
      for (int m=0; m<NWAVE; ++m) a[IDN] += lem[IDN][m]*du[m];

      a[IVX] = 0.0;
      for (int m=0; m<NWAVE; ++m) a[IVX] += lem[IVX][m]*du[m];

      a[IVY] = 0.0;
      for (int m=0; m<NWAVE; ++m) a[IVY] += lem[IVY][m]*du[m];

      a[IVZ] = 0.0;
      for (int m=0; m<NWAVE; ++m) a[IVZ] += lem[IVZ][m]*du[m];

      if (NON_BAROTROPIC_EOS) {
        a[IEN] = 0.0;
        for (int m=0; m<NWAVE; ++m) a[IEN] += lem[IEN][m]*du[m];
      }

//--- Step 7.  Check that the density and pressure in the intermediate states are
// positive.  If not, use LLF fluxes instead.

      int llf_flag = 0;
      u[IDN] = wli[IDN];
      u[IVX] = wli[IDN]*wli[IVX];
      u[IVY] = wli[IDN]*wli[IVY];
      u[IVZ] = wli[IDN]*wli[IVZ];
      if (NON_BAROTROPIC_EOS) u[IEN] = el;

      // jump across wave[0]
      for (int m=0; m<NWAVE; ++m) u[m] += a[0]*rem[m][0];
      if (u[IDN] < 0.0) llf_flag=1;
      if (NON_BAROTROPIC_EOS) {
        Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
        if (p < 0.0) llf_flag=2;
      }

      // jump across wave[1]
      for (int m=0; m<NWAVE; ++m) u[m] += a[1]*rem[m][1];
      if (u[IDN] < 0.0) llf_flag=1;
      if (NON_BAROTROPIC_EOS) {
        Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
        if (p < 0.0) llf_flag=2;
      }

      // jump across wave[2]
      for (int m=0; m<NWAVE; ++m) u[m] += a[2]*rem[m][2];
      if (u[IDN] < 0.0) llf_flag=1;
      if (NON_BAROTROPIC_EOS) {
        Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
        if (p < 0.0) llf_flag=2;
      }

      if (NON_BAROTROPIC_EOS) {
        // jump across wave[3]
        for (int m=0; m<NWAVE; ++m) u[m] += a[3]*rem[m][3];
        if (u[IDN] < 0.0) llf_flag=1;
        if (NON_BAROTROPIC_EOS) {
          Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN];
          if (p < 0.0) llf_flag=2;
        }
      }

      if (llf_flag != 0) {
        Real sl=ev[0];
        Real sr=ev[NWAVE-1];

        flxi[IDN] = (sr*fl[IDN] - sl*fr[IDN] + sl*sr*du[IDN])/(sr - sl);
        flxi[IVX] = (sr*fl[IVX] - sl*fr[IVX] + sl*sr*du[IVX])/(sr - sl);
        flxi[IVY] = (sr*fl[IVY] - sl*fr[IVY] + sl*sr*du[IVY])/(sr - sl);
        flxi[IVZ] = (sr*fl[IVZ] - sl*fr[IVZ] + sl*sr*du[IVZ])/(sr - sl);
        if (NON_BAROTROPIC_EOS) {
          flxi[IEN] = (sr*fl[IEN] - sl*fr[IEN] + sl*sr*du[IEN])/(sr - sl);
        }
      } else {

//--- Step 8.  Compute Roe flux

        for (int m=0; m<NWAVE; ++m) coeff[m] = 0.5*fabs(ev[m])*a[m];

        flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]);
        for (int m=0; m<NWAVE; ++m) flxi[IDN] -= coeff[m]*rem[IDN][m];

        flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]);
        for (int m=0; m<NWAVE; ++m) flxi[IVX] -= coeff[m]*rem[IVX][m];

        flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]);
        for (int m=0; m<NWAVE; ++m) flxi[IVY] -= coeff[m]*rem[IVY][m];

        flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]);
        for (int m=0; m<NWAVE; ++m) flxi[IVZ] -= coeff[m]*rem[IVZ][m];

        if (NON_BAROTROPIC_EOS) {
          flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]);
          for (int m=0; m<NWAVE; ++m) flxi[IEN] -= coeff[m]*rem[IEN][m];
        }
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

void RoeEigensystem(const Real v1, const Real v2, const Real v3, const Real h,
  Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVE)], Real left_eigenmatrix[][(NWAVE)])
{

// Adiabatic hydrodynamics

  if (NON_BAROTROPIC_EOS) {
    Real vsq = v1*v1 + v2*v2 + v3*v3;
    Real asq = gm1*std::max((h-0.5*vsq), TINY_NUMBER);
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
