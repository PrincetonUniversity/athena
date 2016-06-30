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
//! \file  roe_mhd.cpp
//  \brief Roe's linearized Riemann solver for MHD.
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
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"

// this class header
#include "../../hydro.hpp"

// function to compute eigenvalues and eigenvectors of Roe's matrix A
inline static void RoeEigensystem(const Real wroe[], const Real b1, 
  const Real x, const Real y, Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVE)], Real left_eigenmatrix[][(NWAVE)]);

// (gamma-1) and isothermal sound speed made global so can be shared with eigensystem
static Real gm1, iso_cs;

void Hydro::RiemannSolver(const int k,const int j, const int il, const int iu,
  const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;
  Real wli[NWAVE],wri[NWAVE],wroe[NWAVE],fl[NWAVE],fr[NWAVE],flxi[NWAVE];
  gm1 = pmy_block->peos->GetGamma() - 1.0;
  iso_cs = pmy_block->peos->GetIsoSoundSpeed();

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
    if (NON_BAROTROPIC_EOS) wli[IPR]=wl(IPR,i);
    wli[IBY]=wl(IBY,i);
    wli[IBZ]=wl(IBZ,i);

    wri[IDN]=wr(IDN,i);
    wri[IVX]=wr(ivx,i);
    wri[IVY]=wr(ivy,i);
    wri[IVZ]=wr(ivz,i);
    if (NON_BAROTROPIC_EOS) wri[IPR]=wr(IPR,i);
    wri[IBY]=wr(IBY,i);
    wri[IBZ]=wr(IBZ,i);

    Real bxi = bx(k,j,i);

//--- Step 2.  Compute Roe-averaged data from left- and right-states

    Real sqrtdl = sqrt(wli[IDN]);
    Real sqrtdr = sqrt(wri[IDN]);
    Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

    wroe[IDN] = sqrtdl*sqrtdr;
    wroe[IVX] = (sqrtdl*wli[IVX] + sqrtdr*wri[IVX])*isdlpdr;
    wroe[IVY] = (sqrtdl*wli[IVY] + sqrtdr*wri[IVY])*isdlpdr;
    wroe[IVZ] = (sqrtdl*wli[IVZ] + sqrtdr*wri[IVZ])*isdlpdr;
    // Note Roe average of magnetic field is different
    wroe[IBY] = (sqrtdr*wli[IBY] + sqrtdl*wri[IBY])*isdlpdr;
    wroe[IBZ] = (sqrtdr*wli[IBZ] + sqrtdl*wri[IBZ])*isdlpdr;
    Real x = 0.5*(SQR(wli[IBY]-wri[IBY]) + SQR(wli[IBZ]-wri[IBZ]))/(SQR(sqrtdl+sqrtdr));
    Real y = 0.5*(wli[IDN] + wri[IDN])/wroe[IDN];

    // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
    // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    Real pbl = 0.5*(bxi*bxi + SQR(wli[IBY]) + SQR(wli[IBZ]));
    Real pbr = 0.5*(bxi*bxi + SQR(wri[IBY]) + SQR(wri[IBZ]));
    Real el,er;
    if (NON_BAROTROPIC_EOS) {
      el = wli[IPR]/gm1 + 0.5*wli[IDN]*(SQR(wli[IVX])+SQR(wli[IVY])+SQR(wli[IVZ])) +pbl;
      er = wri[IPR]/gm1 + 0.5*wri[IDN]*(SQR(wri[IVX])+SQR(wri[IVY])+SQR(wri[IVZ])) +pbr;
      wroe[IPR] = ((el + wli[IPR] + pbl)/sqrtdl + (er + wri[IPR] + pbr)/sqrtdr)*isdlpdr;
    }

//--- Step 3.  Compute eigenvalues and eigenmatrices using Roe-averaged values

    RoeEigensystem(wroe,bxi,x,y,ev,rem,lem);

//--- Step 4.  Compute L/R fluxes 

    Real mxl = wli[IDN]*wli[IVX];
    Real mxr = wri[IDN]*wri[IVX];

    fl[IDN] = mxl;
    fr[IDN] = mxr;

    fl[IVX] = mxl*wli[IVX] + pbl - SQR(bxi);
    fr[IVX] = mxr*wri[IVX] + pbr - SQR(bxi);

    fl[IVY] = mxl*wli[IVY] - bxi*wli[IBY];
    fr[IVY] = mxr*wri[IVY] - bxi*wri[IBY];

    fl[IVZ] = mxl*wli[IVZ] - bxi*wli[IBZ];
    fr[IVZ] = mxr*wri[IVZ] - bxi*wri[IBZ];

    if (NON_BAROTROPIC_EOS) {
      fl[IVX] += wli[IPR];
      fr[IVX] += wri[IPR];
      fl[IEN] = (el + wli[IPR] + pbl - bxi*bxi)*wli[IVX];
      fr[IEN] = (er + wri[IPR] + pbr - bxi*bxi)*wri[IVX];
      fl[IEN] -= bxi*(wli[IBY]*wli[IVY] + wli[IBZ]*wli[IVZ]);
      fr[IEN] -= bxi*(wri[IBY]*wri[IVY] + wri[IBZ]*wri[IVZ]);
    } else {
      fl[IVX] += (iso_cs*iso_cs)*wli[IDN];
      fr[IVX] += (iso_cs*iso_cs)*wri[IDN];
    }

    fl[IBY] = wli[IBY]*wli[IVX] - bxi*wli[IVY];
    fr[IBY] = wri[IBY]*wri[IVX] - bxi*wri[IVY];

    fl[IBZ] = wli[IBZ]*wli[IVX] - bxi*wli[IVZ];
    fr[IBZ] = wri[IBZ]*wri[IVX] - bxi*wri[IVZ];

//--- Step 5.  Compute projection of dU onto L eigenvectors ("vector A")

    du[IDN] = wri[IDN]          - wli[IDN];
    du[IVX] = wri[IDN]*wri[IVX] - wli[IDN]*wli[IVX];
    du[IVY] = wri[IDN]*wri[IVY] - wli[IDN]*wli[IVY];
    du[IVZ] = wri[IDN]*wri[IVZ] - wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) du[IEN] = er - el;
    du[IBY] = wri[IBY] - wli[IBY];
    du[IBZ] = wri[IBZ] - wli[IBZ];

    a[IDN]  = lem[IDN][IDN]*du[IDN];
    a[IDN] += lem[IDN][IVX]*du[IVX];
    a[IDN] += lem[IDN][IVY]*du[IVY];
    a[IDN] += lem[IDN][IVZ]*du[IVZ];
    a[IDN] += lem[IDN][IBY]*du[IBY];
    a[IDN] += lem[IDN][IBZ]*du[IBZ];

    a[IVX]  = lem[IVX][IDN]*du[IDN];
    a[IVX] += lem[IVX][IVX]*du[IVX];
    a[IVX] += lem[IVX][IVY]*du[IVY];
    a[IVX] += lem[IVX][IVZ]*du[IVZ];
    a[IVX] += lem[IVX][IBY]*du[IBY];
    a[IVX] += lem[IVX][IBZ]*du[IBZ];

    a[IVY]  = lem[IVY][IDN]*du[IDN];
    a[IVY] += lem[IVY][IVX]*du[IVX];
    a[IVY] += lem[IVY][IVY]*du[IVY];
    a[IVY] += lem[IVY][IVZ]*du[IVZ];
    a[IVY] += lem[IVY][IBY]*du[IBY];
    a[IVY] += lem[IVY][IBZ]*du[IBZ];

    a[IVZ]  = lem[IVZ][IDN]*du[IDN];
    a[IVZ] += lem[IVZ][IVX]*du[IVX];
    a[IVZ] += lem[IVZ][IVY]*du[IVY];
    a[IVZ] += lem[IVZ][IVZ]*du[IVZ];
    a[IVZ] += lem[IVZ][IBY]*du[IBY];
    a[IVZ] += lem[IVZ][IBZ]*du[IBZ];

    a[IBY]  = lem[IBY][IDN]*du[IDN];
    a[IBY] += lem[IBY][IVX]*du[IVX];
    a[IBY] += lem[IBY][IVY]*du[IVY];
    a[IBY] += lem[IBY][IVZ]*du[IVZ];
    a[IBY] += lem[IBY][IBY]*du[IBY];
    a[IBY] += lem[IBY][IBZ]*du[IBZ];

    a[IBZ]  = lem[IBZ][IDN]*du[IDN];
    a[IBZ] += lem[IBZ][IVX]*du[IVX];
    a[IBZ] += lem[IBZ][IVY]*du[IVY];
    a[IBZ] += lem[IBZ][IVZ]*du[IVZ];
    a[IBZ] += lem[IBZ][IBY]*du[IBY];
    a[IBZ] += lem[IBZ][IBZ]*du[IBZ];

    if (NON_BAROTROPIC_EOS) {
      a[IDN] += lem[IDN][IEN]*du[IEN];
      a[IVX] += lem[IVX][IEN]*du[IEN];
      a[IVY] += lem[IVY][IEN]*du[IEN];
      a[IVZ] += lem[IVZ][IEN]*du[IEN];
      a[IBY] += lem[IBY][IEN]*du[IEN];
      a[IBZ] += lem[IBZ][IEN]*du[IEN];

      a[IEN]  = lem[IEN][IDN]*du[IDN];
      a[IEN] += lem[IEN][IVX]*du[IVX];
      a[IEN] += lem[IEN][IVY]*du[IVY];
      a[IEN] += lem[IEN][IVZ]*du[IVZ];
      a[IEN] += lem[IEN][IEN]*du[IEN];
      a[IEN] += lem[IEN][IBY]*du[IBY];
      a[IEN] += lem[IEN][IBZ]*du[IBZ];
    }

//--- Step 6.  Check that the density and pressure in the intermediate states are
// positive.  If not, set a flag that will be checked below.

    int llf_flag = 0;
    u[IDN] = wli[IDN];
    u[IVX] = wli[IDN]*wli[IVX];
    u[IVY] = wli[IDN]*wli[IVY];
    u[IVZ] = wli[IDN]*wli[IVZ];
    if (NON_BAROTROPIC_EOS) u[IEN] = el;
    u[IBY] = wli[IBY];
    u[IBZ] = wli[IBZ];

    // jump across wave[0]
    u[IDN] += a[0]*rem[IDN][0];
    if (u[IDN] < 0.0) llf_flag=1;
    if (NON_BAROTROPIC_EOS) {
      u[IVX] += a[0]*rem[IVX][0];
      u[IVY] += a[0]*rem[IVY][0];
      u[IVZ] += a[0]*rem[IVZ][0];
      u[IEN] += a[0]*rem[IEN][0];
      u[IBY] += a[0]*rem[IBY][0];
      u[IBZ] += a[0]*rem[IBZ][0];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN]
                      - 0.5*(SQR(bxi)+SQR(u[IBY])+SQR(u[IBZ]));
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
      u[IBY] += a[1]*rem[IBY][1];
      u[IBZ] += a[1]*rem[IBZ][1];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN]
                      - 0.5*(SQR(bxi)+SQR(u[IBY])+SQR(u[IBZ]));
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
      u[IBY] += a[2]*rem[IBY][2];
      u[IBZ] += a[2]*rem[IBZ][2];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN]
                      - 0.5*(SQR(bxi)+SQR(u[IBY])+SQR(u[IBZ]));
      if (p < 0.0) llf_flag=2;
    }

    // jump across wave[3]
    u[IDN] += a[3]*rem[IDN][3];
    if (u[IDN] < 0.0) llf_flag=1;
    if (NON_BAROTROPIC_EOS) {
      u[IVX] += a[3]*rem[IVX][3];
      u[IVY] += a[3]*rem[IVY][3];
      u[IVZ] += a[3]*rem[IVZ][3];
      u[IEN] += a[3]*rem[IEN][3];
      u[IBY] += a[3]*rem[IBY][3];
      u[IBZ] += a[3]*rem[IBZ][3];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN]
                      - 0.5*(SQR(bxi)+SQR(u[IBY])+SQR(u[IBZ]));
      if (p < 0.0) llf_flag=2;
    }

    // jump across wave[4]
    u[IDN] += a[4]*rem[IDN][4];
    if (u[IDN] < 0.0) llf_flag=1;
    if (NON_BAROTROPIC_EOS) {
      u[IVX] += a[4]*rem[IVX][4];
      u[IVY] += a[4]*rem[IVY][4];
      u[IVZ] += a[4]*rem[IVZ][4];
      u[IEN] += a[4]*rem[IEN][4];
      u[IBY] += a[4]*rem[IBY][4];
      u[IBZ] += a[4]*rem[IBZ][4];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN]
                      - 0.5*(SQR(bxi)+SQR(u[IBY])+SQR(u[IBZ]));
      if (p < 0.0) llf_flag=2;
    }

    if (NON_BAROTROPIC_EOS) {
      // jump across wave[5]
      u[IDN] += a[5]*rem[IDN][5];
      if (u[IDN] < 0.0) llf_flag=1;
      u[IVX] += a[5]*rem[IVX][5];
      u[IVY] += a[5]*rem[IVY][5];
      u[IVZ] += a[5]*rem[IVZ][5];
      u[IEN] += a[5]*rem[IEN][5];
      u[IBY] += a[5]*rem[IBY][5];
      u[IBZ] += a[5]*rem[IBZ][5];
      Real p = u[IEN] - 0.5*(SQR(u[IVX])+SQR(u[IVY])+SQR(u[IVZ]))/u[IDN]
                      - 0.5*(SQR(bxi)+SQR(u[IBY])+SQR(u[IBZ]));
      if (p < 0.0) llf_flag=2;
    }

//--- Step 7.  Compute Roe flux

    coeff[IDN] = 0.5*fabs(ev[IDN])*a[IDN];
    coeff[IVX] = 0.5*fabs(ev[IVX])*a[IVX];
    coeff[IVY] = 0.5*fabs(ev[IVY])*a[IVY];
    coeff[IVZ] = 0.5*fabs(ev[IVZ])*a[IVZ];
    coeff[IBY] = 0.5*fabs(ev[IBY])*a[IBY];
    coeff[IBZ] = 0.5*fabs(ev[IBZ])*a[IBZ];

    flxi[IDN] = 0.5*(fl[IDN] + fr[IDN]);
    flxi[IDN] -= coeff[IDN]*rem[IDN][IDN];
    flxi[IDN] -= coeff[IVX]*rem[IDN][IVX];
    flxi[IDN] -= coeff[IVY]*rem[IDN][IVY];
    flxi[IDN] -= coeff[IVZ]*rem[IDN][IVZ];
    flxi[IDN] -= coeff[IBY]*rem[IDN][IBY];
    flxi[IDN] -= coeff[IBZ]*rem[IDN][IBZ];

    flxi[IVX] = 0.5*(fl[IVX] + fr[IVX]);
    flxi[IVX] -= coeff[IDN]*rem[IVX][IDN];
    flxi[IVX] -= coeff[IVX]*rem[IVX][IVX];
    flxi[IVX] -= coeff[IVY]*rem[IVX][IVY];
    flxi[IVX] -= coeff[IVZ]*rem[IVX][IVZ];
    flxi[IVX] -= coeff[IBY]*rem[IVX][IBY];
    flxi[IVX] -= coeff[IBZ]*rem[IVX][IBZ];

    flxi[IVY] = 0.5*(fl[IVY] + fr[IVY]);
    flxi[IVY] -= coeff[IDN]*rem[IVY][IDN];
    flxi[IVY] -= coeff[IVX]*rem[IVY][IVX];
    flxi[IVY] -= coeff[IVY]*rem[IVY][IVY];
    flxi[IVY] -= coeff[IVZ]*rem[IVY][IVZ];
    flxi[IVY] -= coeff[IBY]*rem[IVY][IBY];
    flxi[IVY] -= coeff[IBZ]*rem[IVY][IBZ];

    flxi[IVZ] = 0.5*(fl[IVZ] + fr[IVZ]);
    flxi[IVZ] -= coeff[IDN]*rem[IVZ][IDN];
    flxi[IVZ] -= coeff[IVX]*rem[IVZ][IVX];
    flxi[IVZ] -= coeff[IVY]*rem[IVZ][IVY];
    flxi[IVZ] -= coeff[IVZ]*rem[IVZ][IVZ];
    flxi[IVZ] -= coeff[IBY]*rem[IVZ][IBY];
    flxi[IVZ] -= coeff[IBZ]*rem[IVZ][IBZ];

    flxi[IBY] = 0.5*(fl[IBY] + fr[IBY]);
    flxi[IBY] -= coeff[IDN]*rem[IBY][IDN];
    flxi[IBY] -= coeff[IVX]*rem[IBY][IVX];
    flxi[IBY] -= coeff[IVY]*rem[IBY][IVY];
    flxi[IBY] -= coeff[IVZ]*rem[IBY][IVZ];
    flxi[IBY] -= coeff[IBY]*rem[IBY][IBY];
    flxi[IBY] -= coeff[IBZ]*rem[IBY][IBZ];

    flxi[IBZ] = 0.5*(fl[IBZ] + fr[IBZ]);
    flxi[IBZ] -= coeff[IDN]*rem[IBZ][IDN];
    flxi[IBZ] -= coeff[IVX]*rem[IBZ][IVX];
    flxi[IBZ] -= coeff[IVY]*rem[IBZ][IVY];
    flxi[IBZ] -= coeff[IVZ]*rem[IBZ][IVZ];
    flxi[IBZ] -= coeff[IBY]*rem[IBZ][IBY];
    flxi[IBZ] -= coeff[IBZ]*rem[IBZ][IBZ];

    if (NON_BAROTROPIC_EOS) {
      coeff[IEN] = 0.5*fabs(ev[IEN])*a[IEN];
      flxi[IDN] -= coeff[IEN]*rem[IDN][IEN];
      flxi[IVX] -= coeff[IEN]*rem[IVX][IEN];
      flxi[IVY] -= coeff[IEN]*rem[IVY][IEN];
      flxi[IVZ] -= coeff[IEN]*rem[IVZ][IEN];
      flxi[IBY] -= coeff[IEN]*rem[IBY][IEN];
      flxi[IBZ] -= coeff[IEN]*rem[IBZ][IEN];

      flxi[IEN] = 0.5*(fl[IEN] + fr[IEN]);
      flxi[IEN] -= coeff[IDN]*rem[IEN][IDN];
      flxi[IEN] -= coeff[IVX]*rem[IEN][IVX];
      flxi[IEN] -= coeff[IVY]*rem[IEN][IVY];
      flxi[IEN] -= coeff[IVZ]*rem[IEN][IVZ];
      flxi[IEN] -= coeff[IEN]*rem[IEN][IEN];
      flxi[IEN] -= coeff[IBY]*rem[IEN][IBY];
      flxi[IEN] -= coeff[IBZ]*rem[IEN][IBZ];
    }

//--- Step 8.  Overwrite with upwind flux if flow is supersonic

    if(ev[0] >= 0.0){
      flxi[IDN] = fl[IDN];
      flxi[IVX] = fl[IVX];
      flxi[IVY] = fl[IVY];
      flxi[IVZ] = fl[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fl[IEN];
      flxi[IBY] = fl[IBY];
      flxi[IBZ] = fl[IBZ];
    }
    if(ev[NWAVE-1] <= 0.0){
      flxi[IDN] = fr[IDN];
      flxi[IVX] = fr[IVX];
      flxi[IVY] = fr[IVY];
      flxi[IVZ] = fr[IVZ];
      if (NON_BAROTROPIC_EOS) flxi[IEN] = fr[IEN];
      flxi[IBY] = fr[IBY];
      flxi[IBZ] = fr[IBZ];
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
      flxi[IBY] = 0.5*(fl[IBY] + fr[IBY]) - a*du[IBY];
      flxi[IBZ] = 0.5*(fl[IBZ] + fr[IBZ]) - a*du[IBZ];
    }

    flx(IDN,i) = flxi[IDN];
    flx(ivx,i) = flxi[IVX];
    flx(ivy,i) = flxi[IVY];
    flx(ivz,i) = flxi[IVZ];
    if (NON_BAROTROPIC_EOS) flx(IEN,i) = flxi[IEN];
    flx(IBY,i) = flxi[IBY];
    flx(IBZ,i) = flxi[IBZ];
  }

  return;
}

//--------------------------------------------------------------------------------------
// \!fn RoeEigensystem()
// \brief computes eigenvalues and eigenvectors for MHD
//
// PURPOSE: Functions to evaluate the eigenvalues, and left- and right-eigenvectors of
// "Roe's matrix A" for the linearized system in the CONSERVED variables, i.e.
// U,t = AU,x, where U=(d,d*vx,d*vy,d*vz,[E],[By,Bz]). The eigenvalues are returned
// through the argument list as a vector of length NWAVE.  The eigenvectors are returned
// as matrices of size (NWAVE)x(NWAVE), with right-eigenvectors stored as COLUMNS
// (so R_i = right_eigenmatrix[*][i]), and left-eigenvectors stored as ROWS
// (so L_i = left_eigenmatrix[i][*]).
//     - Input: d,v1,v2,v3,h,b1,b2,b3=Roe averaged density, velocities, enthalpy, B
//          x,y = numerical factors (see eqn XX)
//     - Output: eigenvalues[], right_eigenmatrix[][], left_eigenmatrix[][];
//
// REFERENCES:
// - P. Cargo & G. Gallice, "Roe matrices for ideal MHD and systematic construction of
//   Roe matrices for systems of conservation laws", JCP, 136, 446 (1997)
//
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix B  Equation numbers refer to this paper.
//--------------------------------------------------------------------------------------

inline static void RoeEigensystem(const Real wroe[], const Real b1, 
  const Real x, const Real y, Real eigenvalues[],
  Real right_eigenmatrix[][(NWAVE)], Real left_eigenmatrix[][(NWAVE)])
{
  Real d  = wroe[IDN];
  Real v1 = wroe[IVX];
  Real v2 = wroe[IVY];
  Real v3 = wroe[IVZ];
  Real b2 = wroe[IBY];
  Real b3 = wroe[IBZ];

// Adiabatic MHD

  if (NON_BAROTROPIC_EOS) {
    Real vsq = v1*v1 + v2*v2 + v3*v3;
    Real btsq = b2*b2 + b3*b3;
    Real bt_starsq = (gm1 - (gm1 - 1.0)*y)*btsq;
    Real vaxsq = b1*b1/d;
    Real hp = wroe[IPR] - (vaxsq + btsq/d);
    Real twid_asq = std::max((gm1*(hp-0.5*vsq)-(gm1-1.0)*x), TINY_NUMBER);

    // Compute fast- and slow-magnetosonic speeds (eq. B18)
    Real ct2 = bt_starsq/d;
    Real tsum = vaxsq + ct2 + twid_asq;
    Real tdif = vaxsq + ct2 - twid_asq;
    Real cf2_cs2 = sqrt(tdif*tdif + 4.0*twid_asq*ct2);

    Real cfsq = 0.5*(tsum + cf2_cs2);
    Real cf = sqrt(cfsq);

    Real cssq = twid_asq*vaxsq/cfsq;
    Real cs = sqrt(cssq);

    // Compute beta(s) (eqs. A17, B20, B28)
    Real bt = sqrt(btsq);
    Real bt_star = sqrt(bt_starsq);
    Real bet2=1.0;
    Real bet3=0.0;
    if (bt != 0.0) {
      bet2 = b2/bt;
      bet3 = b3/bt;
    }
    Real bet2_star = bet2/sqrt(gm1 - (gm1-1.0)*y);
    Real bet3_star = bet3/sqrt(gm1 - (gm1-1.0)*y);
    Real bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
    Real vbet = v2*bet2_star + v3*bet3_star;

    // Compute alpha(s) (eq. A16)
    Real alpha_f=1.0;
    Real alpha_s=0.0;
    if ((twid_asq - cssq) <= 0.0) {
      alpha_f = 0.0;
      alpha_s = 1.0;
    } else if ((cfsq - twid_asq) <= 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
    } else if ((cfsq-cssq) != 0.0){
      alpha_f = sqrt((twid_asq - cssq)/(cfsq - cssq));
      alpha_s = sqrt((cfsq - twid_asq)/(cfsq - cssq));
    }

    // Compute Q(s) and A(s) (eq. A14-15), etc.
    Real sqrtd = sqrt(d);
    Real isqrtd = 1.0/sqrtd;
    Real s = SIGN(b1);
    Real twid_a = sqrt(twid_asq);
    Real qf = cf*alpha_f*s;
    Real qs = cs*alpha_s*s;
    Real af_prime = twid_a*alpha_f*isqrtd;
    Real as_prime = twid_a*alpha_s*isqrtd;
    Real afpbb = af_prime*bt_star*bet_starsq;
    Real aspbb = as_prime*bt_star*bet_starsq;

    // Compute eigenvalues (eq. B17)
    Real vax = sqrt(vaxsq);
    eigenvalues[0] = v1 - cf;
    eigenvalues[1] = v1 - vax;
    eigenvalues[2] = v1 - cs;
    eigenvalues[3] = v1;
    eigenvalues[4] = v1 + cs;
    eigenvalues[5] = v1 + vax;
    eigenvalues[6] = v1 + cf;
  
    // Right-eigenvectors, stored as COLUMNS (eq. B21) */
    right_eigenmatrix[0][0] = alpha_f;
    right_eigenmatrix[0][1] = 0.0; 
    right_eigenmatrix[0][2] = alpha_s;
    right_eigenmatrix[0][3] = 1.0;
    right_eigenmatrix[0][4] = alpha_s;
    right_eigenmatrix[0][5] = 0.0;
    right_eigenmatrix[0][6] = alpha_f;

    right_eigenmatrix[1][0] = alpha_f*eigenvalues[0];
    right_eigenmatrix[1][1] = 0.0;
    right_eigenmatrix[1][2] = alpha_s*eigenvalues[2];
    right_eigenmatrix[1][3] = v1;
    right_eigenmatrix[1][4] = alpha_s*eigenvalues[4];
    right_eigenmatrix[1][5] = 0.0;
    right_eigenmatrix[1][6] = alpha_f*eigenvalues[6];

    Real qa = alpha_f*v2;
    Real qb = alpha_s*v2;
    Real qc = qs*bet2_star;
    Real qd = qf*bet2_star;
    right_eigenmatrix[2][0] = qa + qc;
    right_eigenmatrix[2][1] = -bet3;
    right_eigenmatrix[2][2] = qb - qd;
    right_eigenmatrix[2][3] = v2;
    right_eigenmatrix[2][4] = qb + qd;
    right_eigenmatrix[2][5] = bet3;
    right_eigenmatrix[2][6] = qa - qc;

    qa = alpha_f*v3;
    qb = alpha_s*v3;
    qc = qs*bet3_star;
    qd = qf*bet3_star;
    right_eigenmatrix[3][0] = qa + qc;
    right_eigenmatrix[3][1] = bet2;
    right_eigenmatrix[3][2] = qb - qd;
    right_eigenmatrix[3][3] = v3;
    right_eigenmatrix[3][4] = qb + qd;
    right_eigenmatrix[3][5] = -bet2;
    right_eigenmatrix[3][6] = qa - qc;

    right_eigenmatrix[4][0] = alpha_f*(hp - v1*cf) + qs*vbet + aspbb;
    right_eigenmatrix[4][1] = -(v2*bet3 - v3*bet2);
    right_eigenmatrix[4][2] = alpha_s*(hp - v1*cs) - qf*vbet - afpbb;
    right_eigenmatrix[4][3] = 0.5*vsq + (gm1-1.0)*x/gm1;
    right_eigenmatrix[4][4] = alpha_s*(hp + v1*cs) + qf*vbet - afpbb;
    right_eigenmatrix[4][5] = -right_eigenmatrix[4][1];
    right_eigenmatrix[4][6] = alpha_f*(hp + v1*cf) - qs*vbet + aspbb;
    
    right_eigenmatrix[5][0] = as_prime*bet2_star;
    right_eigenmatrix[5][1] = -bet3*s*isqrtd;
    right_eigenmatrix[5][2] = -af_prime*bet2_star;
    right_eigenmatrix[5][3] = 0.0; 
    right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
    right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
    right_eigenmatrix[5][6] = right_eigenmatrix[5][0];

    right_eigenmatrix[6][0] = as_prime*bet3_star;
    right_eigenmatrix[6][1] = bet2*s*isqrtd;
    right_eigenmatrix[6][2] = -af_prime*bet3_star;
    right_eigenmatrix[6][3] = 0.0;
    right_eigenmatrix[6][4] = right_eigenmatrix[6][2];
    right_eigenmatrix[6][5] = right_eigenmatrix[6][1];
    right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

    // Left-eigenvectors, stored as ROWS (eq. B29)
    // Normalize by 1/2a^{2}: quantities denoted by \hat{f}
    Real norm = 0.5/twid_asq;
    Real cff = norm*alpha_f*cf;
    Real css = norm*alpha_s*cs;
    qf *= norm;
    qs *= norm;
    Real af = norm*af_prime*d;
    Real as = norm*as_prime*d;
    Real afpb = norm*af_prime*bt_star;
    Real aspb = norm*as_prime*bt_star;
  
    // Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f}
    norm *= gm1;
    alpha_f *= norm;
    alpha_s *= norm;
    Real q2_star = bet2_star/bet_starsq;
    Real q3_star = bet3_star/bet_starsq;
    Real vqstr = (v2*q2_star + v3*q3_star);
    norm *= 2.0;
  
    left_eigenmatrix[0][0] = alpha_f*(vsq-hp) + cff*(cf+v1) - qs*vqstr - aspb;
    left_eigenmatrix[0][1] = -alpha_f*v1 - cff;
    left_eigenmatrix[0][2] = -alpha_f*v2 + qs*q2_star;
    left_eigenmatrix[0][3] = -alpha_f*v3 + qs*q3_star;
    left_eigenmatrix[0][4] = alpha_f;
    left_eigenmatrix[0][5] = as*q2_star - alpha_f*b2;
    left_eigenmatrix[0][6] = as*q3_star - alpha_f*b3;

    left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
    left_eigenmatrix[1][1] = 0.0;
    left_eigenmatrix[1][2] = -0.5*bet3;
    left_eigenmatrix[1][3] = 0.5*bet2;
    left_eigenmatrix[1][4] = 0.0; 
    left_eigenmatrix[1][5] = -0.5*sqrtd*bet3*s;
    left_eigenmatrix[1][6] = 0.5*sqrtd*bet2*s;
  
    left_eigenmatrix[2][0] = alpha_s*(vsq-hp) + css*(cs+v1) + qf*vqstr + afpb;
    left_eigenmatrix[2][1] = -alpha_s*v1 - css;
    left_eigenmatrix[2][2] = -alpha_s*v2 - qf*q2_star;
    left_eigenmatrix[2][3] = -alpha_s*v3 - qf*q3_star;
    left_eigenmatrix[2][4] = alpha_s;
    left_eigenmatrix[2][5] = -af*q2_star - alpha_s*b2;
    left_eigenmatrix[2][6] = -af*q3_star - alpha_s*b3;

    left_eigenmatrix[3][0] = 1.0 - norm*(0.5*vsq - (gm1-1.0)*x/gm1); 
    left_eigenmatrix[3][1] = norm*v1;
    left_eigenmatrix[3][2] = norm*v2;
    left_eigenmatrix[3][3] = norm*v3;
    left_eigenmatrix[3][4] = -norm;
    left_eigenmatrix[3][5] = norm*b2;
    left_eigenmatrix[3][6] = norm*b3;

    left_eigenmatrix[4][0] = alpha_s*(vsq-hp) + css*(cs-v1) - qf*vqstr + afpb;
    left_eigenmatrix[4][1] = -alpha_s*v1 + css;
    left_eigenmatrix[4][2] = -alpha_s*v2 + qf*q2_star;
    left_eigenmatrix[4][3] = -alpha_s*v3 + qf*q3_star;
    left_eigenmatrix[4][4] = alpha_s;
    left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
    left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

    left_eigenmatrix[5][0] = -left_eigenmatrix[1][0];
    left_eigenmatrix[5][1] = 0.0;
    left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
    left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
    left_eigenmatrix[5][4] = 0.0;
    left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
    left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

    left_eigenmatrix[6][0] = alpha_f*(vsq-hp) + cff*(cf-v1) + qs*vqstr - aspb;
    left_eigenmatrix[6][1] = -alpha_f*v1 + cff;
    left_eigenmatrix[6][2] = -alpha_f*v2 - qs*q2_star;
    left_eigenmatrix[6][3] = -alpha_f*v3 - qs*q3_star;
    left_eigenmatrix[6][4] = alpha_f;
    left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
    left_eigenmatrix[6][6] = left_eigenmatrix[0][6];

// Isothermal MHD

  } else {
    Real di = 1.0/d;
    Real btsq = b2*b2 + b3*b3;
    Real bt_starsq = btsq*y;
    Real vaxsq = b1*b1*di;
    Real twid_csq = (iso_cs*iso_cs) + x;

    // Compute fast- and slow-magnetosonic speeds (eq. B39)
    Real ct2 = bt_starsq*di;
    Real tsum = vaxsq + ct2 + twid_csq;
    Real tdif = vaxsq + ct2 - twid_csq;
    Real cf2_cs2 = sqrt(tdif*tdif + 4.0*twid_csq*ct2);
  
    Real cfsq = 0.5*(tsum + cf2_cs2);
    Real cf = sqrt(cfsq);
   
    Real cssq = twid_csq*vaxsq/cfsq;
    Real cs = sqrt(cssq);
  
    // Compute beta's (eqs. A17, B28, B40)
    Real bt = sqrt(btsq);
    Real bt_star = sqrt(bt_starsq);
    Real bet2 = 1.0;
    Real bet3 = 0.0;
    if (bt != 0.0) {
      bet2 = b2/bt;
      bet3 = b3/bt;
    }
    Real bet2_star = bet2/sqrt(y);
    Real bet3_star = bet3/sqrt(y);
    Real bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;

    // Compute alpha's (eq. A16)
    Real alpha_f = 1.0;
    Real alpha_s = 0.0;
    if ((twid_csq - cssq) <= 0.0) {
      alpha_f = 0.0;
      alpha_s = 1.0;
    } else if ((cfsq - twid_csq) <= 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
    } else if ((cfsq-cssq) != 0.0) {
      alpha_f = sqrt((twid_csq - cssq)/(cfsq - cssq));
      alpha_s = sqrt((cfsq - twid_csq)/(cfsq - cssq));
    }

    // Compute Q's (eq. A14-15), etc.
    Real sqrtd = sqrt(d);
    Real s = SIGN(b1);
    Real twid_c = sqrt(twid_csq);
    Real qf = cf*alpha_f*s;
    Real qs = cs*alpha_s*s;
    Real af_prime = twid_c*alpha_f/sqrtd;
    Real as_prime = twid_c*alpha_s/sqrtd;

    // Compute eigenvalues (eq. B38)
    Real vax  = sqrt(vaxsq);
    eigenvalues[0] = v1 - cf;
    eigenvalues[1] = v1 - vax;
    eigenvalues[2] = v1 - cs;
    eigenvalues[3] = v1 + cs;
    eigenvalues[4] = v1 + vax;
    eigenvalues[5] = v1 + cf;

    // Right-eigenvectors, stored as COLUMNS (eq. B21)
    right_eigenmatrix[0][0] = alpha_f;
    right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
    right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
    right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
    right_eigenmatrix[4][0] = as_prime*bet2_star;
    right_eigenmatrix[5][0] = as_prime*bet3_star;

    right_eigenmatrix[0][1] = 0.0;
    right_eigenmatrix[1][1] = 0.0;
    right_eigenmatrix[2][1] = -bet3;
    right_eigenmatrix[3][1] = bet2;
    right_eigenmatrix[4][1] = -bet3*s/sqrtd;
    right_eigenmatrix[5][1] = bet2*s/sqrtd;

    right_eigenmatrix[0][2] = alpha_s;
    right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
    right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
    right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
    right_eigenmatrix[4][2] = -af_prime*bet2_star;
    right_eigenmatrix[5][2] = -af_prime*bet3_star;

    right_eigenmatrix[0][3] = alpha_s;
    right_eigenmatrix[1][3] = alpha_s*(v1 + cs);
    right_eigenmatrix[2][3] = alpha_s*v2 + qf*bet2_star;
    right_eigenmatrix[3][3] = alpha_s*v3 + qf*bet3_star;
    right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
    right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

    right_eigenmatrix[0][4] = 0.0;
    right_eigenmatrix[1][4] = 0.0; 
    right_eigenmatrix[2][4] = bet3;
    right_eigenmatrix[3][4] = -bet2;
    right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
    right_eigenmatrix[5][4] = right_eigenmatrix[5][1];
  
    right_eigenmatrix[0][5] = alpha_f;
    right_eigenmatrix[1][5] = alpha_f*(v1 + cf);
    right_eigenmatrix[2][5] = alpha_f*v2 - qs*bet2_star;
    right_eigenmatrix[3][5] = alpha_f*v3 - qs*bet3_star;
    right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
    right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

    // Left-eigenvectors, stored as ROWS (eq. B41)
    // Normalize by 1/2a^{2}: quantities denoted by \hat{f}
    Real norm = 0.5/twid_csq;
    Real cff = norm*alpha_f*cf;
    Real css = norm*alpha_s*cs;
    qf *= norm;
    qs *= norm;
    Real af = norm*af_prime*d;
    Real as = norm*as_prime*d;
    Real afpb = norm*af_prime*bt_star;
    Real aspb = norm*as_prime*bt_star;

    Real q2_star = bet2_star/bet_starsq;
    Real q3_star = bet3_star/bet_starsq;
    Real vqstr = (v2*q2_star + v3*q3_star);
  
    left_eigenmatrix[0][0] = cff*(cf+v1) - qs*vqstr - aspb;
    left_eigenmatrix[0][1] = -cff;
    left_eigenmatrix[0][2] = qs*q2_star;
    left_eigenmatrix[0][3] = qs*q3_star;
    left_eigenmatrix[0][4] = as*q2_star;
    left_eigenmatrix[0][5] = as*q3_star;
  
    left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
    left_eigenmatrix[1][1] = 0.0;
    left_eigenmatrix[1][2] = -0.5*bet3;
    left_eigenmatrix[1][3] = 0.5*bet2;
    left_eigenmatrix[1][4] = -0.5*sqrtd*bet3*s;
    left_eigenmatrix[1][5] = 0.5*sqrtd*bet2*s;
  
    left_eigenmatrix[2][0] = css*(cs+v1) + qf*vqstr + afpb;
    left_eigenmatrix[2][1] = -css;
    left_eigenmatrix[2][2] = -qf*q2_star;
    left_eigenmatrix[2][3] = -qf*q3_star;
    left_eigenmatrix[2][4] = -af*q2_star;
    left_eigenmatrix[2][5] = -af*q3_star;
  
    left_eigenmatrix[3][0] = css*(cs-v1) - qf*vqstr + afpb;
    left_eigenmatrix[3][1] = css;
    left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
    left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
    left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
    left_eigenmatrix[3][5] = left_eigenmatrix[2][5];
  
    left_eigenmatrix[4][0] = -left_eigenmatrix[1][0];
    left_eigenmatrix[4][1] = 0.0;
    left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
    left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
    left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
    left_eigenmatrix[4][5] = left_eigenmatrix[1][5];
  
    left_eigenmatrix[5][0] = cff*(cf-v1) + qs*vqstr - aspb;
    left_eigenmatrix[5][1] = cff;
    left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
    left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
    left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
    left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
  }
}
