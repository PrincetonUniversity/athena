//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"

//======================================================================================
/*! \file  roe.cpp
 *  \brief Roe's linearized Riemann solver.
 *
 * Computes 1D fluxes using Roe's linearization.  When Roe's method fails because of
 * negative density or pressure in the intermediate states, the fluxes are computed with
 * the HLLE solver instead.
 *
 * REFERENCES:
 * - P. Roe, "Approximate Riemann solvers, parameter vectors, and difference schemes",
 *   JCP, 43, 357 (1981).
 *====================================================================================*/

/* maximum wavespeed used by H-correction, value passed from integrator */
Real etah=0.0;

/*
void RiemannSolver::Roe(const int il, const int iu,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  Real sqrtdl,sqrtdr,isdlpdr,droe,v1roe,v2roe,v3roe,pbl=0.0,pbr=0.0;
  Real hroe;
  Real coeff[NWAVE];
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real dU[NWAVE],a[NWAVE];
  Real *pUl, *pUr, *pFl, *pFr, *pF;
  Cons1DS Fl,Fr;
  int n,m,hlle_flag;

  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }

//--- Step 2. ------------------------------------------------------------------
// Compute Roe-averaged data from left- and right-states
///

  sqrtdl = sqrt((double)Wl.d);
  sqrtdr = sqrt((double)Wr.d);
  isdlpdr = 1.0/(sqrtdl + sqrtdr);

  droe  = sqrtdl*sqrtdr;
  v1roe = (sqrtdl*Wl.Vx + sqrtdr*Wr.Vx)*isdlpdr;
  v2roe = (sqrtdl*Wl.Vy + sqrtdr*Wr.Vy)*isdlpdr;
  v3roe = (sqrtdl*Wl.Vz + sqrtdr*Wr.Vz)*isdlpdr;

//
// Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
// rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
///

  hroe  = ((Ul.E + Wl.P + pbl)/sqrtdl + (Ur.E + Wr.P + pbr)/sqrtdr)*isdlpdr;

//--- Step 3. ------------------------------------------------------------------
// Compute eigenvalues and eigenmatrices using Roe-averaged values
///

  esys_roe_adb_hyd(v1roe,v2roe,v3roe,hroe,ev,rem,lem);


//--- Step 4. ------------------------------------------------------------------
// Compute L/R fluxes 
///

  Fl.d  = Ul.Mx;
  Fr.d  = Ur.Mx;

  Fl.Mx = Ul.Mx*Wl.Vx;
  Fr.Mx = Ur.Mx*Wr.Vx;

  Fl.My = Ul.Mx*Wl.Vy;
  Fr.My = Ur.Mx*Wr.Vy;

  Fl.Mz = Ul.Mx*Wl.Vz;
  Fr.Mz = Ur.Mx*Wr.Vz;

  Fl.Mx += Wl.P;
  Fr.Mx += Wr.P;

  Fl.E  = (Ul.E + Wl.P)*Wl.Vx;
  Fr.E  = (Ur.E + Wr.P)*Wr.Vx;

//--- Step 5. ------------------------------------------------------------------
// Return upwind flux if flow is supersonic
///

  if(ev[0] >= 0.0){
    *pFlux = Fl;
    return;
  }

  if(ev[NWAVE-1] <= 0.0){
    *pFlux = Fr;
    return;
  }

//--- Step 6. ------------------------------------------------------------------
// Compute projection of dU onto L eigenvectors ("vector A")
///

  pUr = (Real *)&(Ur);
  pUl = (Real *)&(Ul);

  for (n=0; n<NWAVE; n++) dU[n] = pUr[n] - pUl[n];
  for (n=0; n<NWAVE; n++) {
    a[n] = 0.0;
    for (m=0; m<NWAVE; m++) a[n] += lem[n][m]*dU[m];
  }

//--- Step 7. ------------------------------------------------------------------
// Check that the density and pressure in the intermediate states are positive.
// If not, set hlle_flag=1 if d_inter<0; hlle_flag=2 if p_inter<0, get HLLE
// fluxes, and return
///

  hlle_flag = 0;
#ifdef TEST_INTERMEDIATE_STATES

  for (n=0; n<NWAVE; n++) u_inter[n] = pUl[n];
  for (n=0; n<NWAVE-1; n++) {
    for (m=0; m<NWAVE; m++) u_inter[m] += a[n]*rem[m][n];
    if(ev[n+1] > ev[n]) {
      if (u_inter[0] <= 0.0) {
	hlle_flag=1;
	break;
      }
      p_inter = u_inter[4] - 0.5*
	(SQR(u_inter[1])+SQR(u_inter[2])+SQR(u_inter[3]))/u_inter[0];
      if (p_inter < 0.0) {
	hlle_flag=2;
	break;
      }
    }
  }

  if (hlle_flag != 0) {
    flux_hlle(Ul,Ur,Wl,Wr,Bxi,pFlux);
    return;
  }

#endif // TEST_INTERMEDIATE_STATES

//--- Step 8. ------------------------------------------------------------------
// Compute Roe flux

  pFl = (Real *)&(Fl);
  pFr = (Real *)&(Fr);
  pF  = (Real *)(pFlux); 

  for (m=0; m<NWAVE; m++) {
    coeff[m] = 0.5*MAX(fabs(ev[m]),etah)*a[m];
  }
  for (n=0; n<NWAVE; n++) {
    pF[n] = 0.5*(pFl[n] + pFr[n]);
    for (m=0; m<NWAVE; m++) {
      pF[n] -= coeff[m]*rem[n][m];
    }
  }

  return;
}

*/

/*============================================================================*/
/*! \file esystem_roe.c
 *  \brief Functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors of "Roe's matrix A" for the linearized system in the 
 * CONSERVED variables.
 *
 * PURPOSE: Functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors of "Roe's matrix A" for the linearized system in the 
 * CONSERVED variables, i.e. U,t = AU,x, where U=(d,d*vx,d*vy,d*vz,[E],[By,Bz]).
 * The eigenvalues are returned through the argument list as a vector of length
 * NWAVE.  The eigenvectors are returned as matrices of size (NWAVE)x(NWAVE),
 * with right-eigenvectors stored as COLUMNS (so R_i = right_eigenmatrix[*][i]),
 * and left-eigenvectors stored as ROWS (so L_i = left_eigenmatrix[i][*]).
 *
 * To improve performance components of the eigenvectors which are zero
 * are not set here (eigenmatrices must be initialized to zero in calling
 * routine).  However, for completeness, statements which set these values
 * are included but are commented out.
 *
 * The "Roe-averaging" of the L/R states must be performed in the calling funct
 *
 * REFERENCES:
 * - P. Cargo & G. Gallice, "Roe matrices for ideal MHD and systematic
 *   construction of Roe matrices for systems of conservation laws",
 *   JCP, 136, 446 (1997)
 *
 * - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new
 *   code for astrophysical MHD", ApJS, (2008), Appendix B
 *   Equation numbers refer to this paper.
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \fn void esys_roe_adb_hyd(const Real v1, const Real v2, const Real v3, 
 *			      const Real h, Real right_eigenmatrix[][5], 
 *			      Real left_eigenmatrix[][5])
 *  \brief ADIABATIC HYDRO
 *
 * - Input: v1,v2,v3,h = Roe averaged velocities and enthalpy
 * - Output: eigenvalues[5], right_eigenmatrix[5,5], left_eigenmatrix[5,5];
 */

/*
void esys_roe_adb_hyd(const Real v1, const Real v2, const Real v3, const Real h,
  Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5])
{
  Real vsq,asq,a,na,qa;
  vsq = v1*v1 + v2*v2 + v3*v3;
  asq = Gamma_1*MAX((h-0.5*vsq), TINY_NUMBER);
  a = sqrt(asq);

// Compute eigenvalues (eq. B2)

  eigenvalues[0] = v1 - a;
  eigenvalues[1] = v1;
  eigenvalues[2] = v1;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + a;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

// Right-eigenvectors, stored as COLUMNS (eq. B3)

  right_eigenmatrix[0][0] = 1.0;
  right_eigenmatrix[1][0] = v1 - a;
  right_eigenmatrix[2][0] = v2;
  right_eigenmatrix[3][0] = v3;
  right_eigenmatrix[4][0] = h - v1*a;

//right_eigenmatrix[0][1] = 0.0;
//right_eigenmatrix[1][1] = 0.0;
  right_eigenmatrix[2][1] = 1.0;
//right_eigenmatrix[3][1] = 0.0;
  right_eigenmatrix[4][1] = v2;

//right_eigenmatrix[0][2] = 0.0;
//right_eigenmatrix[1][2] = 0.0;
//right_eigenmatrix[2][2] = 0.0; 
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

  na = 0.5/asq;
  left_eigenmatrix[0][0] = na*(0.5*Gamma_1*vsq + v1*a);
  left_eigenmatrix[0][1] = -na*(Gamma_1*v1 + a);
  left_eigenmatrix[0][2] = -na*Gamma_1*v2;
  left_eigenmatrix[0][3] = -na*Gamma_1*v3;
  left_eigenmatrix[0][4] = na*Gamma_1;

  left_eigenmatrix[1][0] = -v2;
//left_eigenmatrix[1][1] = 0.0;
  left_eigenmatrix[1][2] = 1.0;
//left_eigenmatrix[1][3] = 0.0;
//left_eigenmatrix[1][4] = 0.0;

  left_eigenmatrix[2][0] = -v3;
//left_eigenmatrix[2][1] = 0.0;
//left_eigenmatrix[2][2] = 0.0;
  left_eigenmatrix[2][3] = 1.0;
//left_eigenmatrix[2][4] = 0.0; 

  qa = Gamma_1/asq;
  left_eigenmatrix[3][0] = 1.0 - na*Gamma_1*vsq;
  left_eigenmatrix[3][1] = qa*v1;
  left_eigenmatrix[3][2] = qa*v2;
  left_eigenmatrix[3][3] = qa*v3;
  left_eigenmatrix[3][4] = -qa;

  left_eigenmatrix[4][0] = na*(0.5*Gamma_1*vsq - v1*a);
  left_eigenmatrix[4][1] = -na*(Gamma_1*v1 - a);
  left_eigenmatrix[4][2] = left_eigenmatrix[0][2];
  left_eigenmatrix[4][3] = left_eigenmatrix[0][3];
  left_eigenmatrix[4][4] = left_eigenmatrix[0][4];
}
*/
