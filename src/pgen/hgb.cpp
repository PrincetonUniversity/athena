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
//! \file hgb.cpp
/*! \file hgb.c
 *  \brief Problem generator for 3D shearing sheet.
 *
 * PURPOSE:  Problem generator for 3D shearing sheet.  Based on the initial
 *   conditions described in "Local Three-dimensional Magnetohydrodynamic
 *   Simulations of Accretion Disks" by Hawley, Gammie & Balbus, or HGB.
 *
 * Several different field configurations and perturbations are possible:
 *
 *- ifield = 0 - uses field set by choice of ipert flag
 *- ifield = 1 - Bz=B0sin(kx*x1) field with zero-net-flux [default] (kx input)
 *- ifield = 2 - uniform Bz
 *- ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *- ifield = 4 - B=(0,B0/sqrt(2),B0/sqrt(2))= net toroidal+vertical field
 *- ifield = 5 - uniform By
 *
 *- ipert = 1 - random perturbations to P and V [default, used by HGB]
 *- ipert = 2 - uniform Vx=amp (epicyclic wave test)
 *- ipert = 3 - J&G vortical shwave (hydro test)
 *- ipert = 4 - nonlinear density wave test of Fromang & Papaloizou
 *- ipert = 5 - 2nd MHD shwave test of JGG (2008) -- their figure 9
 *- ipert = 6 - 3rd MHD shwave test of JGG (2008) -- their figure 11
 *- ipert = 7 - nonlinear shearing wave test of Heinemann & Papaloizou (2008)
 *
 * To run simulations of stratified disks (including vertical gravity), use the
 * strat.c problem generator.
 *
 * Code must be configured using --enable-shearing-box
 *
 * REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).
 *            Johnson, Guan, & Gammie, ApJSupp, (2008)			      */
/*============================================================================*/


// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstdlib>    // exit()
#include <cfloat>     // DBL_EPSILON

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief field loop advection problem generator for 2D/3D problems.
//======================================================================================

Real Lx,Ly,Lz; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);
/*
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real expr_Jsq(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif
#ifdef MHD
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k);
static Real hst_dEw2(const GridS *pG, const int i, const int j, const int k);
static Real hst_dBy(const GridS *pG, const int i, const int j, const int k);
#endif
*/
/*----------------------------------------------------------------------------*/

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  FILE *fp;
  Real xFP[160],dFP[160],vxFP[160],vyFP[160];
  static int frst=1; // flag so new history variables enrolled only once

  if (pmy_mesh->mesh_size.nx2 == 1){
	  std::cout << "[hgb.cpp]: HGB only works on a 2D or 3D grid" << std::endl;
  }

// Read problem parameters
  Real Omega_0 = pin->GetOrAddReal("problem","Omega0",1.0e-3);
  Real qshear  = pin->GetOrAddReal("problem","qshear",1.5);
  Real amp = pin->GetReal("problem","amp");
  int ipert = pin->GetOrAddInteger("problem","ipert", 1);

  Real beta, dir_sgn;
  int ifield, Bdir;
  if (MAGNETIC_FIELDS_ENABLED) {
    beta = pin->GetReal("problem","beta");
    ifield = pin->GetOrAddInteger("problem","ifield", 1);
    // For net-flux calc, provide the direction of the B field
    Bdir = pin->GetOrAddInteger("problem","Bdir",1);
    if (Bdir > 0)
      dir_sgn = 1.0;
    else
      dir_sgn = -1.0;
  }

// Compute pressure based on the EOS.
  Real den = 1.0, pres =1.0, gamma=1.0, iso_cs=1.0;
  if (NON_BAROTROPIC_EOS) {
    Real gamma = peos->GetGamma();
    Real pres = pin->GetReal("problem","pres");
  } else {
	Real iso_cs =peos->GetIsoSoundSpeed();
    Real pres = den*SQR(iso_cs);
  }
// Compute field strength based on beta.
  Real B0  = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
    B0 = sqrt((double)(2.0*pres/beta));

// Ensure a different initial random seed for each meshblock.
  long int iseed = -1 - gid;

// Initialize boxsize
  Lx = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  Ly = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

// initialize wavenumbers
  int nwx = pin->GetOrAddInteger("problem","nwx",1);
  int nwy = pin->GetOrAddInteger("problem","nwy",1);
  int nwz = pin->GetOrAddInteger("problem","nwz",1);
  Real kx = (2.0*PI/Lx)*((double)nwx);// nxw=-ve for leading wave
  Real ky = (2.0*PI/Ly)*((double)nwy);
  Real kz = (2.0*PI/Lz)*((double)nwz);

// For PF density wave test, read data from file: not implemented yet.


/* Rescale amp to sound speed for ipert 2,3 */
  if (NON_BAROTROPIC_EOS) {
    if (ipert == 2 || ipert == 3)
	  amp *= sqrt(gamma*pres/den);
  } else {
    if (ipert == 2 || ipert == 3)
	  amp *= iso_cs;
  }

  Real x1,x2,x3,xmin,xmax;
  Real x1f,x2f,x3f;
  Real rd,rp,rvx,rvy,rvz,rbx,rby,rbz;
  Real rval;
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      x1 = pcoord->x1v(i);
	  x2 = pcoord->x2v(j);
	  x3 = pcoord->x3v(k);
      x1f = pcoord->x1f(i);
	  x2f = pcoord->x2f(j);
	  x3f = pcoord->x3f(k);

//Initialize perturbations
// ipert = 1 - random perturbations to P and V [default, used by HGB]
// ipert = 2 - uniform Vx=amp (epicyclic wave test)
// ipert = 3 - vortical shwave (hydro test)
// ipert = 4 - Fromang & Papaloizou nonlinear density wave (hydro test)
// ipert = 5 & 6 - JGG MHD shwave tests
// ipert = 7 - Heinemann & Papaloizou (2008) nonlinear shwave (hydro test)
      if (ipert == 1) {
        rval = amp*(ran2(&iseed) - 0.5);
        if (NON_BAROTROPIC_EOS) {
          rp = pres*(1.0 + 2.0*rval);
          rd = den;
        } else {
          rd = den*(1.0 + 2.0*rval);
		}
        // Follow HGB: the perturbations to V/Cs are
		// (1/5)amp/sqrt(gamma)
        rval = amp*(ran2(&iseed) - 0.5);
        rvx = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvy = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvz = 0.4*rval*sqrt(pres/den);
      }
      if (ipert == 2) {
        rp = pres;
        rd = den;
        rvx = amp;
        rvy = 0.0;
        rvz = 0.0;
      }
      if (ipert == 3) {
        rp = pres;
        rd = den;
        rvx = amp*sin((double)(kx*x1 + ky*x2));
        rvy = -amp*(kx/ky)*sin((double)(kx*x1 + ky*x2));
        rvz = 0.0;
      }
      if (ipert == 4) {
		std::cout << "[hgb.cpp]: ipert=4 not implemented yet!" << std::endl;
		exit(0);
      }
      // Note: ICs in JGG for this test are incorrect.
      if (ipert == 5) {
        ifield = 0;
        rd = den + 8.9525e-10*cos((double)(kx*x1 + ky*x2 + kz*x3 - PI/4.));
        rvx = 8.16589e-8*cos((double)(kx*x1 + ky*x2 + kz*x3 + PI/4.));
        rvy = 8.70641e-8*cos((double)(kx*x1 + ky*x2 + kz*x3 + PI/4.));
        rvz = 0.762537e-8*cos((double)(kx*x1 + ky*x2 + kz*x3 + PI/4.));
        rbx = -1.08076e-7;
        rbx *= cos((double)(kx*(x1-0.5*pcoord->dx1f(i)) + ky*x2 + kz*x3 - PI/4.));
        rby = 1.04172e-7;
        rby *= cos((double)(kx*x1 + ky*(x2-0.5*pcoord->dx2f(j)) + kz*x3 - PI/4.));
        rbz = -0.320324e-7;
        rbz *= cos((double)(kx*x1 + ky*x2 + kz*(x3-0.5*pcoord->dx3f(k)) - PI/4.));;
        rbz += (sqrt(15.0)/16.0)*(Omega_0/kz);
      }
      if (ipert == 6) {
        ifield = 0;
        rd = den + 5.48082e-6*cos((double)(kx*x1 + ky*x2 + kz*x3));
        rvx = -4.5856e-6*cos((double)(kx*x1 + ky*x2 + kz*x3));
        rvy = 2.29279e-6*cos((double)(kx*x1 + ky*x2 + kz*x3));
        rvz = 2.29279e-6*cos((double)(kx*x1 + ky*x2 + kz*x3));
        rbx = 5.48082e-7;
        rbx *= cos((double)(kx*x1f + ky*x2 + kz*x3));
        rbx += (0.1);
        rby = 1.0962e-6;
        rby *= cos((double)(kx*x1 + ky*x2f + kz*x3));
        rby += (0.2);
        rbz = 0.0;
      }
      if (ipert == 7) {
        if (!NON_BAROTROPIC_EOS) {
          Real kappa2 = 2.0*(2.0 - qshear)*Omega_0*Omega_0;
          Real aa = (kx*kx + ky*ky)*SQR(iso_cs) + kappa2;
          Real bb = 2.0*qshear*Omega_0*ky*iso_cs;
          Real denom = aa*aa + bb*bb;
          Real rd_hat =         (ky*iso_cs*bb -2.0*Omega_0*aa)*amp/denom;
          Real px_hat = -iso_cs*(ky*iso_cs*aa +2.0*Omega_0*bb)*amp/denom;
          Real py_hat = (amp + ky*px_hat + (2.0-qshear)*Omega_0*rd_hat)/kx;
          rd  = 1.0 + rd_hat*cos((double)(kx*x1 + ky*x2));
          rvx = px_hat*sin((double)(kx*x1 + ky*x2))/rd;
          rvy = py_hat*sin((double)(kx*x1 + ky*x2))/rd;
		}
        rvz = 0.0;
      }

// Initialize (d, M, P)
// for_the_future: if FARGO do not initialize the bg shear
      phydro->u(IDN,k,j,i) = rd;
      phydro->u(IM1,k,j,i) = rd*rvx;
      phydro->u(IM2,k,j,i) = rd*rvy;
      phydro->u(IM2,k,j,i) -= rd*(qshear*Omega_0*x1);
      phydro->u(IM3,k,j,i) = rd*rvz;
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = rp/(gamma-1.0)
          + 0.5*(SQR(phydro->u(IM1,k,j,i))
			   + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/rd;
      } // Hydro

// Initialize b.  For 3D shearing box B1=Bx, B2=By, B3=Bz
// ifield = 0 - used with ipert=5 or 6
// ifield = 1 - Bz=B0sin(x1) field with zero-net-flux[default]
// ifield = 2 - uniform Bz
// ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))=zero-net flux w helicity
// ifield = 4 - B=(0,B0/sqrt(2),B0/sqrt(2))= net toroidal+vertical field
      if (MAGNETIC_FIELDS_ENABLED) {
        if (ifield == 0) {
          pfield->b.x1f(k,j,i) = rbx;
          pfield->b.x2f(k,j,i) = rby;
          pfield->b.x3f(k,j,i) = rbz;
          if (i==ie) {
            x1f = pcoord->x1f(ie+1);
            rbx = 5.48082e-7;
            rbx *= cos((double)(kx*x1f + ky*x2 + kz*x3));
            rbx += (0.1);
			pfield->b.x1f(k,j,ie+1) =  rbx;
		  }
          if (j==je) {
            x2f = pcoord->x2f(je+1);
            rby = 1.0962e-6;
            rby *= cos((double)(kx*x1 + ky*x2f + kz*x3));
            rby += (0.2);
			pfield->b.x2f(k,je+1,i) =  rby;
		  }
          if (k==ke) {
            x3f = pcoord->x3f(ke+1);
            rbz = 0.0;
			pfield->b.x3f(ke+1,j,i) =  rbz;
		  }
          //pfield->b.x1f(k,j,i) = 0.0;
          //pfield->b.x2f(k,j,i) = 0.0;
          //pfield->b.x3f(k,j,i) = 0.0;
          //if (i==ie) pfield->b.x1f(k,j,ie+1) =  0.0;
          //if (j==je) pfield->b.x2f(k,je+1,i) =  0.0;
          //if (k==ke) pfield->b.x3f(ke+1,j,i) =  0.0;
        }
        if (ifield == 1) {
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = 0.0;
          pfield->b.x3f(k,j,i) = B0*(sin((double)kx*x1));
          if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
          if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*(sin((double)kx*x1));
        }
        if (ifield == 2) {
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = 0.0;
          pfield->b.x3f(k,j,i) = B0*dir_sgn;
          if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
          if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*dir_sgn;
        }
        if (ifield == 3) {
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = B0*(cos((double)kx*x1));
          pfield->b.x3f(k,j,i) = B0*(sin((double)kx*x1));
          if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(k,je+1,i) = B0*(cos((double)kx*x1));
          if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*(sin((double)kx*x1));
        }
        if (ifield == 4) {
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = B0/sqrt(2);
          pfield->b.x3f(k,j,i) = B0/sqrt(2);
          if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(k,je+1,i) = B0/sqrt(2);
          if (k==ke) pfield->b.x3f(ke+1,j,i) = B0/sqrt(2);
        }
        if (ifield == 5) {
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = B0;
          pfield->b.x3f(k,j,i) = 0.0;
          if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(k,je+1,i) = B0;
          if (k==ke) pfield->b.x3f(ke+1,j,i) = 0.0;
        }
      } // MHD
    }
  }}

// calc the cc-B seems unnecessary now
//#ifdef MHD
//  for (k=ks; k<=ke; k++) {
//    for (j=js; j<=je; j++) {
//      for (i=is; i<=ie; i++) {
//        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
//        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
//        pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
//#ifdef ADIABATIC
//      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
//         + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
//#endif
//      }
//    }
//  }
//#endif /* MHD */


// enroll new history variables, only once
// not implemented yet

  return;
}

//------------------------------------------------------------------------------
//  ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 //   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX




