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
//! \file hb3.c
//  \brief Problem generator for 2D MRI simulations using the shearing sheet
//   based on "A powerful local shear instability in weakly magnetized disks.
//
// * PURPOSE: Problem generator for 2D MRI simulations using the shearing sheet
// *   based on "A powerful local shear instability in weakly magnetized disks.
// *   III - Long-term evolution in a shearing sheet" by Hawley & Balbus.  This
// *   is the third of the HB papers on the MRI, thus hb3.
// *
// * Several different perturbations and field configurations are possible:
// * - ipert = 1 - isentropic perturbations to P & d [default]
// * - ipert = 2 - uniform Vx=amp, sinusoidal density
// *
// * - ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
// * - ifield = 2 - uniform Bz
// *
// * PRIVATE FUNCTION PROTOTYPES:
// * - ran2() - random number generator from NR
// *
// * REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).*/
//
//======================================================================================
//
//[JMSHI
//
// C++ headers
#include <iostream>   // cout, endl
#include <stdlib.h>   // exit
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
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
#if !SHEARING_BOX
#error "This problem generator requires shearing box"
#endif

static Real amp, nwx, nwy; // amplitude, Wavenumbers
static int ShBoxCoord, ipert,ifield; // initial pattern
static Real beta,B0,pres;
static Real gm1,iso_cs;
static Real x1size,x2size,x3size;
static Real Omega_0,qshear;
static int nx1,nx2,nvar;
static AthenaArray<Real> ibval,obval; // ghost cells array
static int first_time=1;
AthenaArray<Real> volume; // 1D array of volumes


static double ran2(long int *idum);
static Real hst_BxBy(MeshBlock *pmb, int iout);
void ShearInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ShearOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // initialize global variables
  amp    = pin->GetReal("problem","amp");
  beta   = pin->GetReal("problem","beta");
  nwx = pin->GetOrAddInteger("problem","nwx",1);
  nwy = pin->GetOrAddInteger("problem","nwy",1);
  ShBoxCoord = pin->GetOrAddInteger("problem","shboxcoord",2);
  ipert  = pin->GetOrAddInteger("problem","ipert",1);
  ifield = pin->GetOrAddInteger("problem","ifield",1);
  Omega_0= pin->GetOrAddReal("problem","Omega0",0.001);
  qshear = pin->GetOrAddReal("problem","qshear",1.5);

  // enroll boundary value function pointers
  //EnrollUserBoundaryFunction(INNER_X1, ShearInnerX1);
  //EnrollUserBoundaryFunction(OUTER_X1, ShearOuterX1);

  // enroll new history variables
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, hst_BxBy, "<-BxBy>");
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  if (pmy_mesh->mesh_size.nx2 == 1 || pmy_mesh->mesh_size.nx3 > 1) {
      std::cout << "[hb3.cpp]: only works on 2D grid" << std::endl;
      exit(0);
  }


  if (ShBoxCoord != 2) {
      std::cout << "[hb3.cpp]: only works for x-z plane with ShBoxCoord = 2" << std::endl;
      exit(0);
  }
  // allocate 1D array for cell volume used in usr def history
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  volume.NewAthenaArray(ncells1);

  // Initialize boundary value arrays
  if (first_time) {
    nx1 = (ie-is)+1 + 2*(NGHOST);
    nx2 = (je-js)+1 + 2*(NGHOST);
    nvar = (NHYDRO+NFIELD);  // for now IDN, IVX, IVY, IVZ, NHYDRO,NHYDRO+1,+2
    ibval.NewAthenaArray(nvar,nx2,(NGHOST));
    obval.NewAthenaArray(nvar,nx2,(NGHOST));

    first_time = 0;
  }
  Real d0 = 1.0;
  Real p0 = 1e-5;

  if (NON_BAROTROPIC_EOS) {
    gm1 = (peos->GetGamma() - 1.0);
    iso_cs = sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
  }

  B0 = sqrt((double)(2.0*p0/beta));
  std::cout << "iso_cs = " << iso_cs << std::endl;
  std::cout << "gamma  = " << peos->GetGamma() << std::endl;
  std::cout << "d0     = " << d0     << std::endl;
  std::cout << "p0     = " << p0     << std::endl;
  std::cout << "B0     = " << B0     << std::endl;
  std::cout << "ipert  = " << ipert  << std::endl;
  std::cout << "ifield = " << ifield << std::endl;
  std::cout << "beta   = " << beta   << std::endl;


  x1size = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  x2size = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  x3size = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;
  std::cout << "[hb3.cpp]: [Lx,Lz,Ly] = [" <<x1size <<","<<x2size<<","<<x3size<<"]"<<std::endl;

  Real kx = (2.0*PI/x1size)*((double)nwx);
  Real kz = (2.0*PI/x2size)*((double)nwy);

  Real x1,x2,x3,rd,rp,rval, rvx, rvy, rvz;
  long int iseed = -1-gid; /* Initialize on the first call to ran2 */
// Initialize perturbations
// *  ipert = 1 - isentropic perturbations to P & d [default]
// *  ipert = 2 - uniform Vx=amp, sinusoidal density
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      x1 = pcoord->x1v(i);
      x2 = pcoord->x2v(j);
      rd = d0;
      rp = p0;
      rvx = 0.0;
      rvy = 0.0;
      rvz = 0.0;
      if (ipert == 1) {
        rval = 1.0 + amp*(ran2(&iseed) - 0.5);
        if (NON_BAROTROPIC_EOS) {
          rp = rval*p0;
          rd = d0;
        } else {
          rd = rval*d0;
        }
        rvx = 0.0;
      } else if (ipert == 2) {
        rp = p0;
        rd = d0*(1.0+0.1*sin((double)kx*x1));
        if (NON_BAROTROPIC_EOS) {
          rvx = amp*sqrt((gm1+1.0)*p0/d0);
        } else {
          rvx = amp*sqrt(p0/d0);
        }
      } else {
          std::cout << "[hb3.cpp] ipert = " <<ipert <<" is unrecognized " <<std::endl;
          exit(0);
      }
      phydro->u(IDN,ks,j,i) = rd;
      phydro->u(IM1,ks,j,i) = rd*rvx;
      phydro->u(IM2,ks,j,i) = rd*rvy;
      phydro->u(IM3,ks,j,i) = rd*rvz;
      phydro->u(IM3,ks,j,i) -= rd*qshear*Omega_0*x1;
      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,ks,j,i) = rp/gm1 +
             0.5*(SQR(phydro->u(IM1,ks,j,i)) +
                  SQR(phydro->u(IM2,ks,j,i)) +
                  SQR(phydro->u(IM3,ks,j,i)))/rd;
      }

// Initialize magnetic field.  For 2D shearing box
// B1=Bx, B2=Bz, B3=By
// ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
// ifield = 2 - uniform Bz
// ifield = 3 - sinusiodal modes (Nordita workshop test)
      if (MAGNETIC_FIELDS_ENABLED) {
        if (ifield == 1) {
          pfield->b.x1f(ks,j,i) = 0.0;
          pfield->b.x2f(ks,j,i) = B0*(sin((double)kx*x1));
          pfield->b.x3f(ks,j,i) = 0.0;
          if (i==ie) pfield->b.x1f(ks,j,ie+1) = 0.0;
          if (j==je) pfield->b.x2f(ks,je+1,i) = B0*(sin((double)kx*x1));
        } else if (ifield == 2) {
            pfield->b.x1f(ks,j,i) = 0.0;
            pfield->b.x2f(ks,j,i) = B0;
            pfield->b.x3f(ks,j,i) = 0.0;
            if (i==ie) pfield->b.x1f(ks,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(ks,je+1,i) = B0;
        } else {
            std::cout << "[hb3.cpp] ifield = " <<ifield <<" is unrecognized " <<std::endl;
            exit(0);
        }
          if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,ks,j,i) += 0.5*(
                  SQR(0.5*(pfield->b.x1f(ks,j,i) + pfield->b.x1f(ks,j,i+1))) +
                  SQR(0.5*(pfield->b.x2f(ks,j,i) + pfield->b.x2f(ks,j+1,i))) +
                  SQR(pfield->b.x3f(ks,j,i)));
        }
      }

      }
    }

  return;
}


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================
void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void ShearInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for ssheet problem

void ShearInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  Real qomL = qshear*Omega_0*x1size;

  //// Initialize boundary value arrays
  //int nx1 = (ie-is)+1 + 2*(NGHOST);
  //int nx2 = (je-js)+1 + 2*(NGHOST);
  //int nvar = (NHYDRO+NFIELD);  // for now IDN, IVX, IVY and IVZ, and IEN if non_barotropic
  //AthenaArray<Real> bval;
  //bval.NewAthenaArray(nvar,nx2,(NGHOST));
  //int nyzone = (je-js)+1;

  // set bval variables in inlet ghost zones
  for(int j=0; j<nx2; ++j) {
    for(int i=1; i<=(NGHOST); ++i) {
      ibval(IDN,j,i-1) = a(IDN,ks,j,ie-(NGHOST)+i);
      ibval(IVX,j,i-1) = a(IVX,ks,j,ie-(NGHOST)+i);
      ibval(IVY,j,i-1) = a(IVY,ks,j,ie-(NGHOST)+i);
      ibval(IVZ,j,i-1) = a(IVZ,ks,j,ie-(NGHOST)+i);
      if (NON_BAROTROPIC_EOS) {
        ibval(IEN,j,i-1) = a(IEN,ks,j,ie-(NGHOST)+i);
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        ibval(NHYDRO,j,i-1) = b.x1f(ks,j,ie-(NGHOST)+i);
        ibval(NHYDRO+1,j,i-1) = b.x2f(ks,j,ie-(NGHOST)+i);
        ibval(NHYDRO+2,j,i-1) = b.x3f(ks,j,ie-(NGHOST)+i);
      }
    }
  }

  for(int j=0; j<nx2; ++j) {
    for(int i=1; i<=(NGHOST); ++i) {
      int ib = (NGHOST) - i;
      a(IDN,ks,j,is-i) = ibval(IDN,j,ib);
      a(IVX,ks,j,is-i) = ibval(IVX,j,ib);
      a(IVY,ks,j,is-i) = ibval(IVY,j,ib);
      a(IVZ,ks,j,is-i) = ibval(IVZ,j,ib)+qshear*Omega_0*x1size;
      if (NON_BAROTROPIC_EOS) {
        a(IEN,ks,j,is-i) = ibval(IEN,j,ib);
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        b.x1f(ks,j,is-i) = ibval(NHYDRO,j,ib);
        b.x2f(ks,j,is-i) = ibval(NHYDRO+1,j,ib);
        b.x3f(ks,j,is-i) = ibval(NHYDRO+2,j,ib);
      }
    }}

}

//--------------------------------------------------------------------------------------
//! \fn void ShearOuterX1()
//  \brief Sets boundary condition on right X boundary (oib) for ssheet problem

void ShearOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

//  // Initialize boundary value arrays
//  int nx1 = (ie-is)+1 + 2*(NGHOST);
//  int nx2 = (je-js)+1 + 2*(NGHOST);
//  int nvar = (NHYDRO);  // for now IDN, IVX, IVY and IVZ
//  AthenaArray<Real> bval;
//  bval.NewAthenaArray(nvar,nx2,(NGHOST));
//  AthenaArray<Real> dbval;
//  dbval.NewAthenaArray(nvar,nx2,(NGHOST));
//  int nyzone = (je-js)+1;

  // set primitive variables in inlet ghost zones
  for(int j=0; j<nx2; ++j) {
    for(int i=1; i<=(NGHOST); ++i) {
      obval(IDN,j,i-1) = a(IDN,ks,j,is+i-1);
      obval(IVX,j,i-1) = a(IVX,ks,j,is+i-1);
      obval(IVY,j,i-1) = a(IVY,ks,j,is+i-1);
      obval(IVZ,j,i-1) = a(IVZ,ks,j,is+i-1);
      if (NON_BAROTROPIC_EOS) {
        obval(IEN,j,i-1) = a(IEN,ks,j,is+i-1);
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        obval(NHYDRO,j,i-1) = b.x1f(ks,j,is+i-1);
        obval(NHYDRO+1,j,i-1) = b.x2f(ks,j,is+i-1);
        obval(NHYDRO+2,j,i-1) = b.x3f(ks,j,is+i-1);
      }
    }
  }


  for(int j=0; j<=nx2; ++j) {
    for(int i=1; i<=(NGHOST); ++i) {
      int ib =  i - 1;
      a(IDN,ks,j,ie+i) = obval(IDN,j,ib);
      a(IVX,ks,j,ie+i) = obval(IVX,j,ib);
      a(IVY,ks,j,ie+i) = obval(IVY,j,ib);
      a(IVZ,ks,j,ie+i) = obval(IVZ,j,ib)-qshear*Omega_0*x1size;
      if (NON_BAROTROPIC_EOS) {
        a(IEN,ks,j,ie+i) = obval(IEN,j,ib);
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        b.x1f(ks,j,ie+i) = obval(NHYDRO,j,ib);
        b.x2f(ks,j,ie+i) = obval(NHYDRO+1,j,ib);
        b.x3f(ks,j,ie+i) = obval(NHYDRO+2,j,ib);
      }
    }}

}

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

/*! \fn double ran2(long int *idum)
 *  \brief Extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
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


static Real hst_BxBy(MeshBlock *pmb, int iout)
{
  Real bxby=0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,volume);
      for(int i=is; i<=ie; i++) {
        bxby-=volume(i)*b(IB1,k,j,i)*b(IB3,k,j,i);
      }
    }
  }

  return bxby;
}


//JMSHI]
