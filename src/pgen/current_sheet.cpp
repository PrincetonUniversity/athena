//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file current_sheet.cpp
//  \brief Problem generator for setting up a current sheet.
//
// Can only be run in 2D.  Input parameters are:
//========================================================================================
// C/C++ headers
#include <iostream>   // endl
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

void InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void InnerX2Z0(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

static double ran2(long int *idum);
//static bool first_time=true;
//static long int iseed;
//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  int isheet= pin->GetOrAddInteger("problem","isheet",1);
  // enroll boundary value function pointers
  //if (isheet == 1) {
    EnrollUserBoundaryFunction(INNER_X1, InnerX1);
    EnrollUserBoundaryFunction(OUTER_X1, OuterX1);
    EnrollUserBoundaryFunction(INNER_X2, InnerX2);
    //EnrollUserBoundaryFunction(INNER_X2, InnerX2Z0);
    EnrollUserBoundaryFunction(OUTER_X2, OuterX2);
  //}

  // enroll new history variables
  //AllocateUserHistoryOutput(1);
  //EnrollUserHistoryOutput(0, hst_BxBy, "<-BxBy>");
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief set up a current sheet.
//========================================================================================


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  //GridS *pGrid=(pDomain->Grid);
  //Prim1DS W;
  //Cons1DS U1d;
  //Real x1,x2,x3;
  //Real Bx,uflow,beta;

/* setup uniform ambient medium with current sheet */

  Real uflow = pin->GetOrAddReal("problem","uflow",0.0);
  Real beta  = pin->GetOrAddReal("problem","beta",10.0);
  Real amp   = pin->GetOrAddReal("problem","amp",1e-3);
  Real alpha = pin->GetOrAddReal("problem","alpha",100.0);
  int isheet= pin->GetOrAddInteger("problem","isheet",1);
  Real iso_cs = peos->GetIsoSoundSpeed();

  AthenaArray<Real> phiflux,phiflux0; //flux function
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  phiflux.NewAthenaArray(nx2,nx1);
  phiflux0.NewAthenaArray(nx2,nx1);

  if (isheet == 0) { // striped from athena4
    Real d = 1.0;
    Real P = beta;
    Real E = 0.0; //need to be modified
    Real Vx = 0.0;
    Real Vy = 0.0;
    Real Vz = 0.0;
    Real Bx = 0.0;
    Real By = 1.0;
    Real Bz = 0.0;
    Real rval;

//   Ensure a different initial random seed for each meshblock.
    long int iseed = -1 - gid;

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real x2 = pcoord->x2f(j);
      Real x1 = pcoord->x1f(i);
      Vx = uflow*cos(PI*x2);
      phydro->u(IDN,k,j,i) = d;
      rval = amp*(ran2(&iseed) - 0.5);
      phydro->u(IM1,k,j,i) = d*(Vx+rval);
      rval = amp*(ran2(&iseed) - 0.5);
      phydro->u(IM2,k,j,i) = d*(Vy+rval);
      rval = amp*(ran2(&iseed) - 0.5);
      phydro->u(IM3,k,j,i) = d*(Vz+rval);
      if(NON_BAROTROPIC_EOS)
        phydro->u(IEN,k,j,i) = E;
      pfield->b.x1f(k,j,i) = Bx;
      pfield->b.x2f(k,j,i) = By;
      pfield->b.x3f(k,j,i) = Bz;
      if (i == ie) pfield->b.x1f(k,j,i+1) = Bx;
      if (j == je) pfield->b.x2f(k,j+1,i) = By;
      if (k == ke && ke > ks) pfield->b.x3f(k+1,j,i) = Bz;
      if (x1 > 0.5 && x1 < 1.5) {
        pfield->b.x2f(k,j,i) = -By;
        if (j == je) pfield->b.x2f(k,j+1,i) = -By;
      }
    }}}
  }

  // Huang & Bhattacharjee (2010)
  if(isheet == 1) {
    //calculate the flux function
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie+1; i++) {
      Real x2 = pcoord->x2f(j);
      Real x1 = pcoord->x1f(i);
      phiflux0(j,i) = cos(PI*x1)*sin(2.0*PI*fabs(x2))/2.0/PI;
      phiflux(j,i)  = tanh(alpha*x2)*cos(PI*x1)*sin(2.0*PI*x2)/2.0/PI;
    }}}
    // Initialize interface B
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real dx2f = pcoord->dx2f(j);
      Real dx1f = pcoord->dx1f(i);
      Real Bx = (phiflux(j,i)-phiflux(j+1,i))/dx2f;
      Real Bz = (phiflux(j,i+1)-phiflux(j,i))/dx1f;
      pfield->b.x1f(k,j,i) = Bx;
      pfield->b.x2f(k,j,i) = Bz;
      pfield->b.x3f(k,j,i) = 0.0;
      if (i == ie) pfield->b.x1f(k,j,i+1) = (phiflux(j,ie+1)-phiflux(j+1,ie+1))/dx2f;
      if (j == je) pfield->b.x2f(k,j+1,i) = (phiflux(je+1,i+1)-phiflux(je+1,i))/dx1f;
      if (k == ke && ke > ks) pfield->b.x3f(k+1,j,i) = 0.0;
    }}}

    // Ensure a different initial random seed for each meshblock.
    long int iseed = -1 - gid;

    // Initialize cell-centered density/pressure
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real dx2f = pcoord->dx2f(j);
      //Real dx1f = pcoord->dx1f(i);
      Real Bx = 0.5*(pfield->b.x1f(j,i)+pfield->b.x1f(j,i+1));
      Real Bx0 = 0.5*(phiflux0(j,i)-phiflux0(j+1,i)+phiflux0(j,i+1)-phiflux0(j+1,i+1))/dx2f;
      //Real Bz = 0.5*(pfield->b.x2f(j,i)+pfield->b.x2f(j+1,i));
      //Real Bz0 = 0.5*(phiflux0(j,i+1)-phiflux0(j,i) +
      //                phiflux0(j+1,i+1)-phiflux0(j+1,i))/dx1f;
      Real temp = 0.5*SQR(iso_cs);// 3.0; set cs=sqrt(6) in input file;
      Real pres0 = 2.0*temp+2.5*SQR(0.25*PI*(phiflux0(j,i)+phiflux0(j+1,i)+phiflux0(j,i+1)+phiflux0(j+1,i+1)));
      Real pres = pres0+0.5*(SQR(Bx0)-SQR(Bx));
      Real den  = pres/SQR(iso_cs);
      phydro->u(IDN,k,j,i) = den;
      Real rval = amp*(ran2(&iseed) - 0.5);
      phydro->u(IM1,k,j,i) = den*rval;
      rval = amp*(ran2(&iseed) - 0.5);
      phydro->u(IM2,k,j,i) = den*rval;
      phydro->u(IM3,k,j,i) = 0.0;
    }}}
  }

  // Oishi & MacLow (2015)
  if(isheet == 2) {
  }
}

//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================
//void MeshBlock::UserWorkInLoop(void)
//{
//  // Ensure a different initial random seed for each meshblock.
//  if(first_time){
//    iseed = -1 - gid;
//	first_time = false;
//  }
//
//
//  for (int k=ks; k<=ke; ++k) {
//  for (int j=js; j<=je; ++j) {
//  for (int i=is; i<=ie; ++i) {
//	phydro->u(IM1,k,j,i) +=
//    Real temp = (g-1.0)*prim(IEN,k,j,i)/prim(IDN,k,j,i);
//    cons(IEN,k,j,i) -= dt*prim(IDN,k,j,i)*(temp - 10.0)/tau/(g-1.0);
//      }
//    }
//  }
//
//  return;
//}


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


//--------------------------------------------------------------------------------------
//! \fn void InnerX1()
//  \brief Sets boundary condition on left X boundary (iib)
//  conducting for b-field and non-penetrable+free-slipping for hydro
void InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke){
  // set primitive variables in inlet ghost zones
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
    for(int i=1; i<=(NGHOST); ++i){
      a(IDN,k,j,is-i) = a(IDN,k,j,is+(i-1));
      a(IVX,k,j,is-i) = -a(IVX,k,j,is+(i-1));
      //a(IVY,k,j,is-i) = a(IVY,k,j,is);
      a(IVY,k,j,is-i) = a(IVY,k,j,is+(i-1));
      a(IVZ,k,j,is-i) = 0.0;
    }
  }}

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=(NGHOST); ++i){
        b.x1f(k,j,is) = 0.0;
        b.x1f(k,j,is-i) = -b.x1f(k,j,is+i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,is-i) = b.x2f(k,j,is+(i-1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
          b.x3f(k,j,is-i) = 0.0;
      }
    }}
  }
}

//--------------------------------------------------------------------------------------
//! \fn void OuterX1()
//  \brief Sets boundary condition on right X boundary (oib)
//  conducting for b-field and non-penetrable+free-slipping for hydro
void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke){
  // set primitive variables in inlet ghost zones
  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
    for(int i=1; i<=(NGHOST); ++i){
      a(IDN,k,j,ie+i) = a(IDN,k,j,ie-(i-1));
      a(IVX,k,j,ie+i) = -a(IVX,k,j,ie-(i-1));
      //a(IVY,k,j,ie+i) = a(IVY,k,j,ie);
      a(IVY,k,j,ie+i) = a(IVY,k,j,ie-(i-1));
      a(IVZ,k,j,ie+i) = 0.0;
    }
  }}

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=(NGHOST)-1; ++i){
        b.x1f(k,j,ie+1) = 0.0;
        b.x1f(k,j,ie+i+1) = -b.x1f(k,j,ie-(i-1));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,ie+i) = b.x2f(k,j,ie-(i-1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
          b.x3f(k,j,ie+i) = 0.0;
      }
    }}
  }

}

//--------------------------------------------------------------------------------------
//! \fn void InnerX2()
//  \brief Sets boundary condition on right Y boundary (ijb)
//  conducting for b-field and non-penetrable+free-slipping for hydro
void InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke){
  // set primitive variables in inlet ghost zones
  for(int k=ks; k<=ke; ++k){
  for(int j=1; j<=(NGHOST); ++j){
    for(int i=is; i<=ie; ++i){
      a(IDN,k,js-j,i) = a(IDN,k,js+(j-1),i);
      a(IVX,k,js-j,i) = a(IVX,k,js,i);
      a(IVY,k,js-j,i) = -a(IVY,k,js+(j-1),i);
      a(IVZ,k,js-j,i) = 0.0;
    }
  }}

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=1; j<=(NGHOST); ++j){
      for(int i=is; i<=ie; ++i){
        b.x1f(k,js-j,i) = b.x1f(k,js+(j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,js,i) = 0.0;
        b.x2f(k,js-j,i) = -b.x2f(k,js+j,i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
          b.x3f(k,js-j,i) = 0.0;
      }
    }}
  }

}

//--------------------------------------------------------------------------------------
//! \fn void OuterX2()
//  \brief Sets boundary condition on right Y boundary (ojb)
//  conducting for b-field and non-penetrable+free-slipping for hydro
void OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke){
  // set primitive variables in inlet ghost zones
  for(int k=ks; k<=ke; ++k){
  for(int j=1; j<=(NGHOST); ++j){
    for(int i=is; i<=ie; ++i){
      a(IDN,k,je+j,i) = a(IDN,k,je-(j-1),i);
      a(IVX,k,je+j,i) = a(IVX,k,je,i);
      a(IVY,k,je+j,i) = -a(IVY,k,je-(j-1),i);
      a(IVZ,k,je+j,i) = 0.0;
    }
  }}

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=1; j<=(NGHOST); ++j){
      for(int i=is; i<=ie; ++i){
        b.x1f(k,je+j,i) = b.x1f(k,je-(j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST)-1; ++j) {
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,je+1,i) = 0.0;
        b.x2f(k,je+j+1,i) = -b.x2f(k,je-(j-1),i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
          b.x3f(k,je+j,i) = 0.0;
      }
    }}
  }

}

//--------------------------------------------------------------------------------------
//! \fn void InnerX2Z0()
//  \brief Sets boundary condition on left Y boundary (ijb)
//  apply symmetry along z=0 plane: Bx=0(reflecting),dBz/dz=0(copy); dvx/dz=0(copy),dvz/dz(reflecting)
void InnerX2Z0(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                   Real time, Real dt, int is, int ie, int js, int je, int ks, int ke){
  // set primitive variables in inlet ghost zones
  for(int k=ks; k<=ke; ++k){
  for(int j=1; j<=(NGHOST); ++j){
    for(int i=is; i<=ie; ++i){
      a(IDN,k,js-j,i) = a(IDN,k,js+(j-1),i);
      a(IVX,k,js-j,i) = a(IVX,k,js+(j-1),i);
      a(IVY,k,js-j,i) = -a(IVY,k,js+(j-1),i);
      a(IVZ,k,js-j,i) = 0.0;
    }
  }}

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=1; j<=(NGHOST); ++j){
      for(int i=is; i<=ie; ++i){
        b.x1f(k,js-j,i) = -b.x1f(k,js+(j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
        b.x2f(k,js,i) = 2.0*b.x2f(k,js+1,i)-b.x2f(k,js+2,i);
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,js-j,i) = -b.x2f(k,js+j,i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
          b.x3f(k,js-j,i) = 0.0;
      }
    }}
  }

}
