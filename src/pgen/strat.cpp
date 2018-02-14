//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
/*============================================================================*/
/*! \file strat.c
 *  \brief Problem generator for stratified 3D shearing sheet.
 *
 * PURPOSE:  Problem generator for stratified 3D shearing sheet.  Based on the
 *   initial conditions described in "Three-dimensional Magnetohydrodynamic
 *   Simulations of Vertically Stratified Accretion Disks" by Stone, Hawley,
 *   Gammie & Balbus.
 *
 * Several different field configurations and perturbations are possible:
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *  ifield = 4 - uniform By, but only for |z|<2
 *  ifield = 5 - By with constant \beta versus z
 *  ifield = 6 - zero field everywhere
 *
 * - ipert = 1 - random perturbations to P and V [default, used by HGB]
 *
 * Code must be configured using -sh
 *
 * REFERENCE: Stone, J., Hawley, J. & Balbus, S. A., ApJ 463, 656-673 (1996)
 *            Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992)	      */
/*============================================================================*/

// C/C++ headers
#include <iostream>
#include <cfloat>
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


static double ran2(long int *idum);
void VertGrav(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons);
void StratOutflowInnerX3(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &a, 
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke);
void StratOutflowOuterX3(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &a,
                  FaceField &b, Real time, Real dt, 
                  int is, int ie, int js, int je, int ks, int ke);
static Real hst_BxBy(MeshBlock *pmb, int iout);
static Real hst_dVxVy(MeshBlock *pmb, int iout);

static Real Omega_0,qshear;

/* Apply a density floor - useful for large |z| regions */
static Real dfloor,pfloor;


//====================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, hst_BxBy, "-BxBy");
  EnrollUserHistoryOutput(1, hst_dVxVy, "dVxVy");
// Read problem parameters
  Omega_0 = pin->GetOrAddReal("problem","Omega0",1.0e-3);
  qshear  = pin->GetOrAddReal("problem","qshear",1.5);

// Enroll user-defined physical source terms
//   vertical external gravitational potential
  EnrollUserExplicitSourceFunction(VertGrav);

// enroll user-defined boundary conditions
  if(mesh_bcs[INNER_X3] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(INNER_X3, StratOutflowInnerX3);
  }
  if(mesh_bcs[OUTER_X3] == GetBoundaryFlag("user")) {
      EnrollUserBoundaryFunction(OUTER_X3, StratOutflowOuterX3);
  }

  return;
}



//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief stratified disk problem generator for 3D problems.
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  int ifield, ipert;
  Real beta, amp, pres;
  Real iso_cs=1.0;
  Real B0 = 0.0;

  Real SumRvx=0.0, SumRvy=0.0, SumRvz=0.0;
  Real x1,x2,x3,xmin,xmax;
  Real x1f,x2f,x3f;
  Real rd,rp,rvx,rvy,rvz,rbx,rby,rbz;
  Real rval;

  // initialize density
  const Real den=1.0;

  // Initialize boxsize
  Real Lx = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  Real Ly = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  Real Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  // initialize wavenumbers
  int nwx = pin->GetOrAddInteger("problem","nwx",1);
  int nwy = pin->GetOrAddInteger("problem","nwy",1);
  int nwz = pin->GetOrAddInteger("problem","nwz",1);
  Real kx = (2.0*PI/Lx)*((double)nwx);// nxw=-ve for leading wave
  Real ky = (2.0*PI/Ly)*((double)nwy);
  Real kz = (2.0*PI/Lz)*((double)nwz);

  // Ensure a different initial random seed for each meshblock.
  long int iseed = -1 - gid;

  // adiabatic gamma
  Real gam = peos->GetGamma();

  if (pmy_mesh->mesh_size.nx3 == 1){
    std::cout << "[strat.cpp]: Strat only works on a 3D grid" 
      << std::endl;
  }

  // Read problem parameters for initial conditions
  amp = pin->GetReal("problem","amp");
  ipert = pin->GetOrAddInteger("problem","ipert", 1);
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pfloor=pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));

  if (MAGNETIC_FIELDS_ENABLED) {
    ifield = pin->GetOrAddInteger("problem","ifield", 1);
    beta = pin->GetReal("problem","beta");
  }
  // Compute pressure based on the EOS.
  if (NON_BAROTROPIC_EOS) {
    pres  = pin->GetOrAddReal("problem","pres",1.0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    pres = den*SQR(iso_cs);
  }

  //Compute field strength based on beta.
  if (MAGNETIC_FIELDS_ENABLED) {
    B0 = sqrt((double)(2.0*pres/beta));
    std::cout << "B0=" << B0 << std::endl;
  }

  // With viscosity and/or resistivity, read eta_Ohm and nu_V
  // to be filled in
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
        // ipert = 1 - random perturbations to P/d and V 
        // [default, used by HGB]
        if (ipert == 1) {
          rval = amp*(ran2(&iseed) - 0.5);
          rd = den*exp(-x3*x3)*(1.0+2.0*rval);
          if (rd < dfloor) rd = dfloor;
          if (NON_BAROTROPIC_EOS) {
            rp = pres/den*rd;
            if (rp < pfloor) rp = pfloor;
          }
          rval = amp*(ran2(&iseed) - 0.5);
          rvx = (0.4/sqrt(3.0)) *rval*1e-3;
          SumRvx += rvx;

          rval = amp*(ran2(&iseed) - 0.5);
          rvy = (0.4/sqrt(3.0)) *rval*1e-3;
          SumRvy += rvz;

          rval = amp*(ran2(&iseed) - 0.5);
          rvz = 0.4*rval*sqrt(pres/den);
          rvz = (0.4/sqrt(3.0)) *rval*1e-3;
          SumRvz += rvz;
        // no perturbations
        } else {
          rd = den*exp(-x3*x3);
          rvx = 0;
          rvy = 0;
          rvz = 0;
        }

        // Initialize d, M, and P.
        // for_the_future: if FARGO do not initialize the bg shear
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) = rd*rvy;
        phydro->u(IM2,k,j,i) -= rd*(qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = rd*rvz;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = rp/(gam-1.0)
            + 0.5*(SQR(phydro->u(IM1,k,j,i))
                + SQR(phydro->u(IM2,k,j,i))
                + SQR(phydro->u(IM3,k,j,i)))/rd;
        } // Hydro

        // Initialize magnetic field.  For 3D shearing box B1=Bx, B2=By, B3=Bz
        //  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
        //  ifield = 2 - uniform Bz
        //  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
        //  ifield = 4 - uniform By, but only for |z|<2
        //  ifield = 5 - By with constant \beta versus z
        //  ifield = 6 - zero field everywhere
        if (MAGNETIC_FIELDS_ENABLED) {
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
            pfield->b.x3f(k,j,i) = B0;
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
            if (k==ke) pfield->b.x3f(ke+1,j,i) = B0;
          }
          if (ifield == 3) {
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = B0*(cos((double)kx*x1));
            pfield->b.x3f(k,j,i) = B0*(sin((double)kx*x1));
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = B0*(cos((double)kx*x1));
            if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*(sin((double)kx*x1));
          }
          if (ifield == 4 && fabs(x3) < 2.0) {
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = B0;
            pfield->b.x3f(k,j,i) = 0.0;
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = B0;
            if (k==ke) pfield->b.x3f(ke+1,j,i) = 0.0;
          }
          if (ifield == 5) {
            /* net toroidal field with constant \beta with height */
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = sqrt(den*exp(-x3*x3)*SQR(Omega_0)/beta);
            pfield->b.x3f(k,j,i) = 0.0;
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = sqrt(den*exp(-x3*x3)*
                                                    SQR(Omega_0)/beta);
            if (k==ke) pfield->b.x3f(ke+1,j,i) = 0.0;
          }
          if (ifield == 6) {
            /* zero field everywhere */
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = 0.0;
            pfield->b.x3f(k,j,i) = 0.0;
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
            if (k==ke) pfield->b.x3f(ke+1,j,i) = 0.0;
          }
        } // MHD
      }
    }
  }

  /* For random perturbations as in HGB, ensure net momentum is zero by
   * subtracting off mean of perturbations */

  if (ipert == 1) {
    int cell_num = block_size.nx1*block_size.nx2*block_size.nx3;
    SumRvx /= cell_num;
    SumRvy /= cell_num;
    SumRvz /= cell_num;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IM1,k,j,i) -= rd*SumRvx;
          phydro->u(IM2,k,j,i) -= rd*SumRvy;
          phydro->u(IM3,k,j,i) -= rd*SumRvz;
        }
      }
    }
  }

  return;
}


void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  return;
}

void MeshBlock::UserWorkInLoop(void)
{
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real& u_d  = phydro->u(IDN,k,j,i);
        u_d = (u_d > dfloor) ?  u_d : dfloor;
        if (NON_BAROTROPIC_EOS) {
          Real gam = peos->GetGamma();
          Real& w_p  = phydro->w(IPR,k,j,i);
          Real& u_e  = phydro->u(IEN,k,j,i);
          Real& u_m1 = phydro->u(IM1,k,j,i);
          Real& u_m2 = phydro->u(IM2,k,j,i);
          Real& u_m3 = phydro->u(IM3,k,j,i);
          w_p = (w_p > pfloor) ?  w_p : pfloor;
          Real di = 1.0/u_d;
          Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
          u_e = w_p/(gam-1.0)+ke;
      }
    }
  }}
  return;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  return;
}

/*------------------------------------------------------------------------------
 */

#define IMR1 2147483563
#define IMR2 2147483399
#define AM (1.0/IMR1)
#define IMM1 (IMR1-1)
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
      if (*idum < 0) *idum += IMR1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IMR1 without */
  if (*idum < 0) *idum += IMR1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IMR2 likewise*/
  if (idum2 < 0) idum2 += IMR2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IMR1
#undef IMR2
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

void VertGrav(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons)
{
  Real fsmooth, xi, sign;
  Real Lz = pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min;
  Real z0 = Lz/2.0;
  Real lambda = 0.1 / z0;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real x3 = pmb->pcoord->x3v(k);
        //smoothing function
        if (x3 >= 0) {
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        xi = z0/x3;
        fsmooth = SQR( sqrt( SQR(xi+sign) + SQR(xi*lambda) ) + xi*sign );
        //multiply gravitational potential by smoothing function
        cons(IM3,k,j,i) -= dt*den*SQR(Omega_0)*x3*fsmooth;
        if (NON_BAROTROPIC_EOS) {
          cons(IEN,k,j,i) -= dt*den*SQR(Omega_0)*prim(IVZ,k,j,i)*x3*fsmooth;
        }
      }
    }
  }
  return;
}

 //Here is the lower z outflow boundary.
 //  The basic idea is that the pressure and density
 //  are exponentially extrapolated in the ghost zones
 //  assuming a constant temperature there (i.e., an
 //  isothermal atmosphere). The z velocity (NOT the
 //  momentum) are set to zero in the ghost zones in the
 //  case of the last lower physical zone having an inward
 //  flow.  All other variables are extrapolated into the
 //  ghost zones with zero slope.

void StratOutflowInnerX3(MeshBlock *pmb, Coordinates *pco,
    AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js,
    int je, int ks, int ke)
{

  /* Copy field components from last physical zone
   * zero slope boundary for B field*/
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=NGHOST; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          b.x1f(ks-k,j,i) = b.x1f(ks,j,i);
        }
      }
    }
    for (int k=1; k<=NGHOST; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          b.x2f(ks-k,j,i) = b.x2f(ks,j,i);
        }
      }
    }
    for (int k=1; k<=NGHOST; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          b.x3f(ks-k,j,i) = b.x3f(ks,j,i);
        }
      }
    }
  } // MHD

  for (int k=1; k<=NGHOST; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x3 = pco->x3v(ks-k);
        Real x3b = pco->x3v(ks);
        Real den = prim(IDN,ks,j,i);
        /* First calculate the effective gas temperature (Tks=cs^2)
         * in the last physical zone. If isothermal, use H=1 */
        Real Tks = 0.5*SQR(Omega_0);
        if (NON_BAROTROPIC_EOS) {
          Real pressks = prim(IPR,ks,j,i);
          pressks = std::max(pressks,pfloor);
          Tks = pressks/den;
        }
        /* Now extrapolate the density to balance gravity
         * assuming a constant temperature in the ghost zones */
        prim(IDN,ks-k,j,i) = den*exp(-(SQR(x3)-SQR(x3b))/
                                (2.0*Tks/SQR(Omega_0)));
        /* Copy the velocities, but not the momenta --- 
         * important because of the density extrapolation above */
        prim(IVX,ks-k,j,i) = prim(IVX,ks,j,i);
        prim(IVY,ks-k,j,i) = prim(IVY,ks,j,i);
        /* If there's inflow into the grid, set the normal velocity to zero */
        if (prim(IVZ,ks,j,i) >= 0.0) {
          prim(IVZ,ks-k,j,i) = 0.0;
        } else {
          prim(IVZ,ks-k,j,i) = prim(IVZ,ks,j,i);
        }
        if (NON_BAROTROPIC_EOS)
          prim(IPR,ks-k,j,i) = prim(IDN,ks-k,j,i)*Tks;
      }
    }
  }

  return;

}

 //Here is the upper z outflow boundary.
 // The basic idea is that the pressure and density
 // are exponentially extrapolated in the ghost zones
 // assuming a constant temperature there (i.e., an
 // isothermal atmosphere). The z velocity (NOT the
 // momentum) are set to zero in the ghost zones in the
 // case of the last upper physical zone having an inward
 // flow.  All other variables are extrapolated into the
 // ghost zones with zero slope.
void StratOutflowOuterX3(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke)
{
/* Copy field components from last physical zone */
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=NGHOST; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          b.x1f(ke+k,j,i) = b.x1f(ke,j,i);
        }
      }
    }
    for (int k=1; k<=NGHOST; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          b.x2f(ke+k,j,i) = b.x2f(ke,j,i);
        }
      }
    }
    for (int k=1; k<=NGHOST; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          b.x3f(ke+1+k,j,i) = b.x3f(ke+1,j,i);
        }
      }
    }
  } // MHD

  for (int k=1; k<=NGHOST; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x3 = pco->x3v(ke+k);
        Real x3b = pco->x3v(ke);
        Real den = prim(IDN,ke,j,i);
        /* First calculate the effective gas temperature (Tks=cs^2)
         * in the last physical zone. If isothermal, use H=1 */
        Real Tke = 0.5*SQR(Omega_0);
        if (NON_BAROTROPIC_EOS) {
          Real presske = prim(IPR,ke,j,i);
          presske = std::max(presske,pfloor);
          Real Tke = presske/den;
        }
        /* Now extrapolate the density to balance gravity 
         * assuming a constant temperature in the ghost zones */
        prim(IDN,ke+k,j,i) = den*exp(-(SQR(x3)-SQR(x3b))/
                               (2.0*Tke/SQR(Omega_0)));
        /* Copy the velocities, but not the momenta --- 
         * important because of the density extrapolation above */
        prim(IVX,ke+k,j,i) = prim(IVX,ke,j,i);
        prim(IVY,ke+k,j,i) = prim(IVY,ke,j,i);
        /* If there's inflow into the grid, set the normal velocity to zero */
        if (prim(IVZ,ke,j,i) <= 0.0) {
          prim(IVZ,ke+k,j,i) = 0.0;
        } else {
          prim(IVZ,ke+k,j,i) = prim(IVZ,ke,j,i);
        }
        if (NON_BAROTROPIC_EOS)
          prim(IPR,ke+k,j,i) = prim(IDN,ke+k,j,i)*Tke;
      }
    }
  }
  return;
}

static Real hst_BxBy(MeshBlock *pmb, int iout)
{
  Real bxby=0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  volume.NewAthenaArray(ncells1);

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,volume);
      for(int i=is; i<=ie; i++) {
        bxby-=volume(i)*b(IB1,k,j,i)*b(IB2,k,j,i);
      }
    }
  }
  volume.DeleteAthenaArray();

  return bxby;
}

static Real hst_dVxVy(MeshBlock *pmb, int iout)
{
  Real dvxvy=0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vshear=0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  volume.NewAthenaArray(ncells1);

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,volume);
      for(int i=is; i<=ie; i++) {
        vshear = qshear*Omega_0*pmb->pcoord->x1v(i);
        dvxvy+=volume(i)*w(IDN,k,j,i)*w(IVX,k,j,i)*(w(IVY,k,j,i)+vshear);
      }
    }
  }

  volume.DeleteAthenaArray();
  return dvxvy;
}
