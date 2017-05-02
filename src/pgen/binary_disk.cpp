//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file binary_disk.cpp
//  \brief Initializes accretion disk around binary in cartesian cooordinates(for now)

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min
#include <cstdlib>    // srand
#include <cfloat>     // FLT_MIN

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

void solve_u(const Real t, const Real dto, const Real e, Real *ua);
void compute_f(Real ua, const Real e, Real *f_ptr);
void compute_star_loc(const Real dto, const Real t, const Real e,
					  const Real mf, Real *f_ptr, Real *u_ptr, Real *r_ptr);

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
static Real PoverR(const Real rad, const Real phi, const Real z);
static void VelProfileCyl(const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3);
static Real RampFunc(const Real rad, const Real phi, const Real z,
					 const Real v1, const Real v2, const Real v3);

//user defined src term
//void Cooling(MeshBlock *pmb, const Real time, const Real dt,
//             const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
//             AthenaArray<Real> &cons);
void Binary(MeshBlock *pmb, const Real time, const Real dt,
			const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
			AthenaArray<Real> &cons);

// problem parameters which are useful to make global to this file
static Real semia,ecc,qrat,mu,incli,argp;
static Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
static Real dfloor,pfloor;
static Real rsoft,rsink,rin,rout,rbuf1,rbuf2;
static Real tsink; // mass removing time scale
static Real tbin;  // binary period 2pi/(GM/a^3)^1/2
static Real beta_th,alpha;
// parameters for computing binary orbit
static int itermax=100;
static Real dueps=1e-6;
//instantaneous binary separation, ecc and true anormaly in its orbit plane
static Real rsep,uanorm,fanorm;
//real-time binary postion in cartesian
static Real x1s,x2s,x3s,x1p,x2p,x3p;

// debug for binary orbit
static FILE * pf;
static Real tstart;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Get parameters for binary potential
  // gm0:      total mass of binary
  // qrat:     mass ratio secondary/primary (<=1 always)
  // mu:       mass fraction of the primary defined as 1/(1+qrat)
  // semia:    semi-major axis(also the initial separation)
  // incli:    inclination(angle between disk and binary orbit planes)
  // ecc:      eccentrcity
  // phi0_2:   initial phase angle of the secondary (0 as default)
  // phi0_1:   initial phase angle of the primary(PI as default)
  // argp:     argument of periapse, angle between the axis of periapse and line of nodes
  // lon:      angle between line of nodes and reference X-axis.Set to be zero always.
  gm0  = pin->GetOrAddReal("problem","GMb",1.0);
  qrat = pin->GetOrAddReal("problem","ratio",1.0);
  mu = 1.0/(1.0+qrat);
  ecc  = pin->GetOrAddReal("problem","ecc",0.0);
  incli  = pin->GetOrAddReal("problem","inclination",0.0);
  argp   = pin->GetOrAddReal("problem","periapse",0.0);
  semia  = pin->GetOrAddReal("problem","semi_major",1.0);
  tbin   = 2.0*PI/sqrt(gm0/semia/SQR(semia));
  //phi0_2 = pin->GetOrAddReal("problem","phi0",0.0);
  //phi0_1 = phi0_2+PI;
  rsoft = pin->GetOrAddReal("problem","rsoft",0.05);
  rsink = pin->GetOrAddReal("problem","rsink",0.05);
  rsoft *= semia;
  rsink *= semia;
  beta_th = pin->GetOrAddReal("problem","beta_th",0.0);

  // need to initialize binary position or from restart point
  Real dto = 0.01;
  Real torb = fmod(time,tbin);
  rsep= semia; uanorm=0.0; fanorm=0.0; //separation, ecc and true anormaly
  compute_star_loc(dto,torb,ecc,mu,&fanorm,&uanorm,&rsep);
  //std::cout << "[init_cond]: dt,t,e,mu,f,u,r = "<< dto << " " << torb << " " << ecc <<" " << mu << " " << fanorm << " " << uanorm << " " << rsep << " " << std::endl;
  Real rcos = rsep*cos(fanorm+argp);
  Real rsin = rsep*sin(fanorm+argp);
  Real rsincos = rsin*cos(incli);
  Real rsinsin = rsin*sin(incli);
  x1p = - (1.0-mu)*rcos;
  x2p = - (1.0-mu)*rsincos;
  x3p = (1.0-mu)*rsinsin;
  x1s = mu*rcos;
  x2s = mu*rsincos;
  x3s =-mu*rsinsin;

  // Get parameters for circumbinary disk
  r0 = semia; // set length unit
  rin = pin->GetOrAddReal("problem","rin",0.0);
  rout = pin->GetOrAddReal("problem","rout",5.0);
  tsink = pin->GetOrAddReal("problem","tsink",100.0);
  rbuf1 = pin->GetOrAddReal("problem","rbuf1",8.0);
  rbuf2 = pin->GetOrAddReal("problem","rbuf2",10.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters of initial pressure and cooling parameters
  if(NON_BAROTROPIC_EOS){
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  }else{
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pfloor=pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
  int inu = pin->GetOrAddReal("problem","inu",1);
  if (inu != 0) alpha = pin->GetOrAddReal("problem","nuiso",0.1);
  else          alpha = pin->GetOrAddReal("problem","nuiso",2.5e-4)/p0_over_r0;

  // Enroll user-defined physical source terms
  if (qrat != 0.0) EnrollUserExplicitSourceFunction(Binary);

  // debug binary orbit
  pf = fopen ("binary_orbit.tab","w");
  tstart = torb; //first dump
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real rad, phi, z;
  Real v1, v2, v3;

  //  Initialize density and momenta
  for(int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
	  // 1) pure power law with truncation radii rin rout
      //if (rad >= rin && rad <= rout)
      //if (rad >= rin)
      //  phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
      //else
      //  phydro->u(IDN,k,j,i) = dfloor;
	  // 2) power law profile with exponential taper interior to rin
	  phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z)*exp(-pow((rad/rin),-10.0));
      VelProfileCyl(rad,phi,z,v1,v2,v3);
      phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
      phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
      phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
      if (NON_BAROTROPIC_EOS){
        Real p_over_r = PoverR(rad,phi,z);
        phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                   + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k)
{
  rad=sqrt(SQR(pco->x1v(i))+SQR(pco->x2v(j)));
  phi=atan2(pco->x2v(j),pco->x1v(i));
  z=pco->x3v(k);
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

static Real DenProfileCyl(const Real rad, const Real phi, const Real z)
{
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*pow(rad/r0,dslope);
  Real dentem = denmid*exp(gm0/p_over_r*(1./sqrt(SQR(rad)+SQR(rsoft)+SQR(z))-1./sqrt(SQR(rad)+SQR(rsoft))));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates

static Real PoverR(const Real rad, const Real phi, const Real z)
{
  Real poverr;
  // 1) simple power law with index pslope
  //poverr = p0_over_r0*pow(rad/r0, pslope);
  // 2) P/rho = Cs(r,phi,z)^2/gamma
  //          = (h/r)^2*(GM1/r1+GM2/r2)/gamma
  Real x1 = rad*cos(phi);
  Real x2 = rad*sin(phi);
  Real rad1 = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(z-x3p)+SQR(rsoft));
  Real rad2 = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(z-x3s)+SQR(rsoft));
  poverr = p0_over_r0*gm0*(mu/rad1+(1.0-mu)/rad2)/gamma_gas;
  return poverr;
}


//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

static void VelProfileCyl(const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3)
{
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real rad1 = sqrt(SQR(rad)+SQR(rsoft));
  Real rad2 = sqrt(SQR(rad)+SQR(rsoft)+SQR(z));
  Real qrpl = 0.75*SQR(semia/rad1)*qrat/SQR(1.0+qrat);
  Real vel = sqrt((dslope+pslope)*p_over_r + gm0*SQR(rad/rad1)/rad1*(1.0+qrpl)
            +pslope*gm0*(1.0/rad1 -1.0/rad2));

  v1 = -vel*sin(phi);
  v2 = vel*cos(phi);
  v3 = 0.0;
  // correct with radial drift due to accretion
  if (NON_BAROTROPIC_EOS) {
	vel = 1.5*p_over_r*gamma_gas*alpha/sqrt(gm0/rad);
    v1 += -vel*cos(phi);
	v2 += -vel*sin(phi);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes ramp function R/\tau
// R(x) = a*x^2 + b*x; where a=1/(r2^2-r1r2), and b=-a*r1
// tau  = 2PI/Omega(r1)

static Real RampFunc(const Real rad, const Real phi, const Real z,
					 const Real v1, const Real v2, const Real v3)
{
  Real ramp,tau;
  if (rad >= rbuf2) {
	ramp = 10.0;
  } else if (rad >= rbuf1) {
    Real fac = 1.0/(SQR(rbuf2)-rbuf1*rbuf2);
    ramp = fac*SQR(rad)-fac*rbuf1*rad;
  } else {
	  ramp = 0.0;
  }
  tau = 2.0*PI/sqrt(gm0/rbuf1/SQR(rbuf1));
  return ramp/tau;
}

void MeshBlock::UserWorkInLoop(void)
{
  // DEBUG: binary orbit
  if (gid == 0 && pmy_mesh->time >= tstart && pmy_mesh->time <= (2.0*PI)) {
	tstart += 0.01*tbin;
    fprintf(pf,"%19.15e %19.15e %19.15e %19.15e %19.15e %19.15e %19.15e \
			    %19.15e %19.15e %19.15e\n",pmy_mesh->time, fanorm,uanorm,\
				rsep,x1s,x2s,x3s,x1p,x2p,x3p);
  }
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) { //2D
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) { //3D
    kl -= (NGHOST);
    ku += (NGHOST);
  }
  for(int k=kl; k<=ku; k++) {
    for(int j=jl; j<=ju; j++) {
      for(int i=il; i<=iu; i++) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real radp = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p));
        Real rads = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s));

        Real& u_d  = phydro->u(IDN,k,j,i);
        Real& u_m1 = phydro->u(IM1,k,j,i);
        Real& u_m2 = phydro->u(IM2,k,j,i);
        Real& u_m3 = phydro->u(IM3,k,j,i);
        Real& u_e  = phydro->u(IEN,k,j,i);

        Real& w_d  = phydro->w(IDN,k,j,i);
        Real& w_vx = phydro->w(IVX,k,j,i);
        Real& w_vy = phydro->w(IVY,k,j,i);
        Real& w_vz = phydro->w(IVZ,k,j,i);
        Real& w_p  = phydro->w(IEN,k,j,i);

		Real dt = pmy_mesh->dt;

        // apply sink cells condition within rsink( ususally >= rsoft)
        if ((radp <= rsink)||(rads <=(rsink))) {
		  u_d -= dt*u_d/tsink;
          // apply density floor, without changing momentum or energy
          u_d = (u_d > dfloor) ?  u_d : dfloor;
          w_d = u_d;
        }
		// apply wave-killing zone within [rbuf1,rbuf2] to quench m=4 mode
        Real rad, phi, z;
        GetCylCoord(pcoord,rad,phi,z,i,j,k);
		if (rad >= rbuf1) {
		  Real v1, v2, v3;
		  Real den0 = DenProfileCyl(rad,phi,z);
          VelProfileCyl(rad,phi,z,v1,v2,v3);
		  Real ramp = RampFunc(rad,phi,z,v1,v2,v3);
          u_d  -= dt*(u_d-den0)*ramp;
          u_m1 -= dt*(u_m1-den0*v1)*ramp;
          u_m2 -= dt*(u_m2-den0*v2)*ramp;
          u_m3 -= dt*(u_m3-den0*v3)*ramp;

          Real di = 1.0/u_d;
          w_vx = u_m1*di;
          w_vy = u_m2*di;
          w_vz = u_m3*di;

          if (NON_BAROTROPIC_EOS){
            Real p_over_r = PoverR(rad,phi,z);
            Real pres0 = p_over_r*den0;
            w_p -= dt*(w_p-pres0)*ramp;
            w_p = (w_p > pfloor) ?  w_p : pfloor;
            Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
            u_e = w_p/(gamma_gas-1.0)+ke;
          }
		}
		// apply extremely short cooling
        if (NON_BAROTROPIC_EOS){
		  Real pres0 = u_d*PoverR(rad,phi,z);
          w_p = pres0;
          w_p = (w_p > pfloor) ?  w_p : pfloor;
          Real di = 1.0/u_d;
          Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
          u_e = w_p/(gamma_gas-1.0)+ke;
		}
      }
    }
  }
}


void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  fclose(pf);
  return;
}


//----------------------------------------------------------------------------------------
//!\fn void compute_star_loc(const Real dto, const Real t, const Real e,
//                            const Real mf, Real *f_ptr, Real *u_ptr, Real *r_ptr)
// \brief computes the binary orbit based on its orbit elements (a,e,i) and time t.
// PURPOSE: Functions to compute the locations of binary star members
//   t:            time in units where 2PI is the orbit period
//   ecc:          eccentricity
//   mu:           mass fraction of the primary (>= 0.5)
//   x1p,x2p,x3p:  position of the primary in cartesian coord
//   x1s,x2s,x3s:  position of the secondary
//   dt:           time step since the last evaluation (need not be accurate)
//   uanorm:       the eccentric anomaly
//   fanorm:       the true anomaly
//   rsep:         the real-time separation

void solve_u(const Real t, const Real dto, const Real e, Real *u_ptr)
{
  Real du=1.0;
  Real ua = *u_ptr;
  ua += dto/(1.0-e*cos(ua));

  for (int i=1; i<=itermax; i++){
    du = - (ua-e*sin(ua)-t)/(1.0-e*cos(ua));
    ua += du;
    if (fabs(du) < dueps) {
      *u_ptr = ua;
      return;
    }
  }

  std::cout << "[solve_u error]: exceed the limit ITERMAX = " << itermax << std::endl;
  exit(1);
}

void compute_f(Real ua, const Real e, Real *f_ptr)
{
  Real asp = sqrt((1.0+e)/(1.0-e));
  ua = fmod(ua,2.0*PI);
  *f_ptr = 2.0*atan2(asp*sin(0.5*ua),cos(0.5*ua));
}

void compute_star_loc(const Real dto, const Real t, const Real e, const Real mf, Real *f_ptr, Real *u_ptr, Real *r_ptr)
{
    solve_u(t, dto, e, u_ptr);
    compute_f(*u_ptr, e, f_ptr);
    *r_ptr = (1.0 - SQR(e))/(1.0 + e*cos(*f_ptr));
}


//----------------------------------------------------------------------------------------
//! \fn void Binary(MeshBlock *pmb, const Real time, const Real dt,
//       const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief Add binary force as a user defined source term
//
void Binary(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Coordinates *pco = pmb->pcoord;
  Real src[NHYDRO];

  Real dto = dt;
  Real torb = fmod((time+0.5*dt),2.0*PI);

  //compute the binary position in its orbital plane
  //then project it to disk plane
  compute_star_loc(dto,torb,ecc,mu,&fanorm,&uanorm,&rsep);
  rsep *=semia;
  Real rcos = rsep*cos(fanorm+argp);
  Real rsin = rsep*sin(fanorm+argp);
  Real rsincos = rsin*cos(incli);
  Real rsinsin = rsin*sin(incli);
  x1p = - (1.0-mu)*rcos;
  x2p = - (1.0-mu)*rsincos;
  x3p = (1.0-mu)*rsinsin;
  x1s = mu*rcos;
  x2s = mu*rsincos;
  x3s =-mu*rsinsin;

  //compute the binary potential according to binary
  //location
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real x1  = pmb->pcoord->x1v(i);
        Real x2  = pmb->pcoord->x2v(j);
        Real x3  = pmb->pcoord->x3v(k);
        Real radp = pow(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p)+SQR(rsoft),1.5);
        //Real rads = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft*qrat),1.5);
        Real rads = pow(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s)+SQR(rsoft),1.5);
        Real srcp = dt*den*gm0*mu/radp;
        Real srcs = dt*den*gm0*(1.0-mu)/rads;
        cons(IM1,k,j,i) -= srcp*(x1-x1p)+srcs*(x1-x1s);
        cons(IM2,k,j,i) -= srcp*(x2-x2p)+srcs*(x2-x2s);
        cons(IM3,k,j,i) -= srcp*(x3-x3p)+srcs*(x3-x3s);
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
          srcp*((x1-x1p)*prim(IVX,k,j,i)+(x2-x2p)*prim(IVY,k,j,i)+(x3-x3p)*prim(IVZ,k,j,i))+
          srcs*((x1-x1s)*prim(IVX,k,j,i)+(x2-x2s)*prim(IVY,k,j,i)+(x3-x3s)*prim(IVZ,k,j,i));
  }}}

  // we delay the super-short cooling in UserWorkInLoop.
  //if (NON_BAROTROPIC_EOS && beta_th != 0.0)
  //	Cooling(pmb,time,dt,prim,bcc,cons);

}

// uncomment the following if the cooling is physical and evaluate at each substep
//
//void Cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//{
//  // apply extremely short cooling
//  //Real gam = pmb->peos->GetGamma();
//  //Real gam1 = pmb->peos->GetGamma()-1.0;
//  Real gam1 = gamma_gas-1.0;
//  Real rad, phi, z;
//  for (int k=pmb->ks; k<=pmb->ke; ++k) {
//    for (int j=pmb->js; j<=pmb->je; ++j) {
//      for (int i=pmb->is; i<=pmb->ie; ++i) {
//        GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
//		Real eth = cons(IEN,k,j,i)-0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i))+
//										SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);
//		Real eth0 = cons(IDN,k,j,i)*PoverR(rad,phi,z)/gam1;
//		Real tcool= std::max(beta_th*2.0*PI/sqrt(gm0/SQR(rad)/rad),dt);
//        cons(IEN,k,j,i) -= (eth-eth0)*dt/tcool;
//      }
//    }
//  }
//  return;
//}
