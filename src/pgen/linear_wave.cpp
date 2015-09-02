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
#include "../mesh.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // EOS
#include "../field/field.hpp"      // Field
#include "../coordinates/coordinates.hpp" // Coordinates

//======================================================================================
//! \file linear_wave.c
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//
// In 1D, the problem is setup along one of the three coordinate axes (specified by
// setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
// automatically sets the wavevector along the domain diagonal.
//======================================================================================

// Parameters which define initial solution -- made global so that they can be shared
// with functions A1,2,3 which compute vector potentials
static Real bx0, by0, bz0, dby, dbz;
static int wave_flag;
static Real ang_2, ang_3; // Rotation angles about the y and z' axis
static Real sin_a2, cos_a2, sin_a3, cos_a3;
static Real amp, lambda, k_par; // amplitude, Wavelength, 2*PI/wavelength
static Real gm1,iso_cs;

// functions to compute vector potential to initialize the solution
static Real A1(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A3(const Real x1, const Real x2, const Real x3);

// function to compute eigenvectors of linear waves
static void Eigensystem(const Real d, const Real v1, const Real v2, const Real v3,
  const Real h, const Real b1, const Real b2, const Real b3, const Real x, const Real y,
  Real eigenvalues[(NWAVE)],
  Real right_eigenmatrix[(NWAVE)][(NWAVE)], Real left_eigenmatrix[(NWAVE)][(NWAVE)]);

//======================================================================================
//! \fn ProblemGenerator
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  Coordinates *pco = pmb->pcoord;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  gm1 = (pfl->pf_eos->GetGamma() - 1.0);
  iso_cs = pfl->pf_eos->GetIsoSoundSpeed();

// Read initial conditions

  wave_flag = pin->GetInteger("problem","wave_flag");
  amp = pin->GetReal("problem","amp");
  Real vflow = pin->GetOrAddReal("problem","vflow",0.0);
  ang_2 = pin->GetOrAddReal("problem","ang_2",-999.9);
  ang_3 = pin->GetOrAddReal("problem","ang_3",-999.9);

  Real x1size = pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min;
  Real x2size = pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min;
  Real x3size = pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min;

// For wavevector along coordinate axes, set desired values of ang_2/ang_3.
//    For example, for 1D problem use ang_2 = ang_3 = 0.0
//    For wavevector along grid diagonal, do not input values for ang_2/ang_3.
// Code below will automatically calculate these imposing periodicity and exactly one
// wavelength along each grid direction

// User should never input -999.9 in angles
  if (ang_3 == -999.9) ang_3 = atan(x1size/x2size);
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);

  if (ang_2 == -999.9) ang_2 = atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
  sin_a2 = sin(ang_2);
  cos_a2 = cos(ang_2);

  Real x1 = x1size*cos_a2*cos_a3;
  Real x2 = x2size*cos_a2*sin_a3;
  Real x3 = x3size*sin_a2;

// For lambda choose the smaller of the 3
  lambda = x1;
  if (pmb->block_size.nx2 > 1 && ang_3 != 0.0) lambda = std::min(lambda,x2);
  if (pmb->block_size.nx3 > 1 && ang_2 != 0.0) lambda = std::min(lambda,x3);

// Initialize k_parallel
  k_par = 2.0*(PI)/lambda;

// Get eigenvectors, where the quantities u0 and bx0 are parallel to the wavevector, and
// v0,w0,by0,bz0 are perpendicular.

  Real d0 = 1.0;
  Real p0 = 0.0;
  Real u0 = vflow;
  Real v0 = 0.0;
  Real w0 = 0.0;
  bx0 = 1.0;
  by0 = sqrt(2.0);
  bz0 = 0.5;
  Real xfact = 0.0;
  Real yfact = 1.0;
  Real h0 = 0.0;

  if (NON_BAROTROPIC_EOS) {
    p0 = 1.0/(pfl->pf_eos->GetGamma());
    h0 = ((p0/gm1 + 0.5*d0*(u0*u0+v0*v0+w0*w0)) + p0)/d0;
    if (MAGNETIC_FIELDS_ENABLED) h0 += (bx0*bx0+by0*by0+bz0*bz0)/d0;
  }

  Real ev[NWAVE], rem[NWAVE][NWAVE], lem[NWAVE][NWAVE];
  Eigensystem(d0,u0,v0,w0,h0,bx0,by0,bz0,xfact,yfact,ev,rem,lem);

// Initialize the magnetic fields

  if (MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> a1,a2,a3;
    int nx1 = (ie-is)+1 + 2*(NGHOST);
    int nx2 = (je-js)+1 + 2*(NGHOST);
    int nx3 = (ke-ks)+1 + 2*(NGHOST);
    a1.NewAthenaArray(nx3,nx2,nx1);
    a2.NewAthenaArray(nx3,nx2,nx1);
    a3.NewAthenaArray(nx3,nx2,nx1);

    // wave amplitudes
    dby = amp*rem[NWAVE-2][wave_flag];
    dbz = amp*rem[NWAVE-1][wave_flag];

    int level=pmb->uid.GetLevel();
    // Initialize components of the vector potential
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          if((pmb->nblevel[1][0][1]>level && j==js) || (pmb->nblevel[1][2][1]>level && j==je+1)
          || (pmb->nblevel[0][1][1]>level && k==ks) || (pmb->nblevel[2][1][1]>level && k==ke+1)
          || (pmb->nblevel[0][0][1]>level && j==js && k==ks)   || (pmb->nblevel[0][2][1]>level && j==je+1 && k==ks)
          || (pmb->nblevel[2][0][1]>level && j==js && k==ke+1) || (pmb->nblevel[2][2][1]>level && j==je+1 && k==ke+1))
            a1(k,j,i) = 0.5*(A1(pco->x1f(i)+0.25*pco->dx1f(i), pco->x2f(j), pco->x3f(k))+
                             A1(pco->x1f(i)+0.75*pco->dx1f(i), pco->x2f(j), pco->x3f(k)));
          else a1(k,j,i) = A1(0.5*(pco->x1f(i)+pco->x1f(i+1)), pco->x2f(j), pco->x3f(k));
          if((pmb->nblevel[1][1][0]>level && i==is) || (pmb->nblevel[1][1][2]>level && i==ie+1)
          || (pmb->nblevel[0][1][1]>level && k==ks) || (pmb->nblevel[2][1][1]>level && k==ke+1)
          || (pmb->nblevel[0][1][0]>level && i==is && k==ks)   || (pmb->nblevel[0][1][2]>level && i==ie+1 && k==ks)
          || (pmb->nblevel[2][1][0]>level && i==is && k==ke+1) || (pmb->nblevel[2][1][2]>level && i==ie+1 && k==ke+1))
            a2(k,j,i) = 0.5*(A2(pco->x1f(i), pco->x2f(j)+0.25*pco->dx2f(j), pco->x3f(k))+
                             A2(pco->x1f(i), pco->x2f(j)+0.75*pco->dx2f(j), pco->x3f(k)));
          else a2(k,j,i) = A2(pco->x1f(i), 0.5*(pco->x2f(j)+pco->x2f(j+1)), pco->x3f(k));
          if((pmb->nblevel[1][1][0]>level && i==is) || (pmb->nblevel[1][1][2]>level && i==ie+1)
          || (pmb->nblevel[1][0][1]>level && j==js) || (pmb->nblevel[1][2][1]>level && j==je+1)
          || (pmb->nblevel[1][0][0]>level && i==is && j==js)   || (pmb->nblevel[1][0][2]>level && i==ie+1 && j==js)
          || (pmb->nblevel[1][2][0]>level && i==is && j==je+1) || (pmb->nblevel[1][2][2]>level && i==ie+1 && j==je+1))
            a3(k,j,i) = 0.5*(A3(pco->x1f(i), pco->x2f(j), pco->x3f(k)+0.25*pco->dx3f(k))+
                             A3(pco->x1f(i), pco->x2f(j), pco->x3f(k)+0.75*pco->dx3f(k)));
          else a3(k,j,i) = A3(pco->x1f(i), pco->x2f(j), 0.5*(pco->x3f(k)+pco->x3f(k+1)));
        }
      }
    }

    // Initialize interface fields
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfd->b.x1f(k,j,i) = (a3(k  ,j+1,i) - a3(k,j,i))/pco->dx2f(j) -
                            (a2(k+1,j  ,i) - a2(k,j,i))/pco->dx3f(k);
      }
    }}

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        pfd->b.x2f(k,j,i) = (a1(k+1,j,i  ) - a1(k,j,i))/pco->dx3f(k) -
                            (a3(k  ,j,i+1) - a3(k,j,i))/pco->dx1f(i);
      }
    }}

    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
       pfd->b.x3f(k,j,i) = (a2(k,j  ,i+1) - a2(k,j,i))/pco->dx1f(i) -
                           (a1(k,j+1,i  ) - a1(k,j,i))/pco->dx2f(j);
      }
    }}
    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
  }

// initialize conserved variables

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real x = cos_a2*(pco->x1v(i)*cos_a3 + pco->x2v(j)*sin_a3) + pco->x3v(k)*sin_a2;
      Real sn = sin(k_par*x);

      pfl->u(IDN,k,j,i) = d0 + amp*sn*rem[0][wave_flag];

      Real mx = d0*vflow + amp*sn*rem[1][wave_flag];
      Real my = amp*sn*rem[2][wave_flag];
      Real mz = amp*sn*rem[3][wave_flag];

      pfl->u(IM1,k,j,i) = mx*cos_a2*cos_a3 - my*sin_a3 - mz*sin_a2*cos_a3;
      pfl->u(IM2,k,j,i) = mx*cos_a2*sin_a3 + my*cos_a3 - mz*sin_a2*sin_a3;
      pfl->u(IM3,k,j,i) = mx*sin_a2                    + mz*cos_a2;

      if (NON_BAROTROPIC_EOS) {
        pfl->u(IEN,k,j,i) = p0/gm1 + 0.5*d0*u0*u0 + amp*sn*rem[4][wave_flag];
        if (MAGNETIC_FIELDS_ENABLED) {
          pfl->u(IEN,k,j,i) += 0.5*(bx0*bx0+by0*by0+bz0*bz0);
        }
      }
    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A1(const Real x1,const Real x2,const Real x3)
//  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0, and Ay,
//  Az are functions of x and y alone.

static Real A1(const Real x1, const Real x2, const Real x3)
{
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay =  bz0*x - (dbz/k_par)*cos(k_par*(x));
  Real Az = -by0*x + (dby/k_par)*cos(k_par*(x)) + bx0*y;

  return -Ay*sin_a3 - Az*sin_a2*cos_a3;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A2(const Real x1,const Real x2,const Real x3)
//  \brief A2: 2-component of vector potential

static Real A2(const Real x1, const Real x2, const Real x3)
{
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay =  bz0*x - (dbz/k_par)*cos(k_par*(x));
  Real Az = -by0*x + (dby/k_par)*cos(k_par*(x)) + bx0*y;

  return Ay*cos_a3 - Az*sin_a2*sin_a3;
}

//--------------------------------------------------------------------------------------
//! \fn static Real A3(const Real x1,const Real x2,const Real x3)
//  \brief A3: 3-component of vector potential

static Real A3(const Real x1, const Real x2, const Real x3)
{
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Az = -by0*x + (dby/k_par)*cos(k_par*(x)) + bx0*y;

  return Az*cos_a2;
}

//--------------------------------------------------------------------------------------
//! \fn static void Eigensystem()
//  \brief computes eigenvectors of linear waves 

static void Eigensystem(const Real d, const Real v1, const Real v2, const Real v3,
  const Real h, const Real b1, const Real b2, const Real b3, const Real x, const Real y,
  Real eigenvalues[(NWAVE)],
  Real right_eigenmatrix[(NWAVE)][(NWAVE)], Real left_eigenmatrix[(NWAVE)][(NWAVE)])
{
  if (MAGNETIC_FIELDS_ENABLED) {

//--- Adiabatic MHD ---

    if (NON_BAROTROPIC_EOS) {
      Real vsq,btsq,bt_starsq,vaxsq,hp,twid_asq,cfsq,cf,cssq,cs;
      Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,vbet,alpha_f,alpha_s;
      Real isqrtd,sqrtd,s,twid_a,qf,qs,af_prime,as_prime,afpbb,aspbb,vax;
      Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
      Real ct2,tsum,tdif,cf2_cs2;
      Real qa,qb,qc,qd;
      vsq = v1*v1 + v2*v2 + v3*v3;
      btsq = b2*b2 + b3*b3;
      bt_starsq = (gm1 - (gm1 - 1.0)*y)*btsq;
      vaxsq = b1*b1/d;
      hp = h - (vaxsq + btsq/d);
      twid_asq = std::max((gm1*(hp-0.5*vsq)-(gm1-1.0)*x), TINY_NUMBER);

      // Compute fast- and slow-magnetosonic speeds (eq. B18)
      ct2 = bt_starsq/d;
      tsum = vaxsq + ct2 + twid_asq;
      tdif = vaxsq + ct2 - twid_asq;
      cf2_cs2 = sqrt((double)(tdif*tdif + 4.0*twid_asq*ct2));

      cfsq = 0.5*(tsum + cf2_cs2);
      cf = sqrt((double)cfsq);

      cssq = twid_asq*vaxsq/cfsq;
      cs = sqrt((double)cssq);

      // Compute beta(s) (eqs. A17, B20, B28)
      bt = sqrt(btsq);
      bt_star = sqrt(bt_starsq);
      if (bt == 0.0) {
        bet2 = 1.0;
        bet3 = 0.0;
      } else {
        bet2 = b2/bt;
        bet3 = b3/bt;
      }
      bet2_star = bet2/sqrt(gm1 - (gm1-1.0)*y);
      bet3_star = bet3/sqrt(gm1 - (gm1-1.0)*y);
      bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
      vbet = v2*bet2_star + v3*bet3_star;

      // Compute alpha(s) (eq. A16)
      if ((cfsq-cssq) == 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else if ( (twid_asq - cssq) <= 0.0) {
        alpha_f = 0.0;
        alpha_s = 1.0;
      } else if ( (cfsq - twid_asq) <= 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else {
        alpha_f = sqrt((twid_asq - cssq)/(cfsq - cssq));
        alpha_s = sqrt((cfsq - twid_asq)/(cfsq - cssq));
      }

      // Compute Q(s) and A(s) (eq. A14-15), etc.
      sqrtd = sqrt(d);
      isqrtd = 1.0/sqrtd;
      s = SIGN(b1);
      twid_a = sqrt(twid_asq);
      qf = cf*alpha_f*s;
      qs = cs*alpha_s*s;
      af_prime = twid_a*alpha_f*isqrtd;
      as_prime = twid_a*alpha_s*isqrtd;
      afpbb = af_prime*bt_star*bet_starsq;
      aspbb = as_prime*bt_star*bet_starsq;

      // Compute eigenvalues (eq. B17)
      vax = sqrt(vaxsq);
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

      qa = alpha_f*v2;
      qb = alpha_s*v2;
      qc = qs*bet2_star;
      qd = qf*bet2_star;
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
      norm = 0.5/twid_asq;
      cff = norm*alpha_f*cf;
      css = norm*alpha_s*cs;
      qf *= norm;
      qs *= norm;
      af = norm*af_prime*d;
      as = norm*as_prime*d;
      afpb = norm*af_prime*bt_star;
      aspb = norm*as_prime*bt_star;
  
      // Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f}
      norm *= gm1;
      alpha_f *= norm;
      alpha_s *= norm;
      q2_star = bet2_star/bet_starsq;
      q3_star = bet3_star/bet_starsq;
      vqstr = (v2*q2_star + v3*q3_star);
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

//--- Isothermal MHD ---

    } else {
      Real btsq,bt_starsq,vaxsq,twid_csq,cfsq,cf,cssq,cs;
      Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,alpha_f,alpha_s;
      Real sqrtd,s,twid_c,qf,qs,af_prime,as_prime,vax;
      Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
      Real ct2,tsum,tdif,cf2_cs2;
      Real di = 1.0/d;
      btsq = b2*b2 + b3*b3;
      bt_starsq = btsq*y;
      vaxsq = b1*b1*di;
      twid_csq = (iso_cs*iso_cs) + x;

      // Compute fast- and slow-magnetosonic speeds (eq. B39)
      ct2 = bt_starsq*di;
      tsum = vaxsq + ct2 + twid_csq;
      tdif = vaxsq + ct2 - twid_csq;
      cf2_cs2 = sqrt((double)(tdif*tdif + 4.0*twid_csq*ct2));
    
      cfsq = 0.5*(tsum + cf2_cs2);
      cf = sqrt((double)cfsq);
    
      cssq = twid_csq*vaxsq/cfsq;
      cs = sqrt((double)cssq);
  
      // Compute beta's (eqs. A17, B28, B40)
      bt = sqrt(btsq);
      bt_star = sqrt(bt_starsq);
      if (bt == 0.0) {
        bet2 = 1.0;
        bet3 = 0.0;
      } 
      else {
        bet2 = b2/bt;
        bet3 = b3/bt;
      }
      bet2_star = bet2/sqrt(y);
      bet3_star = bet3/sqrt(y);
      bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;

      // Compute alpha's (eq. A16)
      if ((cfsq-cssq) == 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else if ((twid_csq - cssq) <= 0.0) {
        alpha_f = 0.0;
        alpha_s = 1.0;
      } else if ((cfsq - twid_csq) <= 0.0) {
        alpha_f = 1.0;
        alpha_s = 0.0;
      } else {
        alpha_f = sqrt((twid_csq - cssq)/(cfsq - cssq));
        alpha_s = sqrt((cfsq - twid_csq)/(cfsq - cssq));
      }

      // Compute Q's (eq. A14-15), etc.
      sqrtd = sqrt(d);
      s = SIGN(b1);
      twid_c = sqrt(twid_csq);
      qf = cf*alpha_f*s;
      qs = cs*alpha_s*s;
      af_prime = twid_c*alpha_f/sqrtd;
      as_prime = twid_c*alpha_s/sqrtd;

      // Compute eigenvalues (eq. B38)
      vax  = sqrt(vaxsq);
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
      norm = 0.5/twid_csq;
      cff = norm*alpha_f*cf;
      css = norm*alpha_s*cs;
      qf *= norm;
      qs *= norm;
      af = norm*af_prime*d;
      as = norm*as_prime*d;
      afpb = norm*af_prime*bt_star;
      aspb = norm*as_prime*bt_star;

      q2_star = bet2_star/bet_starsq;
      q3_star = bet3_star/bet_starsq;
      vqstr = (v2*q2_star + v3*q3_star);
  
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
  } else {

//--- Adiabatic Hydrodynamics ---

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

//--- Isothermal Hydrodynamics ---

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
}
