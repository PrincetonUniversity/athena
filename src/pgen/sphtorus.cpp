/*============================================================================*/
/*! \file sphtorus.cpp
 *  \brief Problem generator for the torus problem (Stone et al. 1999)
 *
 * PURPOSE: Problem generator for the torus problem (Stone et al. 1999)
/*============================================================================*/

// C/C++ headers
#include <iostream>   // endl
#include <cmath>      // sqrt
#include <cstdlib>    // srand

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

using namespace std;
// #ifdef ISOTHERMAL
// #error "Isothermal EOS cannot be used."
// #endif

/*----------------------------------------------------------------------------*/
/* function prototypes and global variables*/
void stbv_iib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on inner-x1 (left edge) of grid.
void stbv_ijb(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on inner-x2 (bottom edge) of grid.
void stbv_oib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on outer-x1 (right edge) of grid.
void stbv_ojb(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on outer-x2 (top edge) of grid.

Real A1(  Real x1,   Real x2,   Real x3);
Real A2(  Real x1,   Real x2,   Real x3);
Real A3(  Real x1,   Real x2,   Real x3);
Real magr(MeshBlock *pmb,   int i,   int j,   int k);
Real magt(MeshBlock *pmb,   int i,   int j,   int k);
Real magp(MeshBlock *pmb,   int i,   int j,   int k);
Real en, cprime, w0, rg, dist, acons, d0, amp,beta;
static Real gm,gmgas;

inline Real CUBE(Real x){
  return ((x)*(x)*(x));
}

inline Real MAX(Real x, Real y){
  return (x>y?x:y);
}

Real A1(Real x1, Real x2,Real x3)
{
  return 0.0;
}

Real A2(Real x1, Real x2,Real x3) 
{
  return 0.0;
}

Real A3(Real x1, Real x2,Real x3)
{
  Real eq29,w;
  Real a=0.0;
  Real dens;
  w=x1*sin(x2);
  eq29 = gm/(w0*(en + 1.))*(w0/x1-0.5*SQR(w0/w) - cprime);
  //cout << "eq29: "<<eq29 << endl;
  if (eq29 > 0.0) {
    dens  = pow(eq29/acons,en);
    if (dens > 100.0*d0)
	//cout<< "Inside torus dens: "<<dens<<endl;
      a = SQR(dens)/(beta);
  }
  //cout << "A3: "<<a << endl;
  return a;
}

#define ND 100

Real magr(MeshBlock *pmb,   int i,   int j,   int k)
{
  Real r,t,p,s,a,d,rd;
  Coordinates *pco = pmb->pcoord;
  int n;
  r = pco->x1f(i);
  t = pco->x2f(j);
  p = pco->x3f(k);
  s=2.0*SQR(r)*sin(t+0.5*pco->dx2f(j))*sin(0.5*pco->dx2f(j))*pco->dx3f(k);
  a=(A3(r,t+pco->dx2f(j),p+0.5*pco->dx3f(k))*sin(t+pco->dx2f(j))-A3(r,t,p+0.5*pco->dx3f(k))*sin(t))*r*pco->dx3f(k);
  return a/s;
}

Real magt(MeshBlock *pmb, int i, int j, int k)
{
  Coordinates *pco = pmb->pcoord;
  Real r,t,p,s,a,d,rd;
  int n;
  r = pco->x1f(i);
  t = pco->x2f(j);
  p = pco->x3f(k);
  s=(r+0.5*pco->dx1f(i))*pco->dx1f(i)*sin(t)*pco->dx3f(k);
  a=(A3(r+pco->dx1f(i),t,p+0.5*pco->dx3f(k))*(r+pco->dx1f(i))-A3(r,t,p+0.5*pco->dx3f(k))*r)*sin(t)*pco->dx3f(k);
  return -a/s;
}

Real magp(MeshBlock *pmb, int i, int j, int k)
{
  return 0.0;
}


//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // read parameters
  gm = pin->GetReal("problem","GM");
  gmgas = pin->GetReal("hydro","gamma");
  dist = pin->GetReal("problem","dist");
  rg = pin->GetReal("problem","rg");
  d0= pin->GetReal("problem","d0");
  amp= pin->GetReal("problem","amp");
  if (MAGNETIC_FIELDS_ENABLED)
    beta = pin->GetReal("problem","beta");
  
  w0 = 1.0;
  cprime = 0.5/dist;
  en = 1.0/(gmgas-1.0);
  acons=0.5*(dist-1.0)/dist/(en+1.0);

  // assign boundary conditions
  EnrollUserBoundaryFunction(INNER_X1, stbv_iib);
  EnrollUserBoundaryFunction(INNER_X2, stbv_ijb);
  EnrollUserBoundaryFunction(OUTER_X1, stbv_oib);
  EnrollUserBoundaryFunction(OUTER_X2, stbv_ojb);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the torus problem (Stone et al. 1999)
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  int ii,jj,ftorus;
  // Real en, cprime, w0, rg, dist, acons, d0, amp,beta
  Real ld, lm, lv, vv, rp, rm, tp, tm, tv, rv, vt, pp,eq29,dens,wt;

  AthenaArray<Real> pr;
  std::srand(gid);       

  /* allocate memory for the gas pressure */
  pr.NewAthenaArray(block_size.nx3+2*(NGHOST),block_size.nx2+2*(NGHOST),block_size.nx1+2*(NGHOST));
  
  /* Background */ 
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i)  = d0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        pr(k,j,i)=d0/pcoord->x1v(i);
      }
    }
  }

  /* Torus */
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        ld=0.0; lm=0.0;
        ftorus=0;
        vt=0.0;
        for(int jj=0;jj<10;jj++) {
          tp=pcoord->x2f(j)+pcoord->dx2f(j)*0.1*(jj+1);
          tm=pcoord->x2f(j)+pcoord->dx2f(j)*0.1*jj;
          tv=0.5*(tp+tm)+(1.0-0.5*(tp-tm)/tan(0.5*(tp-tm)))/tan(0.5*(tp+tm));
          for(int ii=0;ii<10;ii++) {
            rp = pcoord->x1f(i)+pcoord->dx1f(i)*0.1*(ii+1);
            rm = pcoord->x1f(i)+pcoord->dx1f(i)*0.1*ii;
            rv= ((SQR(SQR(rp))-SQR(SQR(rm)))/4.0)/((CUBE(rp)-CUBE(rm))/3.0);
            lv= 1.0/3.0*(CUBE(rp)-CUBE(rm))*(cos(tm)-cos(tp));
            vt+=lv;
            wt = rv*sin(tv);
            eq29 = gm/(w0*(en + 1.))*(w0/rv-0.5*SQR(w0/wt) - cprime);
            if (eq29 > 0.0) {
              dens = pow(eq29/acons,en);
              pp=dens*eq29;
              if (pp > d0/rv) {
                ftorus=1;
                ld+=dens*lv;
                lm+=dens*lv*sqrt(gm*w0)/wt;
              }
              else
                ld+=d0*lv;
            }
            else
              ld+=d0*lv;
          }
        }
        vv= 1.0/3.0*(CUBE(pcoord->x1f(i+1))-CUBE(pcoord->x1f(i)))*(cos(pcoord->x2f(j))-cos(pcoord->x2f(j+1)));
        ld/=vv; lm/=vv;
        if(ftorus==1) {
          phydro->u(IDN,k,j,i) = ld;
          phydro->u(IM3,k,j,i) = lm;
          pr(k,j,i)=MAX(acons*pow(ld,gmgas),d0/pcoord->x1v(i))*(1+amp*((double)rand()/(double)RAND_MAX-0.5));
        }
      }
    }
  }


  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          pfield->b.x1f(k,j,i) = magr(this,i,j,k);
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x2f(k,j,i) = magt(this,i,j,k);
        }
      }
    }
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }
  }

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IEN,k,j,i)=pr(k,j,i)*en+0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        // //Adding the magnetic energy contributions onto the internal energy 
        if (MAGNETIC_FIELDS_ENABLED) {
          Real bx = ((pcoord->x1f(i+1)-pcoord->x1v(i))*pfield->b.x1f(k,j,i)
             +  (pcoord->x1v(i)-pcoord->x1f(i))*pfield->b.x1f(k,j,i+1))/pcoord->dx1f(i);
          Real by = ((pcoord->x2f(j+1)-pcoord->x2v(j))*pfield->b.x2f(k,j,i)
             +  (pcoord->x2v(j)-pcoord->x2f(j))*pfield->b.x2f(k,j+1,i))/pcoord->dx2f(j);
          Real bz = (pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i))*0.5;
          phydro->u(IEN,k,j,i) += 0.5*(SQR(bx)+SQR(by)+SQR(bz));
        }
      }
    }
  }

  pr.DeleteAthenaArray();
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


/*  Boundary Condtions, outflowing, ix1, ox1, ix2, ox2  */
void stbv_iib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=1; i<=(NGHOST); i++) {
        a(IDN,k,j,is-i) = a(IDN,k,j,is);
        a(IM1,k,j,is-i) = 0.0;
        a(IM2,k,j,is-i) = a(IM2,k,j,is); //corotating ghost region
        a(IM3,k,j,is-i) = a(IM3,k,j,is);
        a(IEN,k,j,is-i)=a(IEN,k,j,is-i+1)-gm/SQR(pco->x1f(is-i+1))*a(IDN,k,j,is-i+1)*pco->dx1v(is-i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
      }
    }}
  }
  return;
}


void stbv_oib(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=1; i<=(NGHOST); i++) {
        a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
        a(IM1,k,j,ie+i) = a(IM1,k,j,ie);
        a(IM2,k,j,ie+i) = a(IM2,k,j,ie);
        a(IM3,k,j,ie+i) = a(IM3,k,j,ie);
        a(IEN,k,j,ie+i) = a(IEN,k,j,ie);
        if(a(IM1,k,j,ie+i) < 0.0)
          a(IM1,k,j,ie+i) = 0.0;
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
      }
    }}
  }
  return;
}


void stbv_ijb(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; k++) {
    for (int j=1; j<=(NGHOST); j++) {
      for (int i=is; i<=ie; i++) {
        a(IDN,k,js-j,i) = a(IDN,k,js,i);
        a(IM1,k,js-j,i) = a(IM1,k,js,i);
        a(IM2,k,js-j,i) = a(IM2,k,js,i);
        a(IM3,k,js-j,i) = a(IM3,k,js,i);
        a(IEN,k,js-j,i) = a(IEN,k,js,i);
        if(a(IM2,k,js-j,i) > 0.0)
          a(IM2,k,js-j,i) = 0.0;
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) = b.x1f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = b.x2f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) = b.x3f(k,js,i);
      }
    }}
  }
  return;
}


void stbv_ojb(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
              int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; k++) {
    for (int j=1; j<=(NGHOST); j++) {
      for (int i=is; i<=ie; i++) {
        a(IDN,k,je+j,i) = a(IDN,k,je,i);
        a(IM1,k,je+j,i) = a(IM1,k,je,i);
        a(IM2,k,je+j,i) = a(IM2,k,je,i);
        a(IM3,k,je+j,i) = a(IM3,k,je,i);
        a(IEN,k,je+j,i) = a(IEN,k,je,i);
        if(a(IM2,k,je+j,i) < 0.0)
          a(IM2,k,je+j,i) = 0.0;
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) = b.x1f(k,(je  ),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) = b.x3f(k,(je  ),i);
      }
    }}
  }
  return;
}

