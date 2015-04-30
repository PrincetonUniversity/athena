/*============================================================================*/
/*! \file sphtorus.cpp
 *  \brief Problem generator for the torus problem (Stone et al. 1999)
 *
 * PURPOSE: Problem generator for the torus problem (Stone et al. 1999)
/*============================================================================*/
// Primary header
#include "../fluid/fluid.hpp"

// C++ headers
#include <iostream>   // endl
#include <cmath>      // sqrt
#include <cstdlib>    // srand

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../coordinates/coordinates.hpp"  // PrimToCons()
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/eos/eos.hpp"    // ParameterInput
#include "../bvals/bvals.hpp" // EnrollFluidBValFunction


#ifdef ISOTHERMAL
#error "Isothermal EOS cannot be used."
#endif

static Real gm,gmgas;


inline Real CUBE(Real x){
  return ((x)*(x)*(x));
}

inline Real MAX(Real x, Real y){
	return (x>y?x:y);
}

/*----------------------------------------------------------------------------*/
/* function prototypes and global variables*/
void stbv_iib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on inner-x1 (left edge) of grid.
void stbv_ijb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on inner-x2 (bottom edge) of grid.
void stbv_oib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on outer-x1 (right edge) of grid.
void stbv_ojb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on outer-x2 (top edge) of grid.

//problem generator
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  int i, is = pmb->is, ie = pmb->ie;
  int j, js = pmb->js, je = pmb->je;
  int k, ks = pmb->ks, ke = pmb->ke;
  int ii,jj,ftorus;
  Real en, cprime, w0, rg, dist, acons, d0, amp;
  Real ld, lm, lv, vv, rp, rm, tp, tm, tv, rv, vt, pp,eq29,dens,wt;
  AthenaArray<Real> pr;

  /* initialize a random seed */
  std::srand(0);             
  
  /* read parameters */
  gm = pin->GetReal("problem","GM");
  gmgas = pin->GetReal("fluid","gamma");
  dist = pin->GetReal("problem","dist");
  rg = pin->GetReal("problem","rg");
  d0= pin->GetReal("problem","d0");
  amp= pin->GetReal("problem","amp");
  
  /* calculate some constant values */
  w0 = 1.0;
  cprime = 0.5/dist;
  en = 1.0/(gmgas-1.0);
  acons=0.5*(dist-1.0)/dist/(en+1.0);

  /* assign boundary conditions and gravitational force function*/
  pmb->pbval->EnrollFluidBoundaryFunction(inner_x1, stbv_iib);
  pmb->pbval->EnrollFluidBoundaryFunction(inner_x2, stbv_ijb);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x1, stbv_oib);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x2, stbv_ojb);

  /* allocate memory for the gas pressure */
  pr.NewAthenaArray(pmb->block_size.nx3+2*(NGHOST),pmb->block_size.nx2+2*(NGHOST),pmb->block_size.nx1+2*(NGHOST));
  /* Background */ 
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pfl->u(IDN,k,j,i)  = d0;
        pfl->u(IM1,k,j,i) = 0.0;
        pfl->u(IM2,k,j,i) = 0.0;
        pfl->u(IM3,k,j,i) = 0.0;
        pr(k,j,i)=d0/pmb->x1v(i);
      }
    }
  }

  /* Torus */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ld=0.0; lm=0.0;
        ftorus=0;
        vt=0.0;
        for(jj=0;jj<10;jj++) {
          tp=pmb->x2f(j)+pmb->dx2f(j)*0.1*(jj+1);
          tm=pmb->x2f(j)+pmb->dx2f(j)*0.1*jj;
          tv=0.5*(tp+tm)+(1.0-0.5*(tp-tm)/tan(0.5*(tp-tm)))/tan(0.5*(tp+tm));
          for(ii=0;ii<10;ii++) {
            rp = pmb->x1f(i)+pmb->dx1f(i)*0.1*(ii+1);
            rm = pmb->x1f(i)+pmb->dx1f(i)*0.1*ii;
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
        if(ftorus==1)
        {
          vv= 1.0/3.0*(CUBE(pmb->x1f(i+1))-CUBE(pmb->x1f(i)))*(cos(pmb->x2f(j))-cos(pmb->x2f(j+1)));
          ld/=vv; lm/=vv;
          pfl->u(IDN,k,j,i) = ld;
          pfl->u(IM3,k,j,i) = lm;
          pr(k,j,i)=MAX(acons*pow(ld,gmgas),d0/pmb->x1v(i))*(1+amp*((double)rand()/(double)RAND_MAX-0.5));
        }
        pfl->u(IEN,k,j,i)=pr(k,j,i)*en+0.5*SQR(pfl->u(IM3,k,j,i))/pfl->u(IDN,k,j,i);
      }
    }
  }

  pr.DeleteAthenaArray();
}


/*  Boundary Condtions, outflowing, ix1, ox1, ix2, ox2  */
void stbv_iib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
  Real pg;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=(NGHOST); i++) {
        a(IDN,k,j,is-i) = a(IDN,k,j,is);
        if(a(IM1,k,j,is-i) > 0.0)
        {
          a(IM1,k,j,is-i) = 0.0;
        }
//        a(IM1,k,j,is-i) = 0.0;
        a(IM2,k,j,is-i) = 0.0;
        a(IM3,k,j,is-i) = 0.0;
        pg = (a(IEN,k,j,is-i+1)-0.5*(SQR(a(IM1,k,j,is-i+1))+SQR(a(IM2,k,j,is-i+1))+SQR(a(IM3,k,j,is-i+1)))/a(IDN,k,j,is-i+1))*(gmgas-1.0);
        pg-=gm/SQR(pmb->x1f(is-i+1))*a(IDN,k,j,is-i+1)*pmb->dx1v(is-i);
        a(IEN,k,j,is-i)=pg/(gmgas-1.0)+0.5*SQR(a(IM1,k,j,is-i))/a(IDN,k,j,is-i);
      }
    }
  }

  return;
}


void stbv_oib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke))
{
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=(NGHOST); i++) {
        a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
        a(IM2,k,j,ie+i) = a(IM2,k,j,ie);
        a(IM3,k,j,ie+i) = a(IM3,k,j,ie);
        if(a(IM1,k,j,ie+i) < 0.0)
        {
          a(IEN,k,j,ie+i) -= 0.5*SQR(a(IM1,k,j,ie+i))/a(IDN,k,j,ie+i);
          a(IM1,k,j,ie+i) = 0.0;
        }
      }
    }
  }

  return;
}


void stbv_ijb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=(NGHOST); j++) {
      for (i=is; i<=ie; i++) {
        a(IDN,k,js-j,i) = a(IDN,k,js,i);
        a(IM1,k,js-j,i) = a(IM1,k,js,i);
        a(IM3,k,js-j,i) = a(IM3,k,js,i);
        if(a(IM2,k,js-j,i) > 0.0)
        {
          a(IEN,k,js-j,i) -= 0.5*SQR(a(IM2,k,js-j,i))/a(IDN,k,js-j,i);
          a(IM2,k,js-j,i) = 0.0;
        }
      }
    }
  }
  return;
}


void stbv_ojb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=(NGHOST); j++) {
      for (i=is; i<=ie; i++) {
        a(IDN,k,je+j,i) = a(IDN,k,je,i);
        a(IM1,k,je+j,i) = a(IM1,k,je,i);
        a(IM3,k,je+j,i) = a(IM3,k,je,i);
        if(a(IM2,k,je+j,i) < 0.0)
        {
          a(IEN,k,je+j,i) -= 0.5*SQR(a(IM2,k,je+j,i))/a(IDN,k,je+j,i);
          a(IM2,k,je+j,i) = 0.0;
        }
      }
    }
  }
  return;
}

