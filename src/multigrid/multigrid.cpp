//========================================================================================
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C/C++ headers
#include <iostream>
#include <cmath>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memset

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "./multigrid.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
{
  AthenaArray<Real> &dst=u_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int nsrc=ns+n;
#pragma ivdep
    for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
#pragma ivdep
      for(int j=ngh, mj=js; mj<=je; j++, mj++) {
#pragma ivdep
        for(int i=ngh, mi=is; mi<=ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                  Real fac)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac)
{
  AthenaArray<Real> &dst=src_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int nsrc=ns+n;
#pragma ivdep
    for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
#pragma ivdep
      for(int j=ngh, mj=js; mj<=je; j++, mj++) {
#pragma ivdep
        for(int i=ngh, mi=is; mi<=ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i)*fac;
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the result, including the ghost zone
void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
{
  const AthenaArray<Real> &src=src_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=0;
  ie=nx_+2*ngh_-1, je=ny_+2*ngh_-1, ke=nz_+2*ngh_-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int ndst=ns+n;
#pragma ivdep
    for(int k=ngh-ngh_, mk=ks; mk<=ke; k++, mk++) {
#pragma ivdep
      for(int j=ngh-ngh_, mj=js; mj<=je; j++, mj++) {
#pragma ivdep
        for(int i=ngh-ngh_, mi=is; mi<=ie; i++, mi++)
          dst(ndst,k,j,i)=src(n,mk,mj,mi);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::ZeroClearData(int lev)
//  \brief Clear the data array with zero
void Multigrid::ZeroClearData(int lev)
{
  AthenaArray<Real> &u=u_[lev];
  std::memset(u.data(), 0, u.GetSizeInBytes());
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::Restrict(int clev)
//  \brief Restrict the potential/density to level=clev
void Multigrid::Restrict(int clev)
{
  AthenaArray<Real> &dst=src_[clev];
  AthenaArray<Real> &src=def_[clev+1];
  int ll=nlev_-1-clev;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
        for(int i=is, fi=is; i<=ie; i++, fi+=2)
          dst(n, k, j, i)=0.125*(src(n, fk,   fj,   fi)+src(n, fk,   fj,   fi+1)
                                +src(n, fk,   fj+1, fi)+src(n, fk,   fj+1, fi+1)
                                +src(n, fk+1, fj,   fi)+src(n, fk+1, fj,   fi+1)
                                +src(n, fk+1, fj+1, fi)+src(n, fk+1, fj+1, fi+1));
      }
    }
  }
  return;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::ProlongateAndCorrect(int clev)
//  \brief Prolongate the potential from level=clev using tri-linear interpolation
void Multigrid::ProlongateAndCorrect(int clev)
{
  const AthenaArray<Real> &src=u_[clev];
  AthenaArray<Real> &dst=def_[clev+1];
  int ll=nlev_-1-clev;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
        for(int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj  ,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk  ,fj  ,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk+1,fj+1,fi  )+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i+1)+src(n,k,j+1,i+1)));
          dst(n,fk+1,fj+1,fi+1)+=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i+1)+src(n,k,j+1,i+1)));
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::FMGProlongate(int clev)
//  \brief Prolongate the potential from level=clev for Full Multigrid cycle
void Multigrid::FMGProlongate(int clev)
{
  AthenaArray<Real> &src=u_[clev];
  AthenaArray<Real> &dst=u_[clev+1];
  int ll=nlev_-1-clev;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
        for(int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj  ,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk  ,fj  ,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk+1,fj+1,fi  )=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i-1)
                          +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                          +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i+1)+src(n,k,j+1,i+1)));
          dst(n,fk+1,fj+1,fi+1)=0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i+1)
                          +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                          +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i+1)+src(n,k,j+1,i+1)));
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(int n, int nrm)
{
  AthenaArray<Real> &def=def_[nlev_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

  Real norm=0.0;
  if(nrm==0) { // special case: max norm
#pragma ivdep
    for(int k=ks; k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          norm=std::max(norm,std::fabs(def(n,k,j,i)));
      }
    }
  }
  else if (nrm==1) {
#pragma ivdep
    for(int k=ks; k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          norm+=std::fabs(def(n,k,j,i));
      }
    }
  }
  else { // nrm>1 -> nrm=2
#pragma ivdep
    for(int k=ks; k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          norm+=SQR(def(n,k,j,i));
      }
    }
    norm=std::sqrt(norm);
  }
  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotalSource(int n)
//  \brief calculate the sum of the source function

Real CalculateTotalSource(int n)
{
  AthenaArray<Real> &src=src_[nlev_-1];
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

#pragma ivdep
  for(int k=ks; k<=ke; k++) {
#pragma ivdep
    for(int j=js; j<=je; j++) {
#pragma ivdep
      for(int i=is; i<=ie; i++)
        s+=src(n,k,j,i);
    }
  }
  return s*dx_*dx_*dx_;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractMeanSource(int n, Real ave)
//  \brief subtract the mean value of the source function for periodic boundary cases

void Multigrid::SubtractMeanSource(int n, Real ave)
{
  AthenaArray<Real> &src=src_[nlev_-1];
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

#pragma ivdep
  for(int k=ks; k<=ke; k++) {
#pragma ivdep
    for(int j=js; j<=je; j++) {
#pragma ivdep
      for(int i=is; i<=ie; i++)
        src(n,k,j,i)-=ave;
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::Multigrid(int invar, int nx, int ny, int nz, int ngh, Real dx)
//  \brief Multigrid constructor

void Multigrid::Multigrid(int invar, int nx, int ny, int nz, int ngh, Real dx)
{
  nvar=invar, ngh_=ngh;
  nx_=nx, ny_=ny, nz_=nz;
  dx=dx_;

  nlev_=0;
  int n = std::min(nx,std::min(ny, nz));
  for(int l=0; l<20; l++) {
    if((1<<l) == n) {
      nlev_=l+1;
      break;
    }
  }

   // allocate arrays
  u_ = new AthenaArray<Real>[nlev_];
  src_ = new AthenaArray<Real>[nlev_];
  def_ = new AthenaArray<Real>[nlev_];
  for(int l=nlev_-1; l>=0; l++) {
    int ll=nlev_-1-l;
    int ncx=(nx>>ll)+2*ngh, ncy=(ny>>ll)+2*ngh, ncz=(nz>>ll)+2*ngh;
    u_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
    src_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
  }
}

//----------------------------------------------------------------------------------------
//! \fn  virtual void Multigrid::~Multigrid
//  \brief Multigrid destroctor

virtual void Multigrid::~Multigrid()
{
  for(int l=0; l<nlev_; l++) {
    u_[l].DeleteAthenaArray();
    src_[l].DeleteAthenaArray();
    def_[l].DeleteAthenaArray();
  }
  delete [] u_;
  delete [] src_;
  delete [] def_;
}

