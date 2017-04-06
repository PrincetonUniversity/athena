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
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
{
  AthenaArray<Real> &dst=u_[nlevel_-1];
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
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                  Real fac)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac)
{
  AthenaArray<Real> &dst=fmgsrc_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
  if(fac==1.0) {
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
  }
  else {
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
  }
  current_level_=nlevel_-1;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictFMGSource(void)
//  \brief restrict the source through all the multigrid levels

void Multigrid::RestrictFMGSource(void)
{
  current_level_=nlevel_-1;
  std::memcpy(fmgsrc_[current_level_].data(), src_[current_level_].data(),
              src_[current_level_].GetSizeInBytes());
  for(; current_level_>0; current_level_--) {
    int ll=nlevel_-current_level_;
    ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
    AthenaArray<Real> &csrc=fmgsrc_[current_level_-1];
    AthenaArray<Real> &fsrc=fmgsrc_[current_level_];
#pragma ivdep
    for(int n=0; n<nvar; n++) {
#pragma ivdep
      for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
#pragma ivdep
        for(int j=js, fj=js; j<=je; j++, fj+=2) {
#pragma ivdep
          for(int i=is, fi=is; i<=ie; i++, fi+=2)
            csrc(n, k, j, i)=0.125*(fsrc(n, fk,   fj,   fi)+fsrc(n, fk,   fj,   fi+1)
                                   +fsrc(n, fk,   fj+1, fi)+fsrc(n, fk,   fj+1, fi+1)
                                   +fsrc(n, fk+1, fj,   fi)+fsrc(n, fk+1, fj,   fi+1)
                                   +fsrc(n, fk+1, fj+1, fi)+fsrc(n, fk+1, fj+1, fi+1));
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the result, including the ghost zone
void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
{
  const AthenaArray<Real> &src=src_[nlevel_-1];
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
//! \fn void Multigrid::ZeroClearData(void)
//  \brief Clear the data array with zero
void Multigrid::ZeroClearData(void)
{
  AthenaArray<Real> &u=u_[current_level_];
  std::memset(u.data(), 0, u.GetSizeInBytes());
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ApplyPhysicalBoundaries(void)
//  \brief A

void Multigrid::ApplyPhysicalBoundaries(void)
{
  AthenaArray<Real> &dst=u_[current_level_];
  int ll=nlevel_-1-current_level_;
  int ncx=nx_>>ll, ncy=ny_>>ll, ncz=nz_>>ll;
  int is=ngh_, ie=ncx_+ngh_-1, js=ngh_, je=ncy_+ngh_-1, ks=ngh_, ke=ncz_+ngh_-1;
  int bis=is-ngh_, bie=ie+ngh_, bjs=js, bje=je, bks=ks, bke=ke;
  Real dx, dy, dz;
  Real dx=rdx_/(Real)(1<<ll), dy=rdy_/(Real)(1<<ll), dz=rdz_/(Real)(1<<ll);
  Real x0=size.x1min-((Real)ngh_+0.5)*dx;
  Real y0=size.x2min-((Real)ngh_+0.5)*dy;
  Real z0=size.x3min-((Real)ngh_+0.5)*dz;
  if(MGBoundaryFunction_[INNER_X2]==NULL) bjs=js-ngh_;
  if(MGBoundaryFunction_[OUTER_X2]==NULL) bje=je+ngh_;
  if(MGBoundaryFunction_[INNER_X3]==NULL) bks=ks-ngh_;
  if(MGBoundaryFunction_[OUTER_X3]==NULL) bke=ke*ngh_;

  // Apply boundary function on inner-x1
  if (MGBoundaryFunction_[INNER_X1] != NULL)
    MGBoundaryFunction_[INNER_X1](dst, time, dt, nvar, is, ie, bjs, bje, bks, bke,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x1
  if (MGBoundaryFunction_[OUTER_X1] != NULL)
    MGBoundaryFunction_[OUTER_X1](dst, time, dt, nvar, is, ie, bjs, bje, bks, bke,
                                  x0, y0, z0, dx, dy, dz);

  // Apply boundary function on inner-x2
  if (MGBoundaryFunction_[INNER_X2] != NULL)
    MGBoundaryFunction_[INNER_X2](dst, time, dt, nvar, bis, bie, js, je, bks, bke,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x2
  if (MGBoundaryFunction_[OUTER_X2] != NULL)
    MGBoundaryFunction_[OUTER_X2](dst, time, dt, nvar, bis, bie, js, je, bks, bke,
                                  x0, y0, z0, dx, dy, dz);

  bjs=js-ngh_, bje=je+ngh_;
  // Apply boundary function on inner-x3
  if (MGBoundaryFunction_[INNER_X3] != NULL)
    MGBoundaryFunction_[INNER_X3](dst, time, dt, nvar, bis, bie, bjs, bje, ks, ke,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x3
  if (MGBoundaryFunction_[OUTER_X3] != NULL)
    MGBoundaryFunction_[OUTER_X3](dst, time, dt, nvar, bis, bie, bjs, bje, ks, ke,
                                  x0, y0, z0, dx, dy, dz);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(void)
//  \brief Restrict the defect to the source
void Multigrid::Restrict(void)
{
  AthenaArray<Real> &dst=src_[current_level_-1];
  AthenaArray<Real> &src=def_[current_level_];
  int ll=nlevel_-current_level_;
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
  current_level_--;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(void)
//  \brief Prolongate the potential using tri-linear interpolation
void Multigrid::ProlongateAndCorrect(void)
{
  const AthenaArray<Real> &src=u_[current_level_];
  AthenaArray<Real> &dst=def_[current_level_+1];
  int ll=nlevel_-1-current_level_;
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
  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongate(void)
//  \brief Prolongate the potential for Full Multigrid cycle
void Multigrid::FMGProlongate(void)
{
  AthenaArray<Real> &src=u_[current_level_];
  AthenaArray<Real> &dst=u_[current_level_+1];
  int ll=nlevel_-1-current_level_;
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
  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetFMGSource(void)
//  \brief Copy the restricted Source term to the source array
void Multigrid::SetFMGSource(void)
{
  AthenaArray<Real> &src=fmgsrc_[current_level_];
  AthenaArray<Real> &dst=src_[current_level_];
  std::memcpy(dst.data(), src.data(), src.GetSizeInBytes());
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(int n, int nrm)
{
  AthenaArray<Real> &def=def_[current_level_];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;

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
\  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotalSource(int n)
//  \brief calculate the sum of the source function

Real CalculateTotalSource(int n)
{
  AthenaArray<Real> &src=src_[nlevel_-1];
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
  return s*rdx_*rdy_*rdz_;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverageSource(int type, int n, Real ave)
//  \brief subtract the average value (type: 0=source, 1=u)

void Multigrid::SubtractAverageSource(int type, int n, Real ave)
{
  AthenaArray<Real> dst;
  if(type==0) dst.InitWithShallowCopy(src_[nlevel_-1]);
  else dst.InitWithShallowCopy(u_[nlevel_-1]);
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
        dst(n,k,j,i)-=ave;
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Multigrid(MeshBlock *pmb, int invar, int nx, int ny, int nz,
//                                RegionSize isize, MGBoundaryFunc_t *MGBoundary)
//  \brief Multigrid constructor

void Multigrid::Multigrid(MeshBlock *pmb, int invar, int nx, int ny, int nz,
                          RegionSize isize, MGBoundaryFunc_t *MGBoundary)
{
  pmy_block_=pmb;
  size_=isize;
  nvar_=invar;
  ngh_=1;
  nx_=nx, ny_=ny, nz_=nz;
  rdx_=(size_.x1max-size_.x1min)/(Real)nx;
  rdy_=(size_.x2max-size_.x2min)/(Real)ny;
  rdz_=(size_.x3max-size_.x3min)/(Real)nz;

  nlevel_=0;
  int n = std::min(nx,std::min(ny, nz));
  for(int l=0; l<20; l++) {
    if((1<<l) == n) {
      nlevel_=l+1;
      break;
    }
  }

  for(int i=0; i<6; i++)
    MGBoundaryFunction_[i]=MGBoundary[i];
  if(pmb!=NULL) { // not root grid
      if(pmb->pmy_mesh->multilevel==false)
        ngh_=2;
      for(int i=0; i<6; i++) {
        if(pmb->block_bcs[i]==PERIODIC_BNDRY || pmb->block_bcs[i]==BLOCK_BNDRY)
          MGBoundaryFunction_[i]=NULL;
      }
    }
  }

  // allocate arrays
  u_ = new AthenaArray<Real>[nlevel_];
  src_ = new AthenaArray<Real>[nlevel_];
  fmgsrc_ = new AthenaArray<Real>[nlevel_];
  def_ = new AthenaArray<Real>[nlevel_];
  for(int l=nlevel_-1; l>=0; l++) {
    int ll=nlevel_-1-l;
    int ncx=(nx>>ll)+2*ngh, ncy=(ny>>ll)+2*ngh, ncz=(nz>>ll)+2*ngh;
    u_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
    fmgsrc_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar,ncz,ncy,ncx);
  }
}

//----------------------------------------------------------------------------------------
//! \fn virtual void Multigrid::~Multigrid
//  \brief Multigrid destroctor

virtual void Multigrid::~Multigrid()
{
  for(int l=0; l<nlevel_; l++) {
    u_[l].DeleteAthenaArray();
    src_[l].DeleteAthenaArray();
    fmgsrc_[l].DeleteAthenaArray();
    def_[l].DeleteAthenaArray();
  }
  delete [] u_;
  delete [] src_;
  delete [] fmgsrc_;
  delete [] def_;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X1 direction

void MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1)=dst(n,k,j,ie-i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X1 direction

void MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1)=dst(n,k,j,is+i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X2 direction

void MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=0; j<ngh; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i)=dst(n,k,je-j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X2 direction

void MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ks, k<=ke; k++) {
#pragma ivdep
      for(int j=0; j<ngh; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i)=dst(n,k,js+j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, Real dt,
//                     int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                     Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X3 direction

void MGPeriodicInnerX3(AthenaArray<Real> &dst, Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=0, k<=ngh; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i)=dst(n,ke-k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, Real dt,
//                int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X3 direction

void MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, Real dt,
                       int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=0, k<=ngh; k++) {
#pragma ivdep
      for(int j=js; j<=je; j++) {
#pragma ivdep
        for(int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i)=dst(n,ks+k,j,i);
      }
    }
  }
  return;
}

