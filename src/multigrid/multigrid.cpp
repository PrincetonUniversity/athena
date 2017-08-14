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
//! \fn Multigrid::Multigrid(Mesh *pm, MeshBlock *pmb, int invar, int nx, int ny, int nz,
//                           RegionSize isize, MGBoundaryFunc_t *MGBoundary)
//  \brief Multigrid constructor

Multigrid::Multigrid(Mesh *pm, MeshBlock *pmb, int invar, int nx, int ny, int nz,
                     int nghost, RegionSize isize, MGBoundaryFunc_t *MGBoundary)
{
  pmy_mesh_=pm;
  pmy_block_=pmb;
  ngh_=nghost;
  size_=isize;
  nvar_=invar;
  nx_=nx, ny_=ny, nz_=nz;
  rdx_=(size_.x1max-size_.x1min)/(Real)nx;
  rdy_=(size_.x2max-size_.x2min)/(Real)ny;
  rdz_=(size_.x3max-size_.x3min)/(Real)nz;

  nlevel_=0;
  for(int l=0; l<20; l++) {
    if(nx_%(1<<l)==0 && ny_%(1<<l)==0 && nz_%(1<<l)==0) {
      nlevel_=l+1;
    }
  }
  for(int i=0; i<6; i++)
    MGBoundaryFunction_[i]=MGBoundary[i];
  if(pmb!=NULL) { // not root grid
    for(int i=0; i<6; i++) {
      if(pmb->pbval->block_bcs[i]==PERIODIC_BNDRY
      || pmb->pbval->block_bcs[i]==BLOCK_BNDRY)
        MGBoundaryFunction_[i]=NULL;
    }
  }

  // allocate arrays
  u_ = new AthenaArray<Real>[nlevel_];
  src_ = new AthenaArray<Real>[nlevel_];
  def_ = new AthenaArray<Real>[nlevel_];
  for(int l=0; l<nlevel_; l++) {
    int ll=nlevel_-1-l;
    int ncx=(nx>>ll)+2*ngh_, ncy=(ny>>ll)+2*ngh_, ncz=(nz>>ll)+2*ngh_;
    u_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    src_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
  }
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid::~Multigrid
//  \brief Multigrid destroctor

Multigrid::~Multigrid()
{
  for(int l=0; l<nlevel_; l++) {
    u_[l].DeleteAthenaArray();
    src_[l].DeleteAthenaArray();
    def_[l].DeleteAthenaArray();
  }
  delete [] u_;
  delete [] src_;
  delete [] def_;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
{
  AthenaArray<Real> &dst=u_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
  for(int n=0; n<nvar_; n++) {
    int nsrc=ns+n;
    for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
      for(int j=ngh, mj=js; mj<=je; j++, mj++) {
        for(int i=ngh, mi=is; mi<=ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                 Real fac)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac)
{
  AthenaArray<Real> &dst=src_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;
  if(fac==1.0) {
    for(int n=0; n<nvar_; n++) {
      int nsrc=ns+n;
      for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
        for(int j=ngh, mj=js; mj<=je; j++, mj++) {
          for(int i=ngh, mi=is; mi<=ie; i++, mi++)
            dst(n,mk,mj,mi)=src(nsrc,k,j,i);
        }
      }
    }
  }
  else {
    for(int n=0; n<nvar_; n++) {
      int nsrc=ns+n;
      for(int k=ngh, mk=ks; mk<=ke; k++, mk++) {
        for(int j=ngh, mj=js; mj<=je; j++, mj++) {
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
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  current_level_=nlevel_-1;
  for(; current_level_>0; current_level_--) {
    int ll=nlevel_-current_level_;
    ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
    AthenaArray<Real> &csrc=src_[current_level_-1];
    const AthenaArray<Real> &fsrc=src_[current_level_];
    for(int n=0; n<nvar_; n++) {
      for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
        for(int j=js, fj=js; j<=je; j++, fj+=2) {
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
  const AthenaArray<Real> &src=u_[nlevel_-1];
  int sngh=std::min(ngh_,ngh);
  int ie=nx_+ngh_+sngh-1, je=ny_+ngh_+sngh-1, ke=nz_+ngh_+sngh-1;
  for(int n=0; n<nvar_; n++) {
    int ndst=ns+n;
    for(int k=ngh-sngh, mk=ngh_-sngh; mk<=ke; k++, mk++) {
      for(int j=ngh-sngh, mj=ngh_-sngh; mj<=je; j++, mj++) {
        for(int i=ngh-sngh, mi=ngh_-sngh; mi<=ie; i++, mi++)
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
  int is=ngh_, ie=ncx+ngh_-1, js=ngh_, je=ncy+ngh_-1, ks=ngh_, ke=ncz+ngh_-1;
  int bis=is-ngh_, bie=ie+ngh_, bjs=js, bje=je, bks=ks, bke=ke;
  Real dx=rdx_*(Real)(1<<ll), dy=rdy_*(Real)(1<<ll), dz=rdz_*(Real)(1<<ll);
  Real x0=size_.x1min-((Real)ngh_+0.5)*dx;
  Real y0=size_.x2min-((Real)ngh_+0.5)*dy;
  Real z0=size_.x3min-((Real)ngh_+0.5)*dz;
  Real time=pmy_mesh_->time;
  if(MGBoundaryFunction_[INNER_X2]==NULL) bjs=js-ngh_;
  if(MGBoundaryFunction_[OUTER_X2]==NULL) bje=je+ngh_;
  if(MGBoundaryFunction_[INNER_X3]==NULL) bks=ks-ngh_;
  if(MGBoundaryFunction_[OUTER_X3]==NULL) bke=ke+ngh_;

  // Apply boundary function on inner-x1
  if (MGBoundaryFunction_[INNER_X1] != NULL)
    MGBoundaryFunction_[INNER_X1](dst, time, nvar_, is, ie, bjs, bje, bks, bke, ngh_,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x1
  if (MGBoundaryFunction_[OUTER_X1] != NULL)
    MGBoundaryFunction_[OUTER_X1](dst, time, nvar_, is, ie, bjs, bje, bks, bke, ngh_,
                                  x0, y0, z0, dx, dy, dz);

  // Apply boundary function on inner-x2
  if (MGBoundaryFunction_[INNER_X2] != NULL)
    MGBoundaryFunction_[INNER_X2](dst, time, nvar_, bis, bie, js, je, bks, bke, ngh_,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x2
  if (MGBoundaryFunction_[OUTER_X2] != NULL)
    MGBoundaryFunction_[OUTER_X2](dst, time, nvar_, bis, bie, js, je, bks, bke, ngh_,
                                  x0, y0, z0, dx, dy, dz);

  bjs=js-ngh_, bje=je+ngh_;
  // Apply boundary function on inner-x3
  if (MGBoundaryFunction_[INNER_X3] != NULL)
    MGBoundaryFunction_[INNER_X3](dst, time, nvar_, bis, bie, bjs, bje, ks, ke, ngh_,
                                  x0, y0, z0, dx, dy, dz);
  // Apply boundary function on outer-x3
  if (MGBoundaryFunction_[OUTER_X3] != NULL)
    MGBoundaryFunction_[OUTER_X3](dst, time, nvar_, bis, bie, bjs, bje, ks, ke, ngh_,
                                  x0, y0, z0, dx, dy, dz);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(void)
//  \brief Restrict the defect to the source
void Multigrid::Restrict(void)
{
  AthenaArray<Real> &dst=src_[current_level_-1];
  const AthenaArray<Real> &src=def_[current_level_];
  int ll=nlevel_-current_level_;
  int is, ie, js, je, ks, ke;

  CalculateDefect();
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
  for(int n=0; n<nvar_; n++) {
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
        for(int i=is, fi=is; i<=ie; i++, fi+=2)
          dst(n, k, j, i)=0.125*(src(n, fk,   fj,   fi)+src(n, fk,   fj,   fi+1)
                                +src(n, fk,   fj+1, fi)+src(n, fk,   fj+1, fi+1)
                                +src(n, fk+1, fj,   fi)+src(n, fk+1, fj,   fi+1)
                                +src(n, fk+1, fj+1, fi)+src(n, fk+1, fj+1, fi+1));
      }
    }
  }
  current_level_--;
  ZeroClearData();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(void)
//  \brief Prolongate the potential using tri-linear interpolation
void Multigrid::ProlongateAndCorrect(void)
{
  const AthenaArray<Real> &src=u_[current_level_];
  AthenaArray<Real> &dst=u_[current_level_+1];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
  for(int n=0; n<nvar_; n++) {
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
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
  const AthenaArray<Real> &src=u_[current_level_];
  AthenaArray<Real> &dst=u_[current_level_+1];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
  for(int n=0; n<nvar_; n++) {
    for(int k=ks, fk=ks; k<=ke; k++, fk+=2) {
      for(int j=js, fj=js; j<=je; j++, fj+=2) {
        for(int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj,  fi  )=(
          + 125.0*src(n,k-1,j-1,i-1)+  750.0*src(n,k-1,j-1,i  )-  75.0*src(n,k-1,j-1,i+1)
          + 750.0*src(n,k-1,j,  i-1)+ 4500.0*src(n,k-1,j,  i  )- 450.0*src(n,k-1,j,  i+1)
          -  75.0*src(n,k-1,j+1,i-1)-  450.0*src(n,k-1,j+1,i  )+  45.0*src(n,k-1,j+1,i+1)
          + 750.0*src(n,k,  j-1,i-1)+ 4500.0*src(n,k,  j-1,i  )- 450.0*src(n,k,  j-1,i+1)
          +4500.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )-2700.0*src(n,k,  j,  i+1)
          - 450.0*src(n,k,  j+1,i-1)- 2700.0*src(n,k,  j+1,i  )+ 270.0*src(n,k,  j+1,i+1)
          -  75.0*src(n,k+1,j-1,i-1)-  450.0*src(n,k+1,j-1,i  )+  45.0*src(n,k+1,j-1,i+1)
          - 450.0*src(n,k+1,j,  i-1)- 2700.0*src(n,k+1,j,  i  )+ 270.0*src(n,k+1,j,  i+1)
          +  45.0*src(n,k+1,j+1,i-1)+  270.0*src(n,k+1,j+1,i  )-  27.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk,  fj,  fi+1)=(
          -  75.0*src(n,k-1,j-1,i-1)+  750.0*src(n,k-1,j-1,i  )+ 125.0*src(n,k-1,j-1,i+1)
          - 450.0*src(n,k-1,j,  i-1)+ 4500.0*src(n,k-1,j,  i  )+ 750.0*src(n,k-1,j,  i+1)
          +  45.0*src(n,k-1,j+1,i-1)-  450.0*src(n,k-1,j+1,i  )-  75.0*src(n,k-1,j+1,i+1)
          - 450.0*src(n,k,  j-1,i-1)+ 4500.0*src(n,k,  j-1,i  )+ 750.0*src(n,k,  j-1,i+1)
          -2700.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )+4500.0*src(n,k,  j,  i+1)
          + 270.0*src(n,k,  j+1,i-1)- 2700.0*src(n,k,  j+1,i  )- 450.0*src(n,k,  j+1,i+1)
          +  45.0*src(n,k+1,j-1,i-1)-  450.0*src(n,k+1,j-1,i  )-  75.0*src(n,k+1,j-1,i+1)
          + 270.0*src(n,k+1,j,  i-1)- 2700.0*src(n,k+1,j,  i  )- 450.0*src(n,k+1,j,  i+1)
          -  27.0*src(n,k+1,j+1,i-1)+  270.0*src(n,k+1,j+1,i  )+  45.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk  ,fj+1,fi  )=(
          -  75.0*src(n,k-1,j-1,i-1)-  450.0*src(n,k-1,j-1,i  )+  45.0*src(n,k-1,j-1,i+1)
          + 750.0*src(n,k-1,j,  i-1)+ 4500.0*src(n,k-1,j,  i  )- 450.0*src(n,k-1,j,  i+1)
          + 125.0*src(n,k-1,j+1,i-1)+  750.0*src(n,k-1,j+1,i  )-  75.0*src(n,k-1,j+1,i+1)
          - 450.0*src(n,k,  j-1,i-1)- 2700.0*src(n,k,  j-1,i  )+ 270.0*src(n,k,  j-1,i+1)
          +4500.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )-2700.0*src(n,k,  j,  i+1)
          + 750.0*src(n,k,  j+1,i-1)+ 4500.0*src(n,k,  j+1,i  )- 450.0*src(n,k,  j+1,i+1)
          +  45.0*src(n,k+1,j-1,i-1)+  270.0*src(n,k+1,j-1,i  )-  27.0*src(n,k+1,j-1,i+1)
          - 450.0*src(n,k+1,j,  i-1)- 2700.0*src(n,k+1,j,  i  )+ 270.0*src(n,k+1,j,  i+1)
          -  75.0*src(n,k+1,j+1,i-1)-  450.0*src(n,k+1,j+1,i  )+  45.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk,  fj+1,fi+1)=(
          +  45.0*src(n,k-1,j-1,i-1)-  450.0*src(n,k-1,j-1,i  )-  75.0*src(n,k-1,j-1,i+1)
          - 450.0*src(n,k-1,j,  i-1)+ 4500.0*src(n,k-1,j,  i  )+ 750.0*src(n,k-1,j,  i+1)
          -  75.0*src(n,k-1,j+1,i-1)+  750.0*src(n,k-1,j+1,i  )+ 125.0*src(n,k-1,j+1,i+1)
          + 270.0*src(n,k,  j-1,i-1)- 2700.0*src(n,k,  j-1,i  )- 450.0*src(n,k,  j-1,i+1)
          -2700.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )+4500.0*src(n,k,  j,  i+1)
          - 450.0*src(n,k,  j+1,i-1)+ 4500.0*src(n,k,  j+1,i  )+ 750.0*src(n,k,  j+1,i+1)
          -  27.0*src(n,k+1,j-1,i-1)+  270.0*src(n,k+1,j-1,i  )+  45.0*src(n,k+1,j-1,i+1)
          + 270.0*src(n,k+1,j,  i-1)- 2700.0*src(n,k+1,j,  i  )- 450.0*src(n,k+1,j,  i+1)
          +  45.0*src(n,k+1,j+1,i-1)-  450.0*src(n,k+1,j+1,i  )-  75.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk+1,fj,  fi  )=(
          -  75.0*src(n,k-1,j-1,i-1)-  450.0*src(n,k-1,j-1,i  )+  45.0*src(n,k-1,j-1,i+1)
          - 450.0*src(n,k-1,j,  i-1)- 2700.0*src(n,k-1,j,  i  )+ 270.0*src(n,k-1,j,  i+1)
          +  45.0*src(n,k-1,j+1,i-1)+  270.0*src(n,k-1,j+1,i  )-  27.0*src(n,k-1,j+1,i+1)
          + 750.0*src(n,k,  j-1,i-1)+ 4500.0*src(n,k,  j-1,i  )- 450.0*src(n,k,  j-1,i+1)
          +4500.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )-2700.0*src(n,k,  j,  i+1)
          - 450.0*src(n,k,  j+1,i-1)- 2700.0*src(n,k,  j+1,i  )+ 270.0*src(n,k,  j+1,i+1)
          + 125.0*src(n,k+1,j-1,i-1)+  750.0*src(n,k+1,j-1,i  )-  75.0*src(n,k+1,j-1,i+1)
          + 750.0*src(n,k+1,j,  i-1)+ 4500.0*src(n,k+1,j,  i  )- 450.0*src(n,k+1,j,  i+1)
          -  75.0*src(n,k+1,j+1,i-1)-  450.0*src(n,k+1,j+1,i  )+  45.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk+1,fj,  fi+1)=(
          +  45.0*src(n,k-1,j-1,i-1)-  450.0*src(n,k-1,j-1,i  )-  75.0*src(n,k-1,j-1,i+1)
          + 270.0*src(n,k-1,j,  i-1)- 2700.0*src(n,k-1,j,  i  )- 450.0*src(n,k-1,j,  i+1)
          -  27.0*src(n,k-1,j+1,i-1)+  270.0*src(n,k-1,j+1,i  )+  45.0*src(n,k-1,j+1,i+1)
          - 450.0*src(n,k,  j-1,i-1)+ 4500.0*src(n,k,  j-1,i  )+ 750.0*src(n,k,  j-1,i+1)
          -2700.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )+4500.0*src(n,k,  j,  i+1)
          + 270.0*src(n,k,  j+1,i-1)- 2700.0*src(n,k,  j+1,i  )- 450.0*src(n,k,  j+1,i+1)
          -  75.0*src(n,k+1,j-1,i-1)+  750.0*src(n,k+1,j-1,i  )+ 125.0*src(n,k+1,j-1,i+1)
          - 450.0*src(n,k+1,j,  i-1)+ 4500.0*src(n,k+1,j,  i  )+ 750.0*src(n,k+1,j,  i+1)
          +  45.0*src(n,k+1,j+1,i-1)-  450.0*src(n,k+1,j+1,i  )-  75.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk+1,fj+1,fi  )=(
          +  45.0*src(n,k-1,j-1,i-1)+  270.0*src(n,k-1,j-1,i  )-  27.0*src(n,k-1,j-1,i+1)
          - 450.0*src(n,k-1,j,  i-1)- 2700.0*src(n,k-1,j,  i  )+ 270.0*src(n,k-1,j,  i+1)
          -  75.0*src(n,k-1,j+1,i-1)-  450.0*src(n,k-1,j+1,i  )+  45.0*src(n,k-1,j+1,i+1)
          - 450.0*src(n,k,  j-1,i-1)- 2700.0*src(n,k,  j-1,i  )+ 270.0*src(n,k,  j-1,i+1)
          +4500.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )-2700.0*src(n,k,  j,  i+1)
          + 750.0*src(n,k,  j+1,i-1)+ 4500.0*src(n,k,  j+1,i  )- 450.0*src(n,k,  j+1,i+1)
          -  75.0*src(n,k+1,j-1,i-1)-  450.0*src(n,k+1,j-1,i  )+  45.0*src(n,k+1,j-1,i+1)
          + 750.0*src(n,k+1,j,  i-1)+ 4500.0*src(n,k+1,j,  i  )- 450.0*src(n,k+1,j,  i+1)
          + 125.0*src(n,k+1,j+1,i-1)+  750.0*src(n,k+1,j+1,i  )-  75.0*src(n,k+1,j+1,i+1)
            )/32768.0;
          dst(n,fk+1,fj+1,fi+1)=(
          -  27.0*src(n,k-1,j-1,i-1)+  270.0*src(n,k-1,j-1,i  )+  45.0*src(n,k-1,j-1,i+1)
          + 270.0*src(n,k-1,j,  i-1)- 2700.0*src(n,k-1,j,  i  )- 450.0*src(n,k-1,j,  i+1)
          +  45.0*src(n,k-1,j+1,i-1)-  450.0*src(n,k-1,j+1,i  )-  75.0*src(n,k-1,j+1,i+1)
          + 270.0*src(n,k,  j-1,i-1)- 2700.0*src(n,k,  j-1,i  )- 450.0*src(n,k,  j-1,i+1)
          -2700.0*src(n,k,  j,  i-1)+27000.0*src(n,k,  j,  i  )+4500.0*src(n,k,  j,  i+1)
          - 450.0*src(n,k,  j+1,i-1)+ 4500.0*src(n,k,  j+1,i  )+ 750.0*src(n,k,  j+1,i+1)
          +  45.0*src(n,k+1,j-1,i-1)-  450.0*src(n,k+1,j-1,i  )-  75.0*src(n,k+1,j-1,i+1)
          - 450.0*src(n,k+1,j,  i-1)+ 4500.0*src(n,k+1,j,  i  )+ 750.0*src(n,k+1,j,  i+1)
          -  75.0*src(n,k+1,j+1,i-1)+  750.0*src(n,k+1,j+1,i  )+ 125.0*src(n,k+1,j+1,i+1)
            )/32768.0;
        }
      }
    }
  }
  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck)
//  \brief Load the data from the root grid
void Multigrid::SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck)
{
  current_level_=0;
  AthenaArray<Real> &dst=u_[current_level_];
  for(int n=0; n<nvar_; n++) {
    dst(n,ngh_,ngh_,ngh_)=src(n,ck+ngh_,cj+ngh_,ci+ngh_);
  }
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
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++)
          norm=std::max(norm,std::fabs(def(n,k,j,i)));
      }
    }
  }
  else if (nrm==1) {
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++)
          norm+=std::fabs(def(n,k,j,i));
      }
    }
  }
  else { // nrm>1 -> nrm=2
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++)
          norm+=SQR(def(n,k,j,i));
      }
    }
  }
  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotal(int type, int n)
//  \brief calculate the sum of the array (type: 0=src, 1=u)

Real Multigrid::CalculateTotal(int type, int n)
{
  AthenaArray<Real> src;
  int ll=nlevel_-1-current_level_;
  if(type==0) src.InitWithShallowCopy(src_[current_level_]);
  else src.InitWithShallowCopy(u_[current_level_]);
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
  Real dx=rdx_*(Real)(1<<ll), dy=rdy_*(Real)(1<<ll), dz=rdz_*(Real)(1<<ll);
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++)
        s+=src(n,k,j,i);
    }
  }
  return s*dx*dy*dz;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverage(int type, int n, Real ave)
//  \brief subtract the average value (type: 0=source, 1=u)

void Multigrid::SubtractAverage(int type, int n, Real ave)
{
  AthenaArray<Real> dst;
  if(type==0) dst.InitWithShallowCopy(src_[nlevel_-1]);
  else dst.InitWithShallowCopy(u_[nlevel_-1]);
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+nx_-1, je=js+ny_-1, ke=ks+nz_-1;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++)
        dst(n,k,j,i)-=ave;
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X1 direction

void MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
  for(int n=0; n<nvar; n++) {
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1)=dst(n,k,j,ie-i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X1 direction

void MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
  for(int n=0; n<nvar; n++) {
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1)=dst(n,k,j,is+i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X2 direction

void MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
  for(int n=0; n<nvar; n++) {
    for(int k=ks; k<=ke; k++) {
      for(int j=0; j<ngh; j++) {
        for(int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i)=dst(n,k,je-j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X2 direction

void MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
  for(int n=0; n<nvar; n++) {
    for(int k=ks; k<=ke; k++) {
      for(int j=0; j<ngh; j++) {
        for(int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i)=dst(n,k,js+j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X3 direction

void MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
  for(int n=0; n<nvar; n++) {
    for(int k=0; k<ngh; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i)=dst(n,ke-k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X3 direction

void MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
{
  for(int n=0; n<nvar; n++) {
    for(int k=0; k<ngh; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i)=dst(n,ks+k,j,i);
      }
    }
  }
  return;
}

