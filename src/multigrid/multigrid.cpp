//========================================================================================
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>    // memset
#include <iostream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

//----------------------------------------------------------------------------------------
//! \fn Multigrid::Multigrid(MultigridDriver *pmd, LogicalLocation iloc, int igid,
//                           int ilid, int invar, int nghost, RegionSize isize,
//                           MGBoundaryFunc *MGBoundary, BoundaryFlag *input_bcs,
//                           bool root);
//  \brief Multigrid constructor

Multigrid::Multigrid(
    MultigridDriver *pmd, LogicalLocation iloc, int igid, int ilid,
    int invar, int nghost, RegionSize isize, MGBoundaryFunc *MGBoundary,
    BoundaryFlag *input_bcs, bool root) {
  pmy_driver_=pmd;
  loc_=iloc;
  gid_=igid;
  lid_=ilid;
  ngh_=nghost;
  size_=isize;
  nvar_=invar;
  root_flag_=root;
  rdx_=(size_.x1max-size_.x1min)/static_cast<Real>(size_.nx1);
  rdy_=(size_.x2max-size_.x2min)/static_cast<Real>(size_.nx2);
  rdz_=(size_.x3max-size_.x3min)/static_cast<Real>(size_.nx3);
  prev=nullptr;
  next=nullptr;

  nlevel_=0;
  if (root_flag_ == true) {
    int nbx, nby, nbz;
    for (int l=0; l<20; l++) {
      if (size_.nx1%(1<<l)==0 && size_.nx2%(1<<l)==0 && size_.nx3%(1<<l)==0) {
        nbx=size_.nx1/(1<<l), nby=size_.nx2/(1<<l), nbz=size_.nx3/(1<<l);
        nlevel_=l+1;
      }
    }
    int nmaxr=std::max(nbx, std::max(nby, nbz));
    // int nminr=std::min(nbx, std::min(nby, nbz)); // unused variable
    if (nmaxr!=1 && Globals::my_rank==0) {
      std::cout
          << "### Warning in Multigrid::Multigrid" << std::endl
          << "The root grid can not be reduced to a single cell." << std::endl
          << "Multigrid should still work, but this is not the"
          << " most efficient configuration"
          << " as the coarsest level is not solved exactly but iteratively." << std::endl;
    }
    if (nbx*nby*nbz>100 && Globals::my_rank==0) {
      std::cout << "### Warning in Multigrid::Multigrid" << std::endl
                << "The degrees of freedom on the coarsest level is very large: "
                << nbx << " x " << nby << " x " << nbz << " = " << nbx*nby*nbz<< std::endl
                << "Multigrid should still work, but this is not efficient configuration "
                << "as the coarsest level solver costs considerably." << std::endl
                << "We recommend to reconsider grid configuration." << std::endl;
    }
  } else {
    for (int l=0; l<20; l++) {
      if ((1<<l) == size_.nx1) {
        nlevel_=l+1;
        break;
      }
    }
    if (nlevel_==0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Multigrid::Multigrid" << std::endl
          << "The MeshBlock size must be power of two." << std::endl;
      ATHENA_ERROR(msg);
      return;
    }
    // *** temporary ***
    if (std::fabs(rdx_-rdy_)>1.0e-5 || std::fabs(rdx_-rdz_)>1.0e-5) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Multigrid::Multigrid" << std::endl
          << "The cell size must be cubic." << std::endl;
      ATHENA_ERROR(msg);
      return;
    }
  }

  pmgbval = new MGBoundaryValues(this, input_bcs, MGBoundary);

  // allocate arrays
  u_ = new AthenaArray<Real>[nlevel_];
  src_ = new AthenaArray<Real>[nlevel_];
  def_ = new AthenaArray<Real>[nlevel_];
  for (int l=0; l<nlevel_; l++) {
    int ll=nlevel_-1-l;
    int ncx=(size_.nx1>>ll)+2*ngh_, ncy=(size_.nx2>>ll)+2*ngh_,
        ncz=(size_.nx3>>ll)+2*ngh_;
    u_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    src_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
  }
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid::~Multigrid
//  \brief Multigrid destroctor

Multigrid::~Multigrid() {
  if (prev!=nullptr) prev->next=next;
  if (next!=nullptr) next->prev=prev;

  delete [] u_;
  delete [] src_;
  delete [] def_;
  delete pmgbval;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the inital guess in the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh) {
  AthenaArray<Real> &dst=u_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  for (int n=0; n<nvar_; n++) {
    int nsrc=ns+n;
    for (int k=ngh, mk=ks; mk<=ke; k++, mk++) {
      for (int j=ngh, mj=js; mj<=je; j++, mj++) {
        for (int i=ngh, mi=is; mi<=ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                 Real fac)
//  \brief Fill the source in the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac) {
  AthenaArray<Real> &dst=src_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  if (fac==1.0) {
    for (int n=0; n<nvar_; n++) {
      int nsrc=ns+n;
      for (int k=ngh, mk=ks; mk<=ke; k++, mk++) {
        for (int j=ngh, mj=js; mj<=je; j++, mj++) {
          for (int i=ngh, mi=is; mi<=ie; i++, mi++)
            dst(n,mk,mj,mi)=src(nsrc,k,j,i);
        }
      }
    }
  } else {
    for (int n=0; n<nvar_; n++) {
      int nsrc=ns+n;
      for (int k=ngh, mk=ks; mk<=ke; k++, mk++) {
        for (int j=ngh, mj=js; mj<=je; j++, mj++) {
          for (int i=ngh, mi=is; mi<=ie; i++, mi++)
            dst(n,mk,mj,mi)=src(nsrc,k,j,i)*fac;
        }
      }
    }
  }
  current_level_=nlevel_-1;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictFMGSource()
//  \brief restrict the source through all the multigrid levels

void Multigrid::RestrictFMGSource() {
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  current_level_=nlevel_-1;
  for (; current_level_>0; current_level_--) {
    int ll=nlevel_-current_level_;
    ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
    AthenaArray<Real> &csrc=src_[current_level_-1];
    const AthenaArray<Real> &fsrc=src_[current_level_];
    for (int n=0; n<nvar_; n++) {
      for (int k=ks, fk=ks; k<=ke; k++, fk+=2) {
        for (int j=js, fj=js; j<=je; j++, fj+=2) {
          for (int i=is, fi=is; i<=ie; i++, fi+=2)
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
void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh) {
  const AthenaArray<Real> &src=u_[nlevel_-1];
  int sngh=std::min(ngh_,ngh);
  int ie=size_.nx1+ngh_+sngh-1, je=size_.nx2+ngh_+sngh-1, ke=size_.nx3+ngh_+sngh-1;
  for (int n=0; n<nvar_; n++) {
    int ndst=ns+n;
    for (int k=ngh-sngh, mk=ngh_-sngh; mk<=ke; k++, mk++) {
      for (int j=ngh-sngh, mj=ngh_-sngh; mj<=je; j++, mj++) {
        for (int i=ngh-sngh, mi=ngh_-sngh; mi<=ie; i++, mi++)
          dst(ndst,k,j,i)=src(n,mk,mj,mi);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ZeroClearData()
//  \brief Clear the data array with zero
void Multigrid::ZeroClearData() {
  AthenaArray<Real> &u=u_[current_level_];
  std::memset(u.data(), 0, u.GetSizeInBytes());
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict()
//  \brief Restrict the defect to the source
void Multigrid::Restrict() {
  AthenaArray<Real> &dst=src_[current_level_-1];
  const AthenaArray<Real> &src=def_[current_level_];
  int ll=nlevel_-current_level_;
  int is, ie, js, je, ks, ke;

  CalculateDefect();
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  for (int n=0; n<nvar_; n++) {
    for (int k=ks, fk=ks; k<=ke; k++, fk+=2) {
      for (int j=js, fj=js; j<=je; j++, fj+=2) {
        for (int i=is, fi=is; i<=ie; i++, fi+=2)
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
//! \fn void Multigrid::ProlongateAndCorrect()
//  \brief Prolongate the potential using tri-linear interpolation
void Multigrid::ProlongateAndCorrect() {
  const AthenaArray<Real> &src=u_[current_level_];
  AthenaArray<Real> &dst=u_[current_level_+1];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  for (int n=0; n<nvar_; n++) {
    for (int k=ks, fk=ks; k<=ke; k++, fk+=2) {
      for (int j=js, fj=js; j<=je; j++, fj+=2) {
        for (int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj  ,fi  ) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i-1)
                        +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                        +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk  ,fj  ,fi+1) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j-1,i+1)
                        +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k-1,j,i))
                        +3.0*(src(n,k-1,j-1,i)+src(n,k-1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi  ) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i-1)
                        +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                        +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi  ) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i-1)
                        +9.0*(src(n,k,j,i-1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                        +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i-1)+src(n,k,j-1,i-1)));
          dst(n,fk+1,fj+1,fi  ) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i-1)
                        +9.0*(src(n,k,j,i-1)+src(n,k,j+1,i)+src(n,k+1,j,i))
                        +3.0*(src(n,k+1,j+1,i)+src(n,k+1,j,i-1)+src(n,k,j+1,i-1)));
          dst(n,fk+1,fj  ,fi+1) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j-1,i+1)
                        +9.0*(src(n,k,j,i+1)+src(n,k,j-1,i)+src(n,k+1,j,i))
                        +3.0*(src(n,k+1,j-1,i)+src(n,k+1,j,i+1)+src(n,k,j-1,i+1)));
          dst(n,fk  ,fj+1,fi+1) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k-1,j+1,i+1)
                        +9.0*(src(n,k,j,i+1)+src(n,k,j+1,i)+src(n,k-1,j,i))
                        +3.0*(src(n,k-1,j+1,i)+src(n,k-1,j,i+1)+src(n,k,j+1,i+1)));
          dst(n,fk+1,fj+1,fi+1) +=
              0.015625*(27.0*src(n,k,j,i) + src(n,k+1,j+1,i+1)
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
//! \fn void Multigrid::FMGProlongate()
//  \brief Prolongate the potential for Full Multigrid cycle
void Multigrid::FMGProlongate() {
  const AthenaArray<Real> &src=u_[current_level_];
  AthenaArray<Real> &dst=u_[current_level_+1];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  for (int n=0; n<nvar_; n++) {
    for (int k=ks, fk=ks; k<=ke; k++, fk+=2) {
      for (int j=js, fj=js; j<=je; j++, fj+=2) {
        for (int i=is, fi=is; i<=ie; i++, fi+=2) {
          dst(n,fk  ,fj,  fi  )=(
              + 125.*src(n,k-1,j-1,i-1)+  750.*src(n,k-1,j-1,i  )-  75.*src(n,k-1,j-1,i+1)
              + 750.*src(n,k-1,j,  i-1)+ 4500.*src(n,k-1,j,  i  )- 450.*src(n,k-1,j,  i+1)
              -  75.*src(n,k-1,j+1,i-1)-  450.*src(n,k-1,j+1,i  )+  45.*src(n,k-1,j+1,i+1)
              + 750.*src(n,k,  j-1,i-1)+ 4500.*src(n,k,  j-1,i  )- 450.*src(n,k,  j-1,i+1)
              +4500.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )-2700.*src(n,k,  j,  i+1)
              - 450.*src(n,k,  j+1,i-1)- 2700.*src(n,k,  j+1,i  )+ 270.*src(n,k,  j+1,i+1)
              -  75.*src(n,k+1,j-1,i-1)-  450.*src(n,k+1,j-1,i  )+  45.*src(n,k+1,j-1,i+1)
              - 450.*src(n,k+1,j,  i-1)- 2700.*src(n,k+1,j,  i  )+ 270.*src(n,k+1,j,  i+1)
              +  45.*src(n,k+1,j+1,i-1)+  270.*src(n,k+1,j+1,i  )-  27.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk,  fj,  fi+1)=(
              -  75.*src(n,k-1,j-1,i-1)+  750.*src(n,k-1,j-1,i  )+ 125.*src(n,k-1,j-1,i+1)
              - 450.*src(n,k-1,j,  i-1)+ 4500.*src(n,k-1,j,  i  )+ 750.*src(n,k-1,j,  i+1)
              +  45.*src(n,k-1,j+1,i-1)-  450.*src(n,k-1,j+1,i  )-  75.*src(n,k-1,j+1,i+1)
              - 450.*src(n,k,  j-1,i-1)+ 4500.*src(n,k,  j-1,i  )+ 750.*src(n,k,  j-1,i+1)
              -2700.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )+4500.*src(n,k,  j,  i+1)
              + 270.*src(n,k,  j+1,i-1)- 2700.*src(n,k,  j+1,i  )- 450.*src(n,k,  j+1,i+1)
              +  45.*src(n,k+1,j-1,i-1)-  450.*src(n,k+1,j-1,i  )-  75.*src(n,k+1,j-1,i+1)
              + 270.*src(n,k+1,j,  i-1)- 2700.*src(n,k+1,j,  i  )- 450.*src(n,k+1,j,  i+1)
              -  27.*src(n,k+1,j+1,i-1)+  270.*src(n,k+1,j+1,i  )+  45.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk  ,fj+1,fi  )=(
              -  75.*src(n,k-1,j-1,i-1)-  450.*src(n,k-1,j-1,i  )+  45.*src(n,k-1,j-1,i+1)
              + 750.*src(n,k-1,j,  i-1)+ 4500.*src(n,k-1,j,  i  )- 450.*src(n,k-1,j,  i+1)
              + 125.*src(n,k-1,j+1,i-1)+  750.*src(n,k-1,j+1,i  )-  75.*src(n,k-1,j+1,i+1)
              - 450.*src(n,k,  j-1,i-1)- 2700.*src(n,k,  j-1,i  )+ 270.*src(n,k,  j-1,i+1)
              +4500.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )-2700.*src(n,k,  j,  i+1)
              + 750.*src(n,k,  j+1,i-1)+ 4500.*src(n,k,  j+1,i  )- 450.*src(n,k,  j+1,i+1)
              +  45.*src(n,k+1,j-1,i-1)+  270.*src(n,k+1,j-1,i  )-  27.*src(n,k+1,j-1,i+1)
              - 450.*src(n,k+1,j,  i-1)- 2700.*src(n,k+1,j,  i  )+ 270.*src(n,k+1,j,  i+1)
              -  75.*src(n,k+1,j+1,i-1)-  450.*src(n,k+1,j+1,i  )+  45.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk,  fj+1,fi+1)=(
              +  45.*src(n,k-1,j-1,i-1)-  450.*src(n,k-1,j-1,i  )-  75.*src(n,k-1,j-1,i+1)
              - 450.*src(n,k-1,j,  i-1)+ 4500.*src(n,k-1,j,  i  )+ 750.*src(n,k-1,j,  i+1)
              -  75.*src(n,k-1,j+1,i-1)+  750.*src(n,k-1,j+1,i  )+ 125.*src(n,k-1,j+1,i+1)
              + 270.*src(n,k,  j-1,i-1)- 2700.*src(n,k,  j-1,i  )- 450.*src(n,k,  j-1,i+1)
              -2700.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )+4500.*src(n,k,  j,  i+1)
              - 450.*src(n,k,  j+1,i-1)+ 4500.*src(n,k,  j+1,i  )+ 750.*src(n,k,  j+1,i+1)
              -  27.*src(n,k+1,j-1,i-1)+  270.*src(n,k+1,j-1,i  )+  45.*src(n,k+1,j-1,i+1)
              + 270.*src(n,k+1,j,  i-1)- 2700.*src(n,k+1,j,  i  )- 450.*src(n,k+1,j,  i+1)
              +  45.*src(n,k+1,j+1,i-1)-  450.*src(n,k+1,j+1,i  )-  75.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk+1,fj,  fi  )=(
              -  75.*src(n,k-1,j-1,i-1)-  450.*src(n,k-1,j-1,i  )+  45.*src(n,k-1,j-1,i+1)
              - 450.*src(n,k-1,j,  i-1)- 2700.*src(n,k-1,j,  i  )+ 270.*src(n,k-1,j,  i+1)
              +  45.*src(n,k-1,j+1,i-1)+  270.*src(n,k-1,j+1,i  )-  27.*src(n,k-1,j+1,i+1)
              + 750.*src(n,k,  j-1,i-1)+ 4500.*src(n,k,  j-1,i  )- 450.*src(n,k,  j-1,i+1)
              +4500.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )-2700.*src(n,k,  j,  i+1)
              - 450.*src(n,k,  j+1,i-1)- 2700.*src(n,k,  j+1,i  )+ 270.*src(n,k,  j+1,i+1)
              + 125.*src(n,k+1,j-1,i-1)+  750.*src(n,k+1,j-1,i  )-  75.*src(n,k+1,j-1,i+1)
              + 750.*src(n,k+1,j,  i-1)+ 4500.*src(n,k+1,j,  i  )- 450.*src(n,k+1,j,  i+1)
              -  75.*src(n,k+1,j+1,i-1)-  450.*src(n,k+1,j+1,i  )+  45.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk+1,fj,  fi+1)=(
              +  45.*src(n,k-1,j-1,i-1)-  450.*src(n,k-1,j-1,i  )-  75.*src(n,k-1,j-1,i+1)
              + 270.*src(n,k-1,j,  i-1)- 2700.*src(n,k-1,j,  i  )- 450.*src(n,k-1,j,  i+1)
              -  27.*src(n,k-1,j+1,i-1)+  270.*src(n,k-1,j+1,i  )+  45.*src(n,k-1,j+1,i+1)
              - 450.*src(n,k,  j-1,i-1)+ 4500.*src(n,k,  j-1,i  )+ 750.*src(n,k,  j-1,i+1)
              -2700.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )+4500.*src(n,k,  j,  i+1)
              + 270.*src(n,k,  j+1,i-1)- 2700.*src(n,k,  j+1,i  )- 450.*src(n,k,  j+1,i+1)
              -  75.*src(n,k+1,j-1,i-1)+  750.*src(n,k+1,j-1,i  )+ 125.*src(n,k+1,j-1,i+1)
              - 450.*src(n,k+1,j,  i-1)+ 4500.*src(n,k+1,j,  i  )+ 750.*src(n,k+1,j,  i+1)
              +  45.*src(n,k+1,j+1,i-1)-  450.*src(n,k+1,j+1,i  )-  75.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk+1,fj+1,fi  )=(
              +  45.*src(n,k-1,j-1,i-1)+  270.*src(n,k-1,j-1,i  )-  27.*src(n,k-1,j-1,i+1)
              - 450.*src(n,k-1,j,  i-1)- 2700.*src(n,k-1,j,  i  )+ 270.*src(n,k-1,j,  i+1)
              -  75.*src(n,k-1,j+1,i-1)-  450.*src(n,k-1,j+1,i  )+  45.*src(n,k-1,j+1,i+1)
              - 450.*src(n,k,  j-1,i-1)- 2700.*src(n,k,  j-1,i  )+ 270.*src(n,k,  j-1,i+1)
              +4500.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )-2700.*src(n,k,  j,  i+1)
              + 750.*src(n,k,  j+1,i-1)+ 4500.*src(n,k,  j+1,i  )- 450.*src(n,k,  j+1,i+1)
              -  75.*src(n,k+1,j-1,i-1)-  450.*src(n,k+1,j-1,i  )+  45.*src(n,k+1,j-1,i+1)
              + 750.*src(n,k+1,j,  i-1)+ 4500.*src(n,k+1,j,  i  )- 450.*src(n,k+1,j,  i+1)
              + 125.*src(n,k+1,j+1,i-1)+  750.*src(n,k+1,j+1,i  )-  75.*src(n,k+1,j+1,i+1)
                                 )/32768.0;
          dst(n,fk+1,fj+1,fi+1)=(
              -  27.*src(n,k-1,j-1,i-1)+  270.*src(n,k-1,j-1,i  )+  45.*src(n,k-1,j-1,i+1)
              + 270.*src(n,k-1,j,  i-1)- 2700.*src(n,k-1,j,  i  )- 450.*src(n,k-1,j,  i+1)
              +  45.*src(n,k-1,j+1,i-1)-  450.*src(n,k-1,j+1,i  )-  75.*src(n,k-1,j+1,i+1)
              + 270.*src(n,k,  j-1,i-1)- 2700.*src(n,k,  j-1,i  )- 450.*src(n,k,  j-1,i+1)
              -2700.*src(n,k,  j,  i-1)+27000.*src(n,k,  j,  i  )+4500.*src(n,k,  j,  i+1)
              - 450.*src(n,k,  j+1,i-1)+ 4500.*src(n,k,  j+1,i  )+ 750.*src(n,k,  j+1,i+1)
              +  45.*src(n,k+1,j-1,i-1)-  450.*src(n,k+1,j-1,i  )-  75.*src(n,k+1,j-1,i+1)
              - 450.*src(n,k+1,j,  i-1)+ 4500.*src(n,k+1,j,  i  )+ 750.*src(n,k+1,j,  i+1)
              -  75.*src(n,k+1,j+1,i-1)+  750.*src(n,k+1,j+1,i  )+ 125.*src(n,k+1,j+1,i+1)
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
void Multigrid::SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck) {
  current_level_=0;
  AthenaArray<Real> &dst=u_[current_level_];
  for (int n=0; n<nvar_; n++) {
    for (int k=-1; k<=1; k++) {
      for (int j=-1; j<=1; j++) {
        for (int i=-1; i<=1; i++)
          dst(n,ngh_+k,ngh_+j,ngh_+i)=src(n,ck+k+ngh_,cj+j+ngh_,ci+i+ngh_);
      }
    }
  }
  return;
}




//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(int n, int nrm) {
  AthenaArray<Real> &def=def_[current_level_];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  Real norm=0.0;
  if (nrm==0) { // special case: max norm
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          norm=std::max(norm,std::fabs(def(n,k,j,i)));
      }
    }
  } else if (nrm==1) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          norm+=std::fabs(def(n,k,j,i));
      }
    }
  } else { // nrm>1 -> nrm=2
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          norm+=SQR(def(n,k,j,i));
      }
    }
  }
  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotal(int type, int n)
//  \brief calculate the sum of the array (type: 0=src, 1=u)

Real Multigrid::CalculateTotal(int type, int n) {
  AthenaArray<Real> src;
  int ll=nlevel_-1-current_level_;
  if (type==0)
    src.InitWithShallowCopy(src_[current_level_]);
  else
    src.InitWithShallowCopy(u_[current_level_]);
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  Real dx=rdx_*static_cast<Real>(1<<ll), dy=rdy_*static_cast<Real>(1<<ll),
       dz=rdz_*static_cast<Real>  (1<<ll);
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++)
        s+=src(n,k,j,i);
    }
  }
  return s*dx*dy*dz;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverage(int type, int n, Real ave)
//  \brief subtract the average value (type: 0=source, 1=u)

void Multigrid::SubtractAverage(int type, int n, Real ave) {
  AthenaArray<Real> dst;
  if (type==0)
    dst.InitWithShallowCopy(src_[nlevel_-1]);
  else
    dst.InitWithShallowCopy(u_[nlevel_-1]);
  int is, ie, js, je, ks, ke;
  is=js=ks=0;
  ie=is+size_.nx1+1, je=js+size_.nx2+1, ke=ks+size_.nx3+1;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++)
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
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
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
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
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
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
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
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
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
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
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
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i)=dst(n,ks+k,j,i);
      }
    }
  }
  return;
}
