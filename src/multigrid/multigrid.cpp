//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid.cpp
//  \brief implementation of the functions commonly used in Multigrid

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>    // memset, memcpy
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
//! \fn Multigrid::Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int invar, int nghost)
//  \brief Multigrid constructor

Multigrid::Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int invar, int nghost) :
  pmy_driver_(pmd), pmy_block_(pmb), ngh_(nghost), nvar_(invar), defscale_(1.0) {
  if (pmy_block_ != nullptr) {
    loc_ = pmy_block_->loc;
    size_ = pmy_block_->block_size;
    if (size_.nx1 != size_.nx2 || size_.nx1 != size_.nx3) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Multigrid::Multigrid" << std::endl
          << "The Multigrid solver requires logically cubic MeshBlock." << std::endl;
      ATHENA_ERROR(msg);
      return;
    }
  } else {
    loc_.lx1 = loc_.lx2 = loc_.lx3 = 0;
    loc_.level = 0;
    size_ = pmy_driver_->pmy_mesh_->mesh_size;
    size_.nx1 = pmy_driver_->nrbx1_;
    size_.nx2 = pmy_driver_->nrbx2_;
    size_.nx3 = pmy_driver_->nrbx3_;
  }

  rdx_=(size_.x1max-size_.x1min)/static_cast<Real>(size_.nx1);
  rdy_=(size_.x2max-size_.x2min)/static_cast<Real>(size_.nx2);
  rdz_=(size_.x3max-size_.x3min)/static_cast<Real>(size_.nx3);

  nlevel_=0;
  if (pmy_block_ == nullptr) { // root
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
    if (std::abs(rdx_-rdy_)>1.0e-5 || std::abs(rdx_-rdz_)>1.0e-5) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Multigrid::Multigrid" << std::endl
          << "The cell size must be cubic." << std::endl;
      ATHENA_ERROR(msg);
      return;
    }
  }

  current_level_ = nlevel_-1;

  // allocate arrays
  u_ = new AthenaArray<Real>[nlevel_];
  src_ = new AthenaArray<Real>[nlevel_];
  def_ = new AthenaArray<Real>[nlevel_];
  if (pmy_driver_->ffas_) {
    if (pmy_block_ == nullptr)
      uold_ = new AthenaArray<Real>[nlevel_];
    else
      uold_ = new AthenaArray<Real>[nlevel_-1];
  }
  for (int l=0; l<nlevel_; l++) {
    int ll=nlevel_-1-l;
    int ncx=(size_.nx1>>ll)+2*ngh_, ncy=(size_.nx2>>ll)+2*ngh_,
        ncz=(size_.nx3>>ll)+2*ngh_;
    u_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    src_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    if (pmy_driver_->ffas_ && !((pmy_block_ != nullptr) && (l == nlevel_-1)))
      uold_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
  }
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid::~Multigrid
//  \brief Multigrid destroctor

Multigrid::~Multigrid() {
  delete [] u_;
  delete [] src_;
  delete [] def_;
  if (pmy_driver_->ffas_) delete [] uold_;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the inital guess in the active zone of the finest level

void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh) {
  AthenaArray<Real> &dst=u_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  for (int v=0; v<nvar_; ++v) {
    int nsrc=ns+v;
    for (int k=ngh, mk=ks; mk<=ke; ++k, ++mk) {
      for (int j=ngh, mj=js; mj<=je; ++j, ++mj) {
        for (int i=ngh, mi=is; mi<=ie; ++i, ++mi)
          dst(v,mk,mj,mi)=src(nsrc,k,j,i);
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
    for (int v=0; v<nvar_; ++v) {
      int nsrc=ns+v;
      for (int k=ngh, mk=ks; mk<=ke; ++k, ++mk) {
        for (int j=ngh, mj=js; mj<=je; ++j, ++mj) {
          for (int i=ngh, mi=is; mi<=ie; ++i, ++mi)
            dst(v,mk,mj,mi)=src(nsrc,k,j,i);
        }
      }
    }
  } else {
    for (int v=0; v<nvar_; ++v) {
      int nsrc=ns+v;
      for (int k=ngh, mk=ks; mk<=ke; ++k, ++mk) {
        for (int j=ngh, mj=js; mj<=je; ++j, ++mj) {
          for (int i=ngh, mi=is; mi<=ie; ++i, ++mi)
            dst(v,mk,mj,mi)=src(nsrc,k,j,i)*fac;
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
  for (current_level_=nlevel_-1; current_level_>0; current_level_--) {
    int ll=nlevel_-current_level_;
    ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
    Restrict(src_[current_level_-1], src_[current_level_], is, ie, js, je, ks, ke);
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
  for (int v=0; v<nvar_; ++v) {
    int ndst=ns+v;
    for (int k=ngh-sngh, mk=ngh_-sngh; mk<=ke; ++k, ++mk) {
      for (int j=ngh-sngh, mj=ngh_-sngh; mj<=je; ++j, ++mj) {
        for (int i=ngh-sngh, mi=ngh_-sngh; mi<=ie; ++i, ++mi)
          dst(ndst,k,j,i)=src(v,mk,mj,mi);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the defect, including the ghost zone

void Multigrid::RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh) {
  const AthenaArray<Real> &src=def_[nlevel_-1];
  int sngh=std::min(ngh_,ngh);
  int ie=size_.nx1+ngh_+sngh-1, je=size_.nx2+ngh_+sngh-1, ke=size_.nx3+ngh_+sngh-1;
  for (int v=0; v<nvar_; ++v) {
    int ndst=ns+v;
    for (int k=ngh-sngh, mk=ngh_-sngh; mk<=ke; ++k, ++mk) {
      for (int j=ngh-sngh, mj=ngh_-sngh; mj<=je; ++j, ++mj) {
        for (int i=ngh-sngh, mi=ngh_-sngh; mi<=ie; ++i, ++mi)
          dst(ndst,k,j,i)=src(v,mk,mj,mi)*defscale_;
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ZeroClearData()
//  \brief Clear the data array with zero

void Multigrid::ZeroClearData() {
  u_[current_level_].ZeroClear();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictBlock()
//  \brief Restrict the defect to the source

void Multigrid::RestrictBlock() {
  int ll=nlevel_-current_level_;
  int is, ie, js, je, ks, ke;

  CalculateDefectBlock();
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  Restrict(src_[current_level_-1], def_[current_level_], is, ie, js, je, ks, ke);

  // Full Approximation Scheme - restrict the variable itself
  if (pmy_driver_->ffas_)
    Restrict(u_[current_level_-1], u_[current_level_], is, ie, js, je, ks, ke);

  current_level_--;

  if (!pmy_driver_->ffas_) ZeroClearData();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrectBlock()
//  \brief Prolongate the potential using tri-linear interpolation

void Multigrid::ProlongateAndCorrectBlock() {
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  if (pmy_driver_->ffas_) {
    int size = u_[current_level_].GetSize();
    for (int s=0; s<size; ++s)
      u_[current_level_](s) -= uold_[current_level_](s);
  }

  ProlongateAndCorrect(u_[current_level_+1], u_[current_level_], is, ie, js, je, ks, ke,
                                                                 ngh_, ngh_, ngh_);

  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongateBlock()
//  \brief Prolongate the potential for Full Multigrid cycle

void Multigrid::FMGProlongateBlock() {
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  FMGProlongate(u_[current_level_+1], u_[current_level_], is, ie, js, je, ks, ke,
                                                          ngh_, ngh_, ngh_);

  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::SmoothBlock(int color)
//  \brief Apply Smoother on the Block

void Multigrid::SmoothBlock(int color) {
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;

  Smooth(u_[current_level_], src_[current_level_], -ll, is, ie, js, je, ks, ke, color);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateDefectBlock()
//  \brief calculate the residual

void Multigrid::CalculateDefectBlock() {
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;

  CalculateDefect(def_[current_level_], u_[current_level_], src_[current_level_],
                  -ll, is, ie, js, je, ks, ke);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateFASRHSBlock()
//  \brief calculate the RHS for the Full Approximation Scheme

void Multigrid::CalculateFASRHSBlock() {
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;

  CalculateFASRHS(src_[current_level_], u_[current_level_], -ll, is, ie, js, je, ks, ke);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetFromRootGrid(bool folddata)
//  \brief Load the data from the root grid or octets

void Multigrid::SetFromRootGrid(bool folddata) {
  current_level_=0;
  AthenaArray<Real> &dst = u_[current_level_];
  int lev = loc_.level - pmy_driver_->locrootlevel_;
  if (lev == 0) { // from the root grid
    int ci = static_cast<int>(loc_.lx1);
    int cj = static_cast<int>(loc_.lx2);
    int ck = static_cast<int>(loc_.lx3);
    const AthenaArray<Real> &src=pmy_driver_->mgroot_->GetCurrentData();
    for (int v=0; v<nvar_; ++v) {
      for (int k=0; k<=2; ++k) {
        for (int j=0; j<=2; ++j) {
          for (int i=0; i<=2; ++i)
            dst(v, k, j, i) = src(v, ck+k, cj+j, ci+i);
        }
      }
    }
    if (folddata) {
      AthenaArray<Real> &odst = uold_[current_level_];
      const AthenaArray<Real> &osrc = pmy_driver_->mgroot_->GetCurrentOldData();
      for (int v=0; v<nvar_; ++v) {
        for (int k=0; k<=2; ++k) {
          for (int j=0; j<=2; ++j) {
            for (int i=0; i<=2; ++i)
              odst(v, k, j, i) = osrc(v, ck+k, cj+j, ci+i);
          }
        }
      }
    }
  } else { // from an octet
    LogicalLocation oloc;
    oloc.lx1 = (loc_.lx1 >> 1);
    oloc.lx2 = (loc_.lx2 >> 1);
    oloc.lx3 = (loc_.lx3 >> 1);
    oloc.level = loc_.level - 1;
    int olev = oloc.level - pmy_driver_->locrootlevel_;
    int oid = pmy_driver_->octetmap_[olev][oloc];
    int ci = (static_cast<int>(loc_.lx1)&1);
    int cj = (static_cast<int>(loc_.lx2)&1);
    int ck = (static_cast<int>(loc_.lx3)&1);
    const AthenaArray<Real> &src = pmy_driver_->octets_[olev][oid].u;
    for (int v=0; v<nvar_; ++v) {
      for (int k=0; k<=2; ++k) {
        for (int j=0; j<=2; ++j) {
          for (int i=0; i<=2; ++i)
            dst(v, k, j, i)=src(v, ck+k, cj+j, ci+i);
        }
      }
    }
    if (folddata) {
      AthenaArray<Real> &odst = uold_[current_level_];
      const AthenaArray<Real> &osrc = pmy_driver_->octets_[olev][oid].uold;
      for (int v=0; v<nvar_; ++v) {
        for (int k=0; k<=2; ++k) {
          for (int j=0; j<=2; ++j) {
            for (int i=0; i<=2; ++i)
              odst(v, k, j, i)=osrc(v, ck+k, cj+j, ci+i);
          }
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(MGNormType nrm, int n)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(MGNormType nrm, int n) {
  AthenaArray<Real> &def=def_[current_level_];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  Real dx=rdx_*static_cast<Real>(1<<ll), dy=rdy_*static_cast<Real>(1<<ll),
       dz=rdz_*static_cast<Real>(1<<ll);

  CalculateDefect(def_[current_level_], u_[current_level_], src_[current_level_],
                  -ll, is, ie, js, je, ks, ke);

  Real norm=0.0;
  if (nrm == MGNormType::max) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i)
          norm=std::max(norm,std::fabs(def(n,k,j,i)));
      }
    }
    return norm;
  } else if (nrm == MGNormType::l1) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i)
          norm+=std::fabs(def(n,k,j,i));
      }
    }
  } else { // L2 norm
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i)
          norm+=SQR(def(n,k,j,i));
      }
    }
  }
  return norm*dx*dy*dz*defscale_;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotal(MGVariable type, int n)
//  \brief calculate the sum of the array (type: 0=src, 1=u)

Real Multigrid::CalculateTotal(MGVariable type, int n) {
  AthenaArray<Real> &src =
                    (type == MGVariable::src) ? src_[current_level_] : u_[current_level_];
  int ll = nlevel_ - 1 - current_level_;
  Real s=0.0;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  Real dx=rdx_*static_cast<Real>(1<<ll), dy=rdy_*static_cast<Real>(1<<ll),
       dz=rdz_*static_cast<Real>(1<<ll);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i)
        s+=src(n,k,j,i);
    }
  }
  return s*dx*dy*dz;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverage(MGVariable type, int v, Real ave)
//  \brief subtract the average value (type: 0=src, 1=u)

void Multigrid::SubtractAverage(MGVariable type, int n, Real ave) {
  AthenaArray<Real> &dst = (type == MGVariable::src) ? src_[nlevel_-1] : u_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=0;
  ie=is+size_.nx1+1, je=js+size_.nx2+1, ke=ks+size_.nx3+1;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i)
        dst(n,k,j,i)-=ave;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::StoreOldData()
//  \brief store the old u data in the uold array

void Multigrid::StoreOldData() {
  memcpy(uold_[current_level_].data(), u_[current_level_].data(),
         u_[current_level_].GetSizeInBytes());

  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::GetCoarsestData(MGVariable type, int n)
//  \brief get the value on the coarsest level in the MG block (type: 0=src, 1=u)

Real Multigrid::GetCoarsestData(MGVariable type, int n) {
  AthenaArray<Real> &src = (type == MGVariable::src) ? src_[0] : u_[0];
  return src(n, ngh_, ngh_, ngh_);
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetData(MGVariable type, int n, int k, int j, int i, Real v)
//  \brief set a value to a cell on the current level

void Multigrid::SetData(MGVariable type, int n, int k, int j, int i, Real v) {
  AthenaArray<Real> &dst =
                    (type == MGVariable::src) ? src_[current_level_] : u_[current_level_];
  dst(n, ngh_+k, ngh_+j, ngh_+i) = v;

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
//                               int il, int iu, int jl, int ju, int kl, int ku)
//  \brief Actual implementation of prolongation and correction

void Multigrid::Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                         int il, int iu, int jl, int ju, int kl, int ku) {
  for (int v=0; v<nvar_; ++v) {
    for (int k=kl, fk=kl; k<=ku; ++k, fk+=2) {
      for (int j=jl, fj=jl; j<=ju; ++j, fj+=2) {
        for (int i=il, fi=il; i<=iu; ++i, fi+=2)
          dst(v, k, j, i)=0.125*(src(v, fk,   fj,   fi)+src(v, fk,   fj,   fi+1)
                                +src(v, fk,   fj+1, fi)+src(v, fk,   fj+1, fi+1)
                                +src(v, fk+1, fj,   fi)+src(v, fk+1, fj,   fi+1)
                                +src(v, fk+1, fj+1, fi)+src(v, fk+1, fj+1, fi+1));
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(AthenaArray<Real> &dst,
//      const AthenaArray<Real> &src, int il, int iu, int jl, int ju, int kl, int ku,
//      int fil, int fjl, int fkl)
//  \brief Actual implementation of prolongation and correction

void Multigrid::ProlongateAndCorrect(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
     int il, int iu, int jl, int ju, int kl, int ku, int fil, int fjl, int fkl) {
  for (int v=0; v<nvar_; ++v) {
    for (int k=kl, fk=fkl; k<=ku; ++k, fk+=2) {
      for (int j=jl, fj=fjl; j<=ju; ++j, fj+=2) {
        for (int i=il, fi=fil; i<=iu; ++i, fi+=2) {
          dst(v,fk  ,fj  ,fi  ) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k-1,j-1,i-1)
                        +9.0*(src(v,k,j,i-1)+src(v,k,j-1,i)+src(v,k-1,j,i))
                        +3.0*(src(v,k-1,j-1,i)+src(v,k-1,j,i-1)+src(v,k,j-1,i-1)));
          dst(v,fk  ,fj  ,fi+1) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k-1,j-1,i+1)
                        +9.0*(src(v,k,j,i+1)+src(v,k,j-1,i)+src(v,k-1,j,i))
                        +3.0*(src(v,k-1,j-1,i)+src(v,k-1,j,i+1)+src(v,k,j-1,i+1)));
          dst(v,fk  ,fj+1,fi  ) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k-1,j+1,i-1)
                        +9.0*(src(v,k,j,i-1)+src(v,k,j+1,i)+src(v,k-1,j,i))
                        +3.0*(src(v,k-1,j+1,i)+src(v,k-1,j,i-1)+src(v,k,j+1,i-1)));
          dst(v,fk+1,fj  ,fi  ) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k+1,j-1,i-1)
                        +9.0*(src(v,k,j,i-1)+src(v,k,j-1,i)+src(v,k+1,j,i))
                        +3.0*(src(v,k+1,j-1,i)+src(v,k+1,j,i-1)+src(v,k,j-1,i-1)));
          dst(v,fk+1,fj+1,fi  ) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k+1,j+1,i-1)
                        +9.0*(src(v,k,j,i-1)+src(v,k,j+1,i)+src(v,k+1,j,i))
                        +3.0*(src(v,k+1,j+1,i)+src(v,k+1,j,i-1)+src(v,k,j+1,i-1)));
          dst(v,fk+1,fj  ,fi+1) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k+1,j-1,i+1)
                        +9.0*(src(v,k,j,i+1)+src(v,k,j-1,i)+src(v,k+1,j,i))
                        +3.0*(src(v,k+1,j-1,i)+src(v,k+1,j,i+1)+src(v,k,j-1,i+1)));
          dst(v,fk  ,fj+1,fi+1) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k-1,j+1,i+1)
                        +9.0*(src(v,k,j,i+1)+src(v,k,j+1,i)+src(v,k-1,j,i))
                        +3.0*(src(v,k-1,j+1,i)+src(v,k-1,j,i+1)+src(v,k,j+1,i+1)));
          dst(v,fk+1,fj+1,fi+1) +=
              0.015625*(27.0*src(v,k,j,i) + src(v,k+1,j+1,i+1)
                        +9.0*(src(v,k,j,i+1)+src(v,k,j+1,i)+src(v,k+1,j,i))
                        +3.0*(src(v,k+1,j+1,i)+src(v,k+1,j,i+1)+src(v,k,j+1,i+1)));
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongate(AthenaArray<Real> &dst,
//           const AthenaArray<Real> &src, int il, int iu, int jl, int ju, int kl, int ku
//           int fil, int fjl, int fkl)
//  \brief Actual implementation of FMG prolongation

void Multigrid::FMGProlongate(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                              int il, int iu, int jl, int ju, int kl, int ku,
                              int fil, int fjl, int fkl) {
  for (int v=0; v<nvar_; ++v) {
    for (int k=kl, fk=fkl; k<=ku; ++k, fk+=2) {
      for (int j=jl, fj=fjl; j<=ju; ++j, fj+=2) {
        for (int i=il, fi=fil; i<=iu; ++i, fi+=2) {
          dst(v,fk  ,fj,  fi  )=(
              + 125.*src(v,k-1,j-1,i-1)+  750.*src(v,k-1,j-1,i  )-  75.*src(v,k-1,j-1,i+1)
              + 750.*src(v,k-1,j,  i-1)+ 4500.*src(v,k-1,j,  i  )- 450.*src(v,k-1,j,  i+1)
              -  75.*src(v,k-1,j+1,i-1)-  450.*src(v,k-1,j+1,i  )+  45.*src(v,k-1,j+1,i+1)
              + 750.*src(v,k,  j-1,i-1)+ 4500.*src(v,k,  j-1,i  )- 450.*src(v,k,  j-1,i+1)
              +4500.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )-2700.*src(v,k,  j,  i+1)
              - 450.*src(v,k,  j+1,i-1)- 2700.*src(v,k,  j+1,i  )+ 270.*src(v,k,  j+1,i+1)
              -  75.*src(v,k+1,j-1,i-1)-  450.*src(v,k+1,j-1,i  )+  45.*src(v,k+1,j-1,i+1)
              - 450.*src(v,k+1,j,  i-1)- 2700.*src(v,k+1,j,  i  )+ 270.*src(v,k+1,j,  i+1)
              +  45.*src(v,k+1,j+1,i-1)+  270.*src(v,k+1,j+1,i  )-  27.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk,  fj,  fi+1)=(
              -  75.*src(v,k-1,j-1,i-1)+  750.*src(v,k-1,j-1,i  )+ 125.*src(v,k-1,j-1,i+1)
              - 450.*src(v,k-1,j,  i-1)+ 4500.*src(v,k-1,j,  i  )+ 750.*src(v,k-1,j,  i+1)
              +  45.*src(v,k-1,j+1,i-1)-  450.*src(v,k-1,j+1,i  )-  75.*src(v,k-1,j+1,i+1)
              - 450.*src(v,k,  j-1,i-1)+ 4500.*src(v,k,  j-1,i  )+ 750.*src(v,k,  j-1,i+1)
              -2700.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )+4500.*src(v,k,  j,  i+1)
              + 270.*src(v,k,  j+1,i-1)- 2700.*src(v,k,  j+1,i  )- 450.*src(v,k,  j+1,i+1)
              +  45.*src(v,k+1,j-1,i-1)-  450.*src(v,k+1,j-1,i  )-  75.*src(v,k+1,j-1,i+1)
              + 270.*src(v,k+1,j,  i-1)- 2700.*src(v,k+1,j,  i  )- 450.*src(v,k+1,j,  i+1)
              -  27.*src(v,k+1,j+1,i-1)+  270.*src(v,k+1,j+1,i  )+  45.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk  ,fj+1,fi  )=(
              -  75.*src(v,k-1,j-1,i-1)-  450.*src(v,k-1,j-1,i  )+  45.*src(v,k-1,j-1,i+1)
              + 750.*src(v,k-1,j,  i-1)+ 4500.*src(v,k-1,j,  i  )- 450.*src(v,k-1,j,  i+1)
              + 125.*src(v,k-1,j+1,i-1)+  750.*src(v,k-1,j+1,i  )-  75.*src(v,k-1,j+1,i+1)
              - 450.*src(v,k,  j-1,i-1)- 2700.*src(v,k,  j-1,i  )+ 270.*src(v,k,  j-1,i+1)
              +4500.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )-2700.*src(v,k,  j,  i+1)
              + 750.*src(v,k,  j+1,i-1)+ 4500.*src(v,k,  j+1,i  )- 450.*src(v,k,  j+1,i+1)
              +  45.*src(v,k+1,j-1,i-1)+  270.*src(v,k+1,j-1,i  )-  27.*src(v,k+1,j-1,i+1)
              - 450.*src(v,k+1,j,  i-1)- 2700.*src(v,k+1,j,  i  )+ 270.*src(v,k+1,j,  i+1)
              -  75.*src(v,k+1,j+1,i-1)-  450.*src(v,k+1,j+1,i  )+  45.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk,  fj+1,fi+1)=(
              +  45.*src(v,k-1,j-1,i-1)-  450.*src(v,k-1,j-1,i  )-  75.*src(v,k-1,j-1,i+1)
              - 450.*src(v,k-1,j,  i-1)+ 4500.*src(v,k-1,j,  i  )+ 750.*src(v,k-1,j,  i+1)
              -  75.*src(v,k-1,j+1,i-1)+  750.*src(v,k-1,j+1,i  )+ 125.*src(v,k-1,j+1,i+1)
              + 270.*src(v,k,  j-1,i-1)- 2700.*src(v,k,  j-1,i  )- 450.*src(v,k,  j-1,i+1)
              -2700.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )+4500.*src(v,k,  j,  i+1)
              - 450.*src(v,k,  j+1,i-1)+ 4500.*src(v,k,  j+1,i  )+ 750.*src(v,k,  j+1,i+1)
              -  27.*src(v,k+1,j-1,i-1)+  270.*src(v,k+1,j-1,i  )+  45.*src(v,k+1,j-1,i+1)
              + 270.*src(v,k+1,j,  i-1)- 2700.*src(v,k+1,j,  i  )- 450.*src(v,k+1,j,  i+1)
              +  45.*src(v,k+1,j+1,i-1)-  450.*src(v,k+1,j+1,i  )-  75.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk+1,fj,  fi  )=(
              -  75.*src(v,k-1,j-1,i-1)-  450.*src(v,k-1,j-1,i  )+  45.*src(v,k-1,j-1,i+1)
              - 450.*src(v,k-1,j,  i-1)- 2700.*src(v,k-1,j,  i  )+ 270.*src(v,k-1,j,  i+1)
              +  45.*src(v,k-1,j+1,i-1)+  270.*src(v,k-1,j+1,i  )-  27.*src(v,k-1,j+1,i+1)
              + 750.*src(v,k,  j-1,i-1)+ 4500.*src(v,k,  j-1,i  )- 450.*src(v,k,  j-1,i+1)
              +4500.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )-2700.*src(v,k,  j,  i+1)
              - 450.*src(v,k,  j+1,i-1)- 2700.*src(v,k,  j+1,i  )+ 270.*src(v,k,  j+1,i+1)
              + 125.*src(v,k+1,j-1,i-1)+  750.*src(v,k+1,j-1,i  )-  75.*src(v,k+1,j-1,i+1)
              + 750.*src(v,k+1,j,  i-1)+ 4500.*src(v,k+1,j,  i  )- 450.*src(v,k+1,j,  i+1)
              -  75.*src(v,k+1,j+1,i-1)-  450.*src(v,k+1,j+1,i  )+  45.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk+1,fj,  fi+1)=(
              +  45.*src(v,k-1,j-1,i-1)-  450.*src(v,k-1,j-1,i  )-  75.*src(v,k-1,j-1,i+1)
              + 270.*src(v,k-1,j,  i-1)- 2700.*src(v,k-1,j,  i  )- 450.*src(v,k-1,j,  i+1)
              -  27.*src(v,k-1,j+1,i-1)+  270.*src(v,k-1,j+1,i  )+  45.*src(v,k-1,j+1,i+1)
              - 450.*src(v,k,  j-1,i-1)+ 4500.*src(v,k,  j-1,i  )+ 750.*src(v,k,  j-1,i+1)
              -2700.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )+4500.*src(v,k,  j,  i+1)
              + 270.*src(v,k,  j+1,i-1)- 2700.*src(v,k,  j+1,i  )- 450.*src(v,k,  j+1,i+1)
              -  75.*src(v,k+1,j-1,i-1)+  750.*src(v,k+1,j-1,i  )+ 125.*src(v,k+1,j-1,i+1)
              - 450.*src(v,k+1,j,  i-1)+ 4500.*src(v,k+1,j,  i  )+ 750.*src(v,k+1,j,  i+1)
              +  45.*src(v,k+1,j+1,i-1)-  450.*src(v,k+1,j+1,i  )-  75.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk+1,fj+1,fi  )=(
              +  45.*src(v,k-1,j-1,i-1)+  270.*src(v,k-1,j-1,i  )-  27.*src(v,k-1,j-1,i+1)
              - 450.*src(v,k-1,j,  i-1)- 2700.*src(v,k-1,j,  i  )+ 270.*src(v,k-1,j,  i+1)
              -  75.*src(v,k-1,j+1,i-1)-  450.*src(v,k-1,j+1,i  )+  45.*src(v,k-1,j+1,i+1)
              - 450.*src(v,k,  j-1,i-1)- 2700.*src(v,k,  j-1,i  )+ 270.*src(v,k,  j-1,i+1)
              +4500.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )-2700.*src(v,k,  j,  i+1)
              + 750.*src(v,k,  j+1,i-1)+ 4500.*src(v,k,  j+1,i  )- 450.*src(v,k,  j+1,i+1)
              -  75.*src(v,k+1,j-1,i-1)-  450.*src(v,k+1,j-1,i  )+  45.*src(v,k+1,j-1,i+1)
              + 750.*src(v,k+1,j,  i-1)+ 4500.*src(v,k+1,j,  i  )- 450.*src(v,k+1,j,  i+1)
              + 125.*src(v,k+1,j+1,i-1)+  750.*src(v,k+1,j+1,i  )-  75.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
          dst(v,fk+1,fj+1,fi+1)=(
              -  27.*src(v,k-1,j-1,i-1)+  270.*src(v,k-1,j-1,i  )+  45.*src(v,k-1,j-1,i+1)
              + 270.*src(v,k-1,j,  i-1)- 2700.*src(v,k-1,j,  i  )- 450.*src(v,k-1,j,  i+1)
              +  45.*src(v,k-1,j+1,i-1)-  450.*src(v,k-1,j+1,i  )-  75.*src(v,k-1,j+1,i+1)
              + 270.*src(v,k,  j-1,i-1)- 2700.*src(v,k,  j-1,i  )- 450.*src(v,k,  j-1,i+1)
              -2700.*src(v,k,  j,  i-1)+27000.*src(v,k,  j,  i  )+4500.*src(v,k,  j,  i+1)
              - 450.*src(v,k,  j+1,i-1)+ 4500.*src(v,k,  j+1,i  )+ 750.*src(v,k,  j+1,i+1)
              +  45.*src(v,k+1,j-1,i-1)-  450.*src(v,k+1,j-1,i  )-  75.*src(v,k+1,j-1,i+1)
              - 450.*src(v,k+1,j,  i-1)+ 4500.*src(v,k+1,j,  i  )+ 750.*src(v,k+1,j,  i+1)
              -  75.*src(v,k+1,j+1,i-1)+  750.*src(v,k+1,j+1,i  )+ 125.*src(v,k+1,j+1,i+1)
                                 )/32768.0;
        }
      }
    }
  }

  return;
}

