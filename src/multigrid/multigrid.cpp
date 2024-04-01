//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid.cpp
//! \brief implementation of the functions commonly used in Multigrid

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
//! \fn Multigrid::Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int nghost)
//  \brief Multigrid constructor

Multigrid::Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int nghost) :
  pmy_driver_(pmd), pmy_block_(pmb), ngh_(nghost), nvar_(pmd->nvar_),
  ncoeff_(pmd->ncoeff_), nmatrix_(pmd->nmatrix_), defscale_(1.0) {
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
    for (int i = 0; i < 6; ++i) {
      if (pmy_block_->pbval->block_bcs[i] == BoundaryFlag::block)
        mg_block_bcs_[i] = BoundaryFlag::block;
      else
        mg_block_bcs_[i] = pmy_driver_->mg_mesh_bcs_[i];
    }
  } else {
    loc_.lx1 = loc_.lx2 = loc_.lx3 = 0;
    loc_.level = 0;
    size_ = pmy_driver_->pmy_mesh_->mesh_size;
    size_.nx1 = pmy_driver_->nrbx1_;
    size_.nx2 = pmy_driver_->nrbx2_;
    size_.nx3 = pmy_driver_->nrbx3_;
    for (int i = 0; i < 6; ++i)
      mg_block_bcs_[i] = pmy_driver_->mg_mesh_bcs_[i];
  }
  rdx_ = (size_.x1max-size_.x1min)/static_cast<Real>(size_.nx1);
  rdy_ = (size_.x2max-size_.x2min)/static_cast<Real>(size_.nx2);
  rdz_ = (size_.x3max-size_.x3min)/static_cast<Real>(size_.nx3);

  nlevel_ = 0;
  if (pmy_block_ == nullptr) { // root
    int nbx = 0, nby = 0, nbz = 0;
    for (int l = 0; l < 20; l++) {
      if (size_.nx1%(1<<l) == 0 && size_.nx2%(1<<l) == 0 && size_.nx3%(1<<l) == 0) {
        nbx = size_.nx1/(1<<l), nby = size_.nx2/(1<<l), nbz = size_.nx3/(1<<l);
        nlevel_ = l+1;
      }
    }
    int nmaxr = std::max(nbx, std::max(nby, nbz));
    // int nminr=std::min(nbx, std::min(nby, nbz)); // unused variable
    if (nmaxr != 1 && Globals::my_rank == 0) {
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
    for (int l = 0; l < 20; l++) {
      if ((1<<l) == size_.nx1) {
        nlevel_=l+1;
        break;
      }
    }
    if (nlevel_ == 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Multigrid::Multigrid" << std::endl
          << "The MeshBlock size must be power of two." << std::endl;
      ATHENA_ERROR(msg);
      return;
    }
    // *** temporary ***
    if (std::abs(rdx_-rdy_) > 1.0e-5 || std::abs(rdx_-rdz_) > 1.0e-5) {
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
  coord_ = new MGCoordinates[nlevel_];
  ccoord_ = new MGCoordinates[nlevel_];
  coeff_ = new AthenaArray<Real>[nlevel_];
  matrix_ = new AthenaArray<Real>[nlevel_];
  if (pmy_block_ == nullptr)
    uold_ = new AthenaArray<Real>[nlevel_];
  else
    uold_ = new AthenaArray<Real>[nlevel_];
  for (int l = 0; l < nlevel_; l++) {
    int ll=nlevel_-1-l;
    int ncx=(size_.nx1>>ll)+2*ngh_;
    int ncy=(size_.nx2>>ll)+2*ngh_;
    int ncz=(size_.nx3>>ll)+2*ngh_;
    u_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    src_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    def_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    if (!((pmy_block_ != nullptr) && (l == nlevel_-1)))
      uold_[l].NewAthenaArray(nvar_,ncz,ncy,ncx);
    coord_[l].AllocateMGCoordinates(ncx,ncy,ncz);
    coord_[l].CalculateMGCoordinates(size_, ll, ngh_);
    if (ncoeff_ > 0)
      coeff_[l].NewAthenaArray(ncoeff_,ncz,ncy,ncx);
    if (nmatrix_ > 0)
      matrix_[l].NewAthenaArray(nmatrix_,ncz,ncy,ncx);
    ncx=(size_.nx1>>(ll+1))+2*ngh_;
    ncy=(size_.nx2>>(ll+1))+2*ngh_;
    ncz=(size_.nx3>>(ll+1))+2*ngh_;
    ccoord_[l].AllocateMGCoordinates(ncx,ncy,ncz);
    ccoord_[l].CalculateMGCoordinates(size_, ll+1, ngh_);
  }
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid::~Multigrid
//! \brief Multigrid destroctor

Multigrid::~Multigrid() {
  delete [] u_;
  delete [] src_;
  delete [] def_;
  delete [] uold_;
  delete [] coeff_;
  delete [] matrix_;
  delete [] coord_;
  delete [] ccoord_;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//! \brief Fill the inital guess in the active zone of the finest level

void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh) {
  AthenaArray<Real> &dst=u_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  for (int v=0; v<nvar_; ++v) {
    int nsrc=ns+v;
    for (int mk=ks; mk<=ke; ++mk) {
      int k = mk - ks + ngh;
      for (int mj=js; mj<=je; ++mj) {
        int j = mj - js + ngh;
#pragma omp simd
        for (int mi=is; mi<=ie; ++mi) {
          int i = mi - is + ngh;
          dst(v,mk,mj,mi)=src(nsrc,k,j,i);
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//!                                Real fac)
//! \brief Fill the source in the active zone of the finest level

void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac) {
  AthenaArray<Real> &dst=src_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  if (fac == 1.0) {
    for (int v=0; v<nvar_; ++v) {
      int nsrc=ns+v;
      for (int mk=ks; mk<=ke; ++mk) {
        int k = mk - ks + ngh;
        for (int mj=js; mj<=je; ++mj) {
          int j = mj - js + ngh;
#pragma omp simd
          for (int mi=is; mi<=ie; ++mi) {
            int i = mi - is + ngh;
            dst(v,mk,mj,mi)=src(nsrc,k,j,i);
          }
        }
      }
    }
  } else {
    for (int v=0; v<nvar_; ++v) {
      int nsrc=ns+v;
      for (int mk=ks; mk<=ke; ++mk) {
        int k = mk - ks + ngh;
        for (int mj=js; mj<=je; ++mj) {
          int j = mj - js + ngh;
#pragma omp simd
          for (int mi=is; mi<=ie; ++mi) {
            int i = mi - is + ngh;
            dst(v,mk,mj,mi)=src(nsrc,k,j,i)*fac;
          }
        }
      }
    }
  }
  current_level_ = nlevel_-1;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadCoefficients(const AthenaArray<Real> &coeff, int ngh)
//! \brief Load coefficients of the diffusion and source terms

void Multigrid::LoadCoefficients(const AthenaArray<Real> &coeff, int ngh) {
  AthenaArray<Real> &cm=coeff_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=0;
  ie=size_.nx1+2*ngh_-1, je=size_.nx2+2*ngh_-1, ke=size_.nx3+2*ngh_-1;
  for (int v = 0; v < ncoeff_; ++v) {
    for (int mk=ks; mk<=ke; ++mk) {
      int k = mk + ngh - ngh_;
      for (int mj=js; mj<=je; ++mj) {
        int j = mj + ngh - ngh_;
#pragma omp simd
        for (int mi=is; mi<=ie; ++mi) {
          int i = mi + ngh - ngh_;
          cm(v,mk,mj,mi) = coeff(v,k,j,i);
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ApplyMask()
//  \brief Apply the user-defined source mask function on the finest level

void Multigrid::ApplyMask() {
  int is, ie, js, je, ks, ke;
  is = js = ks = ngh_;
  ie = is + size_.nx1;
  je = js + size_.nx2;
  ke = ks + size_.nx3;
  if (pmy_driver_->srcmask_ != nullptr)
    pmy_driver_->srcmask_(src_[nlevel_-1], is, ie, js, je, ks, ke, coord_[nlevel_-1]);
  if (ncoeff_ > 0 && pmy_driver_->coeffmask_ != nullptr)
    pmy_driver_->coeffmask_(coeff_[nlevel_-1], is, ie, js, je, ks, ke, coord_[nlevel_-1]);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictFMGSource()
//! \brief restrict the source through all the multigrid levels

void Multigrid::RestrictFMGSource() {
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  for (current_level_=nlevel_-1; current_level_>0; current_level_--) {
    int ll=nlevel_-current_level_;
    ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
    Restrict(src_[current_level_-1], src_[current_level_],
             nvar_, is, ie, js, je, ks, ke, false);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictCoefficients()
//! \brief restrict coefficients within a Multigrid object

void Multigrid::RestrictCoefficients() {
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  for (int lev = nlevel_ - 1; lev > 0; lev--) {
    int ll = nlevel_ - lev;
    ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
    Restrict(coeff_[lev-1], coeff_[lev], ncoeff_, is, ie, js, je, ks, ke, false);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
//! \brief Set the result, including the ghost zone

void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh) {
  const AthenaArray<Real> &src=u_[nlevel_-1];
  int sngh=std::min(ngh_,ngh);
  int ie=size_.nx1+ngh_+sngh-1, je=size_.nx2+ngh_+sngh-1, ke=size_.nx3+ngh_+sngh-1;
  for (int v=0; v<nvar_; ++v) {
    int ndst=ns+v;
    for (int mk=ngh_-sngh; mk<=ke; ++mk) {
      int k = mk - ngh_ + ngh;
      for (int mj=ngh_-sngh; mj<=je; ++mj) {
        int j = mj - ngh_ + ngh;
#pragma omp simd
        for (int mi=ngh_-sngh; mi<=ie; ++mi) {
          int i = mi - ngh_ + ngh;
          dst(ndst,k,j,i)=src(v,mk,mj,mi);
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh)
//! \brief Set the defect, including the ghost zone

void Multigrid::RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh) {
  const AthenaArray<Real> &src=def_[nlevel_-1];
  int sngh=std::min(ngh_,ngh);
  int ie=size_.nx1+ngh_+sngh-1, je=size_.nx2+ngh_+sngh-1, ke=size_.nx3+ngh_+sngh-1;
  for (int v=0; v<nvar_; ++v) {
    int ndst=ns+v;
    for (int mk=ngh_-sngh; mk<=ke; ++mk) {
      int k = mk - ngh_ + ngh;
      for (int mj=ngh_-sngh; mj<=je; ++mj) {
        int j = mj - ngh_ + ngh;
#pragma omp simd
        for (int mi=ngh_-sngh; mi<=ie; ++mi) {
          int i = mi - ngh_ + ngh;
          dst(ndst,k,j,i)=src(v,mk,mj,mi)*defscale_;
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ZeroClearData()
//! \brief Clear the data array with zero

void Multigrid::ZeroClearData() {
  u_[current_level_].ZeroClear();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictBlock()
//! \brief Restrict the defect to the source

void Multigrid::RestrictBlock() {
  int ll=nlevel_-current_level_;
  int is, ie, js, je, ks, ke;
  int th = false;
#ifdef OPENMP_PARALLEL
  if (pmy_block_ == nullptr)
    th = true;
#endif

  CalculateDefectBlock();
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  Restrict(src_[current_level_-1], def_[current_level_],
           nvar_, is, ie, js, je, ks, ke, th);

  // Full Approximation Scheme - restrict the variable itself
  if (pmy_driver_->ffas_)
    Restrict(u_[current_level_-1], u_[current_level_],
             nvar_, is, ie, js, je, ks, ke, th);

  current_level_--;

  if (!pmy_driver_->ffas_) ZeroClearData();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrectBlock()
//! \brief Prolongate the potential using tri-linear interpolation

void Multigrid::ProlongateAndCorrectBlock() {
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  int th = false;
#ifdef OPENMP_PARALLEL
  if (pmy_block_ == nullptr)
    th = true;
#endif
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  if (pmy_driver_->ffas_) {
    int size = u_[current_level_].GetSize();
    for (int s=0; s<size; ++s)
      u_[current_level_](s) -= uold_[current_level_](s);
  }

  ProlongateAndCorrect(u_[current_level_+1], u_[current_level_],
                       is, ie, js, je, ks, ke, ngh_, ngh_, ngh_, th);

  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongateBlock()
//! \brief Prolongate the potential for Full Multigrid cycle

void Multigrid::FMGProlongateBlock() {
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  int th = false;
#ifdef OPENMP_PARALLEL
  if (pmy_block_ == nullptr)
    th = true;
#endif
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;

  FMGProlongate(u_[current_level_+1], u_[current_level_],
                is, ie, js, je, ks, ke, ngh_, ngh_, ngh_, th);

  current_level_++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::SmoothBlock(int color)
//! \brief Apply Smoother on the Block

void Multigrid::SmoothBlock(int color) {
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  int th = false;
#ifdef OPENMP_PARALLEL
  if (pmy_block_ == nullptr)
    th = true;
#endif
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;

  Smooth(u_[current_level_], src_[current_level_],  coeff_[current_level_],
         matrix_[current_level_], -ll, is, ie, js, je, ks, ke, color, th);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateDefectBlock()
//! \brief calculate the residual

void Multigrid::CalculateDefectBlock() {
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  int th = false;
#ifdef OPENMP_PARALLEL
  if (pmy_block_ == nullptr)
    th = true;
#endif
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;

  CalculateDefect(def_[current_level_], u_[current_level_], src_[current_level_],
                  coeff_[current_level_], matrix_[current_level_],
                  -ll, is, ie, js, je, ks, ke, th);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateFASRHSBlock()
//! \brief calculate the RHS for the Full Approximation Scheme

void Multigrid::CalculateFASRHSBlock() {
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  int th = false;
#ifdef OPENMP_PARALLEL
  if (pmy_block_ == nullptr)
    th = true;
#endif
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;

  CalculateFASRHS(src_[current_level_], u_[current_level_], coeff_[current_level_],
                  matrix_[current_level_], -ll, is, ie, js, je, ks, ke, th);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateMatrixBlock(Real dt)
//  \brief calculate matrix elements for all the levels

void Multigrid::CalculateMatrixBlock(Real dt) {
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  for (int lev = nlevel_ - 1; lev >= 0; lev--) {
    int ll = nlevel_ - lev - 1;
    ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
    CalculateMatrix(matrix_[lev], coeff_[lev], dt, -ll, is, ie, js, je, ks, ke, false);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetFromRootGrid(bool folddata)
//! \brief Load the data from the root grid or octets

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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
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
#pragma ivdep
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
//! \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(MGNormType nrm, int n) {
  AthenaArray<Real> &def=def_[current_level_];
  int ll=nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
  Real dx=rdx_*static_cast<Real>(1<<ll), dy=rdy_*static_cast<Real>(1<<ll),
       dz=rdz_*static_cast<Real>(1<<ll);

  CalculateDefect(def_[current_level_], u_[current_level_], src_[current_level_],
                  coeff_[current_level_], matrix_[current_level_],
                  -ll, is, ie, js, je, ks, ke, false);

  Real norm=0.0;
  if (nrm == MGNormType::max) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd reduction(max: norm)
        for (int i=is; i<=ie; ++i)
          norm = std::max(norm, std::abs(def(n,k,j,i)));
      }
    }
    return norm;
  } else if (nrm == MGNormType::l1) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd reduction(+: norm)
        for (int i=is; i<=ie; ++i)
          norm += std::abs(def(n,k,j,i));
      }
    }
  } else { // L2 norm
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd reduction(+: norm)
        for (int i=is; i<=ie; ++i)
          norm += SQR(def(n,k,j,i));
      }
    }
  }
  return norm*dx*dy*dz*defscale_;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotal(MGVariable type, int n)
//! \brief calculate the sum of the array (type: 0=src, 1=u)

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
#pragma omp simd reduction(+: s)
      for (int i=is; i<=ie; ++i)
        s+=src(n,k,j,i);
    }
  }
  return s*dx*dy*dz;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverage(MGVariable type, int v, Real ave)
//! \brief subtract the average value (type: 0=src, 1=u)

void Multigrid::SubtractAverage(MGVariable type, int n, Real ave) {
  AthenaArray<Real> &dst = (type == MGVariable::src) ? src_[nlevel_-1] : u_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=0;
  ie=is+size_.nx1+1, je=js+size_.nx2+1, ke=ks+size_.nx3+1;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        dst(n,k,j,i)-=ave;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::StoreOldData()
//! \brief store the old u data in the uold array

void Multigrid::StoreOldData() {
  memcpy(uold_[current_level_].data(), u_[current_level_].data(),
         u_[current_level_].GetSizeInBytes());

  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::GetCoarsestData(MGVariable type, int n)
//  \brief get the value on the coarsest level in the MG block

Real Multigrid::GetCoarsestData(MGVariable type, int n) {
  if (type == MGVariable::src)
    return src_[0](n, ngh_, ngh_, ngh_);
  else if (type == MGVariable::u)
    return u_[0](n, ngh_, ngh_, ngh_);
  else
    return coeff_[0](n, ngh_, ngh_, ngh_);
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetData(MGVariable type, int n, int k, int j, int i, Real v)
//! \brief set a value to a cell on the current level

void Multigrid::SetData(MGVariable type, int n, int k, int j, int i, Real v) {
  if (type == MGVariable::src)
    src_[current_level_](n, ngh_+k, ngh_+j, ngh_+i) = v;
  else if (type == MGVariable::u)
    u_[current_level_](n, ngh_+k, ngh_+j, ngh_+i) = v;
  else
    coeff_[current_level_](n, ngh_+k, ngh_+j, ngh_+i) = v;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
//                      int nvar, int il, int iu, int jl, int ju, int kl, int ku, bool th)
//  \brief Actual implementation of prolongation and correction

void Multigrid::Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                int nvar, int il, int iu, int jl, int ju, int kl, int ku, bool th) {
  for (int v=0; v<nvar; ++v) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl) >= minth_)
    for (int k=kl; k<=ku; ++k) {
      int fk = 2*k - kl;
      for (int j=jl; j<=ju; ++j) {
        int fj = 2*j - jl;
#pragma ivdep
        for (int i=il; i<=iu; ++i) {
          int fi = 2*i - il;
          dst(v, k, j, i)=0.125*(src(v, fk,   fj,   fi)+src(v, fk,   fj,   fi+1)
                                +src(v, fk,   fj+1, fi)+src(v, fk,   fj+1, fi+1)
                                +src(v, fk+1, fj,   fi)+src(v, fk+1, fj,   fi+1)
                                +src(v, fk+1, fj+1, fi)+src(v, fk+1, fj+1, fi+1));
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(AthenaArray<Real> &dst,
//!     const AthenaArray<Real> &src, int il, int iu, int jl, int ju, int kl, int ku,
//!     int fil, int fjl, int fkl, bool th)
//! \brief Actual implementation of prolongation and correction

void Multigrid::ProlongateAndCorrect(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
     int il, int iu, int jl, int ju, int kl, int ku, int fil, int fjl, int fkl, bool th) {
  if (pmy_driver_->fprolongation_ == 1) { // tricubic
    for (int v=0; v<nvar_; ++v) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl)>=minth_/2)
      for (int k=kl; k<=ku; ++k) {
        int fk = 2*(k-kl) + fkl;
        for (int j=jl; j<=ju; ++j) {
          int fj = 2*(j-jl) + fjl;
#pragma ivdep
          for (int i=il; i<=iu; ++i) {
            int fi = 2*(i-il) + fil;
          dst(v,fk  ,fj,  fi  ) += (
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
          dst(v,fk,  fj,  fi+1) += (
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
          dst(v,fk  ,fj+1,fi  ) += (
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
          dst(v,fk,  fj+1,fi+1) += (
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
          dst(v,fk+1,fj,  fi  ) += (
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
          dst(v,fk+1,fj,  fi+1) += (
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
          dst(v,fk+1,fj+1,fi  ) += (
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
          dst(v,fk+1,fj+1,fi+1) += (
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
  } else { // trilinear
    for (int v=0; v<nvar_; ++v) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl)>=minth_/2)
      for (int k=kl; k<=ku; ++k) {
        int fk = 2*(k-kl) + fkl;
        for (int j=jl; j<=ju; ++j) {
          int fj = 2*(j-jl) + fjl;
#pragma ivdep
          for (int i=il; i<=iu; ++i) {
            int fi = 2*(i-il) + fil;
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
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongate(AthenaArray<Real> &dst,
//!          const AthenaArray<Real> &src, int il, int iu, int jl, int ju, int kl, int ku
//!          int fil, int fjl, int fkl, bool th)
//! \brief Actual implementation of FMG prolongation

void Multigrid::FMGProlongate(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                              int il, int iu, int jl, int ju, int kl, int ku,
                              int fil, int fjl, int fkl, bool th) {
  for (int v=0; v<nvar_; ++v) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl)>=minth_/2)
    for (int k=kl; k<=ku; ++k) {
      int fk = 2*(k-kl) + fkl;
      for (int j=jl; j<=ju; ++j) {
        int fj = 2*(j-jl) + fjl;
#pragma ivdep
        for (int i=il; i<=iu; ++i) {
          int fi = 2*(i-il) + fil;
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


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateMultipoleCoefficients(AthenaArray<Real> &mpcoeff)
//! \brief Actual implementation of calculation of multipole expansion coeficients

void Multigrid::CalculateMultipoleCoefficients(AthenaArray<Real> &mpcoeff) {
  AthenaArray<Real> &src = src_[nlevel_-1];
  MGCoordinates &coord = coord_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  // *** Note ***: Currently this calculates coefficients of the zeroth variable only.
  // It is trivial to extend it, but I'm afraid it slows down the code considerably
  // as it requires non-continuous memory access.
  Real vol = (coord.x1f(is+1)-coord.x1f(is)) * (coord.x2f(js+1)-coord.x2f(js))
           * (coord.x3f(ks+1)-coord.x3f(ks));
  Real xorigin = pmy_driver_->mpo_(0);
  Real yorigin = pmy_driver_->mpo_(1);
  Real zorigin = pmy_driver_->mpo_(2);
  if (pmy_driver_->mporder_ == 4) {
    Real m0=0.0, m1=0.0, m2=0.0, m3=0.0, m4=0.0, m5=0.0, m6=0.0, m7=0.0, m8=0.0, m9=0.0,
         m10=0.0, m11=0.0, m12=0.0, m13=0.0, m14=0.0, m15=0.0, m16=0.0, m17=0.0, m18=0.0,
         m19=0.0, m20=0.0, m21=0.0, m22=0.0, m23=0.0, m24=0.0;
    if (pmy_driver_->nodipole_) {
      for (int k = ks; k <= ke; ++k) {
        Real z = coord.x3v(k) - zorigin;
        Real z2 = z*z;
        for (int j = js; j <= je; ++j) {
          Real y = coord.x2v(j) - yorigin;
          Real y2 = y*y, yz = y*z;
#pragma ivdep
          for (int i = is; i <= ie; ++i) {
            Real x = coord.x1v(i) - xorigin;
            Real x2 = x*x, xy = x*y, zx = z*x;
            Real r2 = x2 + y2 + z2;
            Real hx2my2 = 0.5*(x2-y2);
            Real x2mty2 = x2-3.0*y2;
            Real tx2my2 = 3.0*x2-y2;
            Real fz2mr2 = 5.0*z2-r2;
            Real sz2mr2 = 7.0*z2-r2;
            Real sz2mtr2 = 7.0*z2-3.0*r2;
            Real s = src(k,j,i) * vol;
            // Y00
            m0  += s;
            // r^2*(Y2-2, Y2-1, Y20, Y21, Y22)
            m4  += s*xy;
            m5  += s*yz;
            m6  += s*(3.0*z2-r2);
            m7  += s*zx;
            m8  += s*hx2my2;
            // r^3*(Y3-3, Y3-2, Y3-1, Y30, Y31, Y32, Y33)
            m9  += s*y*tx2my2;
            m10 += s*xy*z;
            m11 += s*y*fz2mr2;
            m12 += s*z*(z2-3.0*r2);
            m13 += s*x*fz2mr2;
            m14 += s*z*hx2my2;
            m15 += s*x*x2mty2;
            // r^3*(Y3-3, Y3-2, Y3-1, Y30, Y31, Y32, Y33)
            m16 += s*xy*hx2my2;
            m17 += s*yz*tx2my2;
            m18 += s*xy*sz2mr2;
            m19 += s*yz*sz2mtr2;
            m20 += s*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2);
            m21 += s*zx*sz2mtr2;
            m22 += s*hx2my2*sz2mr2;
            m23 += s*zx*x2mty2;
            m24 += s*0.125*(x2*x2mty2-y2*tx2my2);
          }
        }
      }
    } else {
      for (int k = ks; k <= ke; ++k) {
        Real z = coord.x3v(k) - zorigin;
        Real z2 = z*z;
        for (int j = js; j <= je; ++j) {
          Real y = coord.x2v(j) - yorigin;
          Real y2 = y*y, yz = y*z;
#pragma ivdep
          for (int i = is; i <= ie; ++i) {
            Real x = coord.x1v(i) - xorigin;
            Real x2 = x*x, xy = x*y, zx = z*x;
            Real r2 = x2 + y2 + z2;
            Real hx2my2 = 0.5*(x2-y2);
            Real x2mty2 = x2-3.0*y2;
            Real tx2my2 = 3.0*x2-y2;
            Real fz2mr2 = 5.0*z2-r2;
            Real sz2mr2 = 7.0*z2-r2;
            Real sz2mtr2 = 7.0*z2-3.0*r2;
            Real s = src(k,j,i) * vol;
            // Y00
            m0  += s;
            // r*(Y1-1, Y10, Y11)
            m1  += s*y;
            m2  += s*z;
            m3  += s*x;
            // r^2*(Y2-2, Y2-1, Y20, Y21, Y22)
            m4  += s*xy;
            m5  += s*yz;
            m6  += s*(3.0*z2-r2);
            m7  += s*zx;
            m8  += s*hx2my2;
            // r^3*(Y3-3, Y3-2, Y3-1, Y30, Y31, Y32, Y33)
            m9  += s*y*tx2my2;
            m10 += s*xy*z;
            m11 += s*y*fz2mr2;
            m12 += s*z*(z2-3.0*r2);
            m13 += s*x*fz2mr2;
            m14 += s*z*hx2my2;
            m15 += s*x*x2mty2;
            // r^4*(Y4-4, Y4-3, Y4-2, Y4-1, Y40, Y41, Y42, Y43, Y44)
            m16 += s*xy*hx2my2;
            m17 += s*yz*tx2my2;
            m18 += s*xy*sz2mr2;
            m19 += s*yz*sz2mtr2;
            m20 += s*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2);
            m21 += s*zx*sz2mtr2;
            m22 += s*hx2my2*sz2mr2;
            m23 += s*zx*x2mty2;
            m24 += s*0.125*(x2*x2mty2-y2*tx2my2);
          }
        }
      }
    }
    mpcoeff(0)  += m0;
    mpcoeff(1)  += m1;
    mpcoeff(2)  += m2;
    mpcoeff(3)  += m3;
    mpcoeff(4)  += m4;
    mpcoeff(5)  += m5;
    mpcoeff(6)  += m6;
    mpcoeff(7)  += m7;
    mpcoeff(8)  += m8;
    mpcoeff(9)  += m9;
    mpcoeff(10) += m10;
    mpcoeff(11) += m11;
    mpcoeff(12) += m12;
    mpcoeff(13) += m13;
    mpcoeff(14) += m14;
    mpcoeff(15) += m15;
    mpcoeff(16) += m16;
    mpcoeff(17) += m17;
    mpcoeff(18) += m18;
    mpcoeff(19) += m19;
    mpcoeff(20) += m20;
    mpcoeff(21) += m21;
    mpcoeff(22) += m22;
    mpcoeff(23) += m23;
    mpcoeff(24) += m24;
  } else if (pmy_driver_->mporder_ == 2) {
    Real m0=0.0, m1=0.0, m2=0.0, m3=0.0, m4=0.0, m5=0.0, m6=0.0, m7=0.0, m8=0.0;
    if (pmy_driver_->nodipole_) {
      for (int k = ks; k <= ke; ++k) {
        Real z = coord.x3v(k) - zorigin;
        Real z2 = z*z;
        for (int j = js; j <= je; ++j) {
          Real y = coord.x2v(j) - yorigin;
          Real y2 = y*y, yz = y*z;
#pragma ivdep
          for (int i = is; i <= ie; ++i) {
            Real x = coord.x1v(i) - xorigin;
            Real x2 = x*x, xy = x*y, zx = z*x;
            Real r2 = x2 + y2 + z2;
            Real s = src(k,j,i) * vol;
            // Y00
            m0 += s;
            // r^2*(Y2-2, Y2-1, Y20, Y21, Y22)
            m4 += s*xy;
            m5 += s*yz;
            m6 += s*(3.0*z2-r2);
            m7 += s*zx;
            m8 += s*0.5*(x2-y2);
          }
        }
      }
    } else {
      for (int k = ks; k <= ke; ++k) {
        Real z = coord.x3v(k) - zorigin;
        Real z2 = z*z;
        for (int j = js; j <= je; ++j) {
          Real y = coord.x2v(j) - yorigin;
          Real y2 = y*y, yz = y*z;
#pragma ivdep
          for (int i = is; i <= ie; ++i) {
            Real x = coord.x1v(i) - xorigin;
            Real x2 = x*x, xy = x*y, zx = z*x;
            Real r2 = x2 + y2 + z2;
            Real s = src(k,j,i) * vol;
            // Y00
            m0 += s;
            // r*(Y1-1, Y10, Y11)
            m1 += s*y;
            m2 += s*z;
            m3 += s*x;
            // r^2*(Y2-2, Y2-1, Y20, Y21, Y22)
            m4 += s*xy;
            m5 += s*yz;
            m6 += s*(3.0*z2-r2);
            m7 += s*zx;
            m8 += s*0.5*(x2-y2);
          }
        }
      }
    }
    mpcoeff(0) += m0;
    mpcoeff(1) += m1;
    mpcoeff(2) += m2;
    mpcoeff(3) += m3;
    mpcoeff(4) += m4;
    mpcoeff(5) += m5;
    mpcoeff(6) += m6;
    mpcoeff(7) += m7;
    mpcoeff(8) += m8;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::CalculateCenterOfMass(AthenaArray<Real> &mpcoeff)
//! \brief Calculate the position of the center of mass from the dipole moment

void Multigrid::CalculateCenterOfMass(AthenaArray<Real> &mpcoeff) {
  AthenaArray<Real> &src = src_[nlevel_-1];
  MGCoordinates &coord = coord_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  Real vol = (coord.x1f(is+1)-coord.x1f(is)) * (coord.x2f(js+1)-coord.x2f(js))
           * (coord.x3f(ks+1)-coord.x3f(ks));
  Real m0 = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j);
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i);
        Real s = src(k,j,i) * vol;
        m0 += s;
        m1 += s*y;
        m2 += s*z;
        m3 += s*x;
      }
    }
  }
  mpcoeff(0) += m0;
  mpcoeff(1) += m1;
  mpcoeff(2) += m2;
  mpcoeff(3) += m3;
}

//----------------------------------------------------------------------------------------
//! \fn void MGCoordinates::AllocateMGCoordinates(int nx, int ny, int nz)
//! \brief Allocate coordinate arrays for multigrid

void MGCoordinates::AllocateMGCoordinates(int nx, int ny, int nz) {
  x1f.NewAthenaArray(nx+1);
  x2f.NewAthenaArray(ny+1);
  x3f.NewAthenaArray(nz+1);
  x1v.NewAthenaArray(nx);
  x2v.NewAthenaArray(ny);
  x3v.NewAthenaArray(nz);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCoordinates::CalculateMGCoordinates(const RegionSize &size,
//!                                                int ll, int ngh)
//! \brief Calculate coordinates for Multigrid
//!        Currently uniform Cartesian only

void MGCoordinates::CalculateMGCoordinates(const RegionSize &size, int ll, int ngh) {
  int ncx = (size.nx1>>ll), ncy = (size.nx2>>ll), ncz = (size.nx3>>ll);
  Real dx = (size.x1max-size.x1min)/ncx;
  Real dy = (size.x2max-size.x2min)/ncy;
  Real dz = (size.x3max-size.x3min)/ncz;
  for (int i = 0; i <= ncx+2*ngh; ++i)
    x1f(i) = size.x1min + (i-ngh)*dx;
  x1f(ngh) = size.x1min;
  x1f(ncx+ngh) = size.x1max;
  for (int i = 0; i < ncx+2*ngh; ++i)
    x1v(i) = 0.5*(x1f(i)+x1f(i+1));

  for (int j = 0; j <= ncy+2*ngh; ++j)
    x2f(j) = size.x2min + (j-ngh)*dy;
  x2f(ngh) = size.x2min;
  x2f(ncy+ngh) = size.x2max;
  for (int j = 0; j < ncy+2*ngh; ++j)
    x2v(j) = 0.5*(x2f(j)+x2f(j+1));

  for (int k = 0; k <= ncz+2*ngh; ++k)
    x3f(k) = size.x3min + (k-ngh)*dz;
  x3f(ngh) = size.x3min;
  x3f(ncz+ngh) = size.x3max;
  for (int k = 0; k < ncz+2*ngh; ++k)
    x3v(k) = 0.5*(x3f(k)+x3f(k+1));
}
