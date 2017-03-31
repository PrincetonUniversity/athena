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
#include "../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns)
{
  AthenaArray<Real> &dst=u_[nlev-1];
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int nsrc=ns+n;
#pragma ivdep
    for(int k=pmb->ks, mk=ngh; k<=pmb->ke; k++, mk++) {
#pragma ivdep
      for(int j=pmb->js, mj=ngh; j<=pmb->je; j++, mj++) {
#pragma ivdep
        for(int i=pmb->is, mi=ngh; i<=pmb->ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, Real fac)
//  \brief Fill the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, Real fac)
{
  AthenaArray<Real> &dst=src_[nlev-1];
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int nsrc=ns+n;
#pragma ivdep
    for(int k=pmb->ks, mk=ngh; k<=pmb->ke; k++, mk++) {
#pragma ivdep
      for(int j=pmb->js, mj=ngh; j<=pmb->je; j++, mj++) {
#pragma ivdep
        for(int i=pmb->is, mi=ngh; i<=pmb->ie; i++, mi++)
          dst(n,mk,mj,mi)=src(nsrc,k,j,i)*fac;
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns)
//  \brief Set the result
void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns)
{
  const AthenaArray<Real> &src=src_[nlev-1];
#pragma ivdep
  for(int n=0; n<nvar; n++) {
    int ndst=ns+n;
#pragma ivdep
    for(int k=pmb->ks, mk=ngh; k<=pmb->ke; k++, mk++) {
#pragma ivdep
      for(int j=pmb->js, mj=ngh; j<=pmb->je; j++, mj++) {
#pragma ivdep
        for(int i=pmb->is, mi=ngh; i<=pmb->ie; i++, mi++)
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
  int ns=ngh, ne=ngh+(1<<clev)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ns, fk=ns; k<=ne; k++, fk+=2) {
#pragma ivdep
      for(int j=ns, fj=ns; j<=ne; j++, fj+=2) {
#pragma ivdep
        for(int i=ns, fi=ns; i<=ne; i++, fi+=2)
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
  int ns=ngh, ne=ngh+(1<<clev)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ns, fk=ns; k<=ne; k++, fk+=2) {
#pragma ivdep
      for(int j=ns, fj=ns; j<=ne; j++, fj+=2) {
#pragma ivdep
        for(int i=ns, fi=ns; i<=ne; i++, fi+=2) {
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
  int ns=ngh, ne=ngh+(1<<clev)-1;
#pragma ivdep
  for(int n=0; n<nvar; n++) {
#pragma ivdep
    for(int k=ns, fk=ns; k<=ne; k++, fk+=2) {
#pragma ivdep
      for(int j=ns, fj=ns; j<=ne; j++, fj+=2) {
#pragma ivdep
        for(int i=ns, fi=ns; i<=ne; i++, fi+=2) {
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
//! \fn Real Multigrid::CalculateDefectNorm(int nrm)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(int nrm)
{
  AthenaArray<Real> &def=def_[nlev-1];
  int ns=ngh, ne=ngh+(1<<nlev)-1;
  Real norm=0.0;
  if(nrm==0) { // special case: max norm
#pragma ivdep
    for(int n=0; n<nvar; n++) {
#pragma ivdep
      for(int k=ns; k<=ne; k++) {
#pragma ivdep
        for(int j=ns; j<=ne; j++) {
#pragma ivdep
          for(int i=ns; i<=ne; i++)
            norm=std::max(norm,std::fabs(def(k,j,i)));
        }
      }
    }
  }
  else if (nrm==1) {
#pragma ivdep
    for(int n=0; n<nvar; n++) {
#pragma ivdep
      for(int k=ns; k<=ne; k++) {
#pragma ivdep
        for(int j=ns; j<=ne; j++) {
#pragma ivdep
          for(int i=ns; i<=ne; i++)
            norm+=std::fabs(def(k,j,i));
        }
      }
    }
  }
  else { // nrm>1 -> nrm=2
#pragma ivdep
    for(int n=0; n<nvar; n++) {
#pragma ivdep
      for(int k=ns; k<=ne; k++) {
#pragma ivdep
        for(int j=ns; j<=ne; j++) {
#pragma ivdep
            for(int i=ns; i<=ne; i++)
            norm+=SQR(def(k,j,i));
        }
      }
    }
    norm=std::sqrt(norm);
  }
  return norm;
}


//----------------------------------------------------------------------------------------
//! \fn  void Multigrid::Multigrid(MeshBlock *pmb, int invar)
//  \brief Multigrid constructor

void Multigrid::Multigrid(MeshBlock *pmb, int invar)
{
  pmy_block=pmb;
  pmgd=pmb->pmy_mesh->pmgd;
  nvar=invar;

  if(pmb->block_size.nx1!=pmb->block_size.nx2
  || pmb->block_size.nx1!=pmb->block_size.nx3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Multigrid::Initialize" << std::endl
        << "The Multigrid solver requires cubic MeshBlocks." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pmb->block_size.nx2==1 || pmb->block_size.nx3==1 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Multigrid::Initialize" << std::endl
        << "Currently the Multigrid solver works only in 3D." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pmb->pmy_mesh->use_meshgen_fn_[X1DIR]==true
  || pmb->pmy_mesh->use_meshgen_fn_[X2DIR]==true
  || pmb->pmy_mesh->use_meshgen_fn_[X3DIR]==true) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Multigrid::Initialize" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  dx_=pmb->pcoord->dx1f(0);
  if(dx_!=pmb->pcoord->dx2f(0) || dx_!=pmb->pcoord->dx3f(0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Multigrid::Initialize" << std::endl
        << "The cell size must be cubic." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  // count multigrid levels
  nlev=0;
  for(int l=0; l<20; l++) {
    if((1<<l) == pmb->block_size.nx1) {
      nlev=l+1;
      break;
    }
  }
  if(nlev==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Multigrid::Initialize" << std::endl
        << "The MeshBlock size must be power of two." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  // allocate arrays
  u_ = new AthenaArray<Real>[nlev];
  src_ = new AthenaArray<Real>[nlev];
  def_ = new AthenaArray<Real>[nlev];
  for(int l=0; l<nlev; l++) {
    int nc=(1<<l)+2*ngh;
    u_[l].NewAthenaArray(nvar,nc,nc,nc);
    src_[l].NewAthenaArray(nvar,nc,nc,nc);
    def_[l].NewAthenaArray(nvar,nc,nc,nc);
  }
}

//----------------------------------------------------------------------------------------
//! \fn  virtual void Multigrid::~Multigrid
//  \brief Multigrid destroctor

virtual void Multigrid::~Multigrid()
{
  for(int l=0; l<nlev; l++) {
    u_[l].DeleteAthenaArray();
    src_[l].DeleteAthenaArray();
    def_[l].DeleteAthenaArray();
  }
  delete [] u_;
  delete [] src_;
  delete [] def_;
}

