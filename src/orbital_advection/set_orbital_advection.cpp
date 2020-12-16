//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file set_orbital_advection.cpp
//! \brief functions to put variables into buffers for orbital communication

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs
#include <cstring>    // memcpy
#include <iostream>   // cout, endl
#include <sstream>    //
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "orbital_advection.hpp"


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetOrbitalAdvectionCC(const AthenaArray<Real> &u,
//!                                                  const AthenaArray<Real> &s)
//! \brief put cell-centered variables of this meshblock into orbital buffer

void OrbitalAdvection::SetOrbitalAdvectionCC(const AthenaArray<Real> &u,
                                             const AthenaArray<Real> &s) {
  int is = pmb_->is, ie = pmb_->ie;
  int js = pmb_->js, je = pmb_->je;
  int ks = pmb_->ks, ke = pmb_->ke;

  // set variables of this meshblock into orbital blocks
  if (orbital_direction ==1) { // cartesian or cylindrical
    // hydro variables
    for (int k = ks; k <= ke; k++) {
      for (int i = is; i <= ie; i++) {
        int offset = ofc(k,i);
        const int shift = (offset>0)?onx:0;
        int jl = std::max(js, js-xgh-offset);
        int ju = std::min(je, je+1+xgh-offset);
        for (int nph = 0; nph < NHYDRO; nph++) {
#pragma omp simd
          for (int j = jl; j <= ju; j++) {
            orbital_cons(nph,k,i,j+shift) = u(nph,k,j,i);
          }
        }
      }
    }
    // passive scalars
    if (NSCALARS>0) {
      for (int k = ks; k <= ke; k++) {
        for (int i = is; i <= ie; i++) {
          int offset = ofc(k,i);
          const int shift = (offset>0)?onx:0;
          int jl = std::max(js, js-xgh-offset);
          int ju = std::min(je, je+1+xgh-offset);
          for (int nsc = 0; nsc < NSCALARS; nsc++) {
#pragma omp simd
            for (int j = jl; j <= ju; j++) {
              orbital_scalar(nsc,k,i,j+shift)= s(nsc,k,j,i);
            }
          }
        }
      }
    }
  } else if (orbital_direction == 2) { // spherical_polar
    // hydro
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        int offset = ofc(j,i);
        const int shift = (offset>0)?onx:0;
        int kl = std::max(ks, ks-xgh-offset);
        int ku = std::min(ke, ke+1+xgh-offset);
        for (int nph = 0; nph < NHYDRO; nph++) {
#pragma omp simd
          for (int k = kl; k <= ku; k++) {
            orbital_cons(nph,j,i,k+shift) = u(nph,k,j,i);
          }
        }
      }
    }
    // passive scalars
    if (NSCALARS>0) {
      for (int j = js; j <= je; j++) {
        for (int i = is; i <= ie; i++) {
          int offset = ofc(j,i);
          const int shift = (offset>0)?onx:0;
          int kl = std::max(ks, ks-xgh-offset);
          int ku = std::min(ke, ke+1+xgh-offset);
          for (int nsc = 0; nsc < NSCALARS; nsc++) {
#pragma omp simd
            for (int k = kl; k <= ku; k++) {
              orbital_scalar(nsc,j,i,k+shift)= s(nsc,k,j,i);
            }
          }
        }
      }
    }
  }
  if (orbital_refinement) {
    // restriction for orbital communications with refinement
    pmb_->pmr->RestrictCellCenteredValues(u, u_coarse_send, 0, NHYDRO-1,
                                          pmb_->cis, pmb_->cie, pmb_->cjs,
                                          pmb_->cje, pmb_->cks, pmb_->cke);
    if (NSCALARS>0)
      pmb_->pmr->RestrictCellCenteredValues(s, s_coarse_send, 0, NSCALARS-1,
                                            pmb_->cis, pmb_->cie, pmb_->cjs,
                                            pmb_->cje, pmb_->cks, pmb_->cke);
  }
return;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetOrbitalAdvectionFC(const FaceField &b)
//! \brief put face-centered variables of this meshblock into orbital buffer

void OrbitalAdvection::SetOrbitalAdvectionFC(const FaceField &b) {
  int is = pmb_->is, ie = pmb_->ie;
  int js = pmb_->js, je = pmb_->je;
  int ks = pmb_->ks, ke = pmb_->ke;

  //set variables of this meshblock into orbital blocks
  if (orbital_direction ==1) { // cartesian or cylindrical
    // put b1 into orbital_b1
    for (int k = ks; k <= ke  ; k++) {
      for (int i = is; i <= ie+1; i++) {
        const int shift = (off[0](k,i)>0)?onx:0;
#pragma omp simd
        for (int j = js; j <= je; j++) {
          orbital_b1(k,i,j+shift) = b.x1f(k,j,i);
        }
      }
    }
    // put -b3 into orbital_b2
    for (int k = ks; k <= ke+1; k++) {
      for (int i = is; i <= ie  ; i++) {
        const int shift = (off[1](k,i)>0)?onx:0;
#pragma omp simd
        for (int j = js; j <= je; j++) {
          orbital_b2(k,i,j+shift) = -b.x3f(k,j,i);
        }
      }
    }
  } else if (orbital_direction == 2) { // spherical_polar
    // put -b1 into orbital_b1
    for (int j = js; j <= je  ; j++) {
      for (int i = is; i <= ie+1; i++) {
        const int shift = (off[0](j,i)>0)?onx:0;
#pragma omp simd
        for (int k = ks; k <= ke; k++) {
          orbital_b1(j,i,k+shift) = -b.x1f(k,j,i);
        }
      }
    }
    // put b2 into orbital_b2
    for (int j = js; j <= je+1; j++) {
      for (int i = is; i <= ie  ; i++) {
        const int shift = (off[1](j,i)>0)?onx:0;
#pragma omp simd
        for (int k = ks; k <= ke; k++) {
          orbital_b2(j,i,k+shift) = b.x2f(k,j,i);
        }
      }
    }
  }
  if (orbital_refinement) {
    // restriction for orbital communications with refinement
    pmb_->pmr->RestrictFieldX1(b.x1f, b1_coarse_send, pmb_->cis,
                               pmb_->cie+1, pmb_->cjs, pmb_->cje,
                               pmb_->cks, pmb_->cke);
    if (orbital_direction ==1) { // cartesian or cylindrical
      pmb_->pmr->RestrictFieldX3(b.x3f, b2_coarse_send, pmb_->cis,
                                 pmb_->cie, pmb_->cjs, pmb_->cje,
                                 pmb_->cks, pmb_->cke+1);
    } else if (orbital_direction == 2) { // spherical_polar
      pmb_->pmr->RestrictFieldX2(b.x2f, b2_coarse_send, pmb_->cis,
                                 pmb_->cie, pmb_->cjs, pmb_->cje+1,
                                 pmb_->cks, pmb_->cke);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetOrbitalEdgeCC(const Real dt,
//!                                             int *ssize[2], int *rsize[2])
//! \brief set orbital edge for cell-centered variables

void OrbitalAdvection::SetOrbitalEdgeCC(const Real dt, int *ssize[2], int *rsize[2]) {
  int is = pmb_->is, ie = pmb_->ie;
  int js = pmb_->js, je = pmb_->je;
  int ks = pmb_->ks, ke = pmb_->ke;

  if (orbital_uniform_mesh) { // uniform mesh
    // for communication with meshblock at same level
    int p1=0; int p2=0;
    if (orbital_direction == 1) { // cartesian or cylindrical
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie; i++) {
          Real olen = dt*vKc(k,i)/pco_->h2v(i);
          int offset = static_cast<int>(olen/dx);
          if (olen > 0.0) offset++;
          orc(k,i)  = std::fmod(olen, dx)/dx;
          ofc(k,i)  = offset;
          p1 += std::max(offset+xgh, 0);
          p2 += std::max(1+xgh-offset, 0);
        }
      }
    } else if(orbital_direction == 2) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++) {
          Real olen = dt*vKc(j,i)/(pco_->h2v(i)*pco_->h32v(j));
          int offset = static_cast<int>(olen/dx);
          if (olen>0.0) offset++;
          orc(j,i)  = std::fmod(olen, dx)/dx;
          ofc(j,i)  = offset;
          p1 += std::max(offset+xgh, 0);
          p2 += std::max(1+xgh-offset, 0);
        }
      }
    }
    ssize[0][0] = p1;
    rsize[0][0] = p1;
    ssize[1][0] = p2;
    rsize[1][0] = p2;
    if (orbital_refinement) { // with orbital refinement
      // communication with coarser level
      p1 = 0; p2 = 0;
      Real cdx = 2.0*dx;
      max_ofc_coarse = -onx;
      min_ofc_coarse = onx;
      Coordinates *cpco = pmb_->pmr->pcoarsec;
      if (orbital_direction == 1) {
        for(int k=pmb_->cks; k<=pmb_->cke; k++) {
          for(int i=pmb_->cis; i<=pmb_->cie; i++) {
            Real olen = dt*vKc_coarse(k,i)/cpco->h2v(i);
            int offset = static_cast<int>(olen/cdx);
            if (olen > 0.0) offset++;
            ofc_coarse(k,i)  = offset;
            max_ofc_coarse = std::max(offset+1, max_ofc_coarse);
            min_ofc_coarse = std::min(offset-1, min_ofc_coarse);
            p1 += std::max(offset+xgh, 0);
            p2 += std::max(1+xgh-offset, 0);
          }
        }
      } else if(orbital_direction == 2) { // spherical_polar
        for(int j=pmb_->cjs; j<=pmb_->cje; j++) {
          for(int i=pmb_->cis; i<=pmb_->cie; i++) {
            Real olen = dt*vKc_coarse(j,i)/(cpco->h2v(i)*cpco->h32v(j));
            int offset = static_cast<int>(olen/cdx);
            if (olen > 0.0) offset++;
            ofc_coarse(j,i)  = offset;
            max_ofc_coarse = std::max(offset+1, max_ofc_coarse);
            min_ofc_coarse = std::min(offset-1, min_ofc_coarse);
            p1 += std::max(offset+xgh, 0);
            p2 += std::max(1+xgh-offset, 0);
          }
        }
      }
      ssize[0][1] = p1;
      ssize[1][1] = p2;
      rsize[0][1] = std::max(max_ofc_coarse+xgh+2,0);
      rsize[1][1] = std::max(3+xgh-min_ofc_coarse,0);

      // communication with finner level
      int hnx1 = pmb_->block_size.nx1/2;
      int hnx2 = pmb_->block_size.nx2/2;
      int hnx3 = pmb_->block_size.nx3/2;
      if (orbital_direction == 1) { // cartesian or cylindrical
        if (hnx3!=0) {
          for (int n3=0; n3<2; n3++) {
            for (int n1=0; n1<2; n1++) {
              p1 = 0; p2 = 0;
              int max_ofc = -hnx2; int min_ofc = hnx2;
              for(int k=ks+n3*hnx3; k<=ke-(1-n3)*hnx3; k++) {
                for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1; i++) {
                  int offset = ofc(k,i);
                  max_ofc = std::max(offset+1, max_ofc);
                  min_ofc = std::min(offset-1, min_ofc);
                  p1 += std::max(offset+xgh, 0);
                  p2 += std::max(1+xgh-offset, 0);
                }
              }
              ssize[0][2+(n1+2*n3)] = std::max(max_ofc+xgh+2,0);
              ssize[1][2+(n1+2*n3)] = std::max(3+xgh-min_ofc,0);
              rsize[0][2+(n1+2*n3)] = p1;
              rsize[1][2+(n1+2*n3)] = p2;
            }
          }
        } else {
          for (int n1=0; n1<2; n1++) {
            p1 = 0; p2 = 0;
            int max_ofc = -hnx2; int min_ofc = hnx2;
            int k = ks;
            for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1; i++) {
              int offset = ofc(k,i);
              max_ofc = std::max(offset+1, max_ofc);
              min_ofc = std::min(offset-1, min_ofc);
              p1 += std::max(offset+xgh, 0);
              p2 += std::max(1+xgh-offset, 0);
            }
            ssize[0][2+n1] = std::max(max_ofc+xgh+2,0);
            ssize[1][2+n1] = std::max(3+xgh-min_ofc,0);
            rsize[0][2+n1] = p1;
            rsize[1][2+n1] = p2;
          }
        }
      } else if(orbital_direction == 2) { // spherical_polar
        for (int n2=0; n2<2; n2++) {
          for (int n1=0; n1<2; n1++) {
            p1 = 0; p2 = 0;
            int max_ofc = -hnx3; int min_ofc = hnx3;
            for(int j=js+n2*hnx2; j<=je-(1-n2)*hnx2; j++) {
              for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1; i++) {
                int offset = ofc(j,i);
                max_ofc = std::max(ofc(j,i)+1, max_ofc);
                min_ofc = std::min(ofc(j,i)-1, min_ofc);
                p1 += std::max(offset+xgh, 0);
                p2 += std::max(1+xgh-offset, 0);
              }
            }
            ssize[0][2+(n1+2*n2)] = std::max(max_ofc+xgh+2,0);
            ssize[1][2+(n1+2*n2)] = std::max(3+xgh-min_ofc,0);
            rsize[0][2+(n1+2*n2)] = p1;
            rsize[1][2+(n1+2*n2)] = p2;
          }
        }
      }
    }
  }
  // else { // non-uniform mesh
  // }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::SetOrbitalEdgeFC(const Real dt,
//!                                             int *ssize[2], int *rsize[2])
//! \brief set orbital edge for face-centered variables

void OrbitalAdvection::SetOrbitalEdgeFC(const Real dt, int *ssize[2], int *rsize[2]) {
  int is = pmb_->is, ie = pmb_->ie;
  int js = pmb_->js, je = pmb_->je;
  int ks = pmb_->ks, ke = pmb_->ke;

  if (orbital_uniform_mesh) { // uniform mesh
    // for communication with meshblock at same level
    int p1=0; int p2=0;
    if (orbital_direction == 1) { // cartesian or cylindrical
      for(int k=ks; k<=ke  ; k++) {
        for(int i=is; i<=ie+1; i++) {
          Real olen = dt*vKf[0](k,i)/pco_->h2f(i);
          int offset = static_cast<int>(olen/dx);
          if (olen>0.0) offset++;
          orf[0](k,i) = std::fmod(olen,dx)/dx;
          off[0](k,i) = offset;
          p1+=std::max(offset+xgh, 0);
          p2+=std::max(1+xgh-offset, 0);
        }
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie  ; i++) {
          Real olen = dt*vKf[1](k,i)/pco_->h2v(i);
          int offset = static_cast<int>(olen/dx);
          if (olen>0.0) offset++;
          orf[1](k,i) = std::fmod(olen,dx)/dx;
          off[1](k,i) = offset;
          p1+=std::max(offset+xgh, 0);
          p2+=std::max(1+xgh-offset, 0);
        }
      }
    } else if (orbital_direction == 2) { // spherical_polar
      for(int j=js; j<=je  ; j++) {
        for(int i=is; i<=ie+1; i++) {
          Real olen = dt*vKf[0](j,i)/(pco_->h2f(i)*pco_->h32v(j));
          int offset = static_cast<int>(olen/dx);
          if (olen>0.0) offset++;
          orf[0](j,i) = std::fmod(olen,dx)/dx;
          off[0](j,i) = offset;
          p1 += std::max(offset+xgh, 0);
          p2 += std::max(1+xgh-offset, 0);
        }
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie  ; i++) {
          Real h32f = pco_->h32f(j);
          Real olen = (h32f > 0.0) ? dt*vKf[1](j,i)/(pco_->h2v(i)*h32f):0.0;
          int offset = static_cast<int>(olen/dx);
          if (olen>0.0) offset++;
          orf[1](j,i) = std::fmod(olen,dx)/dx;
          off[1](j,i) = offset;
          p1 += std::max(offset+xgh, 0);
          p2 += std::max(1+xgh-offset, 0);
        }
      }
    }
    ssize[0][0] = p1;
    ssize[1][0] = p2;
    rsize[0][0] = p1;
    rsize[1][0] = p2;
    if (orbital_refinement) { // orbital refinement
      // for communication with meshblock at coarse level
      p1 = 0; p2 = 0;
      Real cdx = 2.0*dx;
      max_off_coarse = -onx;
      min_off_coarse = onx;
      Coordinates *cpco = pmb_->pmr->pcoarsec;
      if (orbital_direction == 1) { // cartesian or cylindrical
        for(int k=pmb_->cks; k<=pmb_->cke  ; k++) {
          for(int i=pmb_->cis; i<=pmb_->cie+1; i++) {
            Real olen = dt*vKf_coarse[0](k,i)/cpco->h2f(i);
            int offset = static_cast<int>(olen/cdx);
            if (olen > 0.0) offset++;
            off_coarse[0](k,i)  = offset;
            max_off_coarse = std::max(offset+1, max_off_coarse);
            min_off_coarse = std::min(offset-1, min_off_coarse);
            p1 += std::max(offset+xgh, 0);
            p2 += std::max(1+xgh-offset, 0);
          }
        }
        for(int k=pmb_->cks; k<=pmb_->cke+1; k++) {
          for(int i=pmb_->cis; i<=pmb_->cie  ; i++) {
            Real olen = dt*vKf_coarse[1](k,i)/cpco->h2v(i);
            int offset = static_cast<int>(olen/cdx);
            if (olen > 0.0) offset++;
            off_coarse[1](k,i)  = offset;
            max_off_coarse = std::max(offset+1, max_off_coarse);
            min_off_coarse = std::min(offset-1, min_off_coarse);
            p1 += std::max(offset+xgh, 0);
            p2 += std::max(1+xgh-offset, 0);
          }
        }
      } else if(orbital_direction == 2) { // spherical_polar
        for(int j=pmb_->cjs; j<=pmb_->cje  ; j++) {
          for(int i=pmb_->cis; i<=pmb_->cie+1; i++) {
            Real olen = dt*vKf_coarse[0](j,i)/(cpco->h2f(i)*cpco->h32v(j));
            int offset = static_cast<int>(olen/cdx);
            if (olen > 0.0) offset++;
            off_coarse[0](j,i)  = offset;
            max_off_coarse = std::max(offset+1, max_off_coarse);
            min_off_coarse = std::min(offset-1, min_off_coarse);
            p1 += std::max(offset+xgh, 0);
            p2 += std::max(1+xgh-offset, 0);
          }
        }
        for(int j=pmb_->cjs; j<=pmb_->cje+1; j++) {
          for(int i=pmb_->cis; i<=pmb_->cie  ; i++) {
            Real olen = dt*vKf_coarse[1](j,i)/(cpco->h2v(i)*cpco->h32f(j));
            int offset = static_cast<int>(olen/cdx);
            if (olen > 0.0) offset++;
            off_coarse[1](j,i)  = offset;
            max_off_coarse = std::max(offset+1, max_off_coarse);
            min_off_coarse = std::min(offset-1, min_off_coarse);
            p1 += std::max(offset+xgh, 0);
            p2 += std::max(1+xgh-offset, 0);
          }
        }
      }
      ssize[0][1] = p1;
      ssize[1][1] = p2;
      rsize[0][1] = std::max(max_off_coarse+xgh,0);
      rsize[1][1] = std::max(1+xgh-min_off_coarse,0);

      // for communication with meshblock at finner level
      int hnx1 = pmb_->block_size.nx1/2;
      int hnx2 = pmb_->block_size.nx2/2;
      int hnx3 = pmb_->block_size.nx3/2;
      if (orbital_direction == 1) { // cartesian or cylindrical
        if (hnx3!=0) { // 3D
          for (int n3=0; n3<2; n3++) {
            for (int n1=0; n1<2; n1++) {
              p1 = 0; p2 = 0;
              int max_off = -hnx2; int min_off = hnx2;
              for(int k=ks+n3*hnx3; k<=ke-(1-n3)*hnx3  ; k++) {
                for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1+1; i++) {
                  int offset = off[0](k,i);
                  max_off = std::max(offset+1, max_off);
                  min_off = std::min(offset-1, min_off);
                  p1 += std::max(offset+xgh, 0);
                  p2 += std::max(1+xgh-offset, 0);
                }
              }

              for(int k=ks+n3*hnx3; k<=ke-(1-n3)*hnx3+1; k++) {
                for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1  ; i++) {
                  int offset = off[1](k,i);
                  max_off = std::max(offset+1, max_off);
                  min_off = std::min(offset-1, min_off);
                  p1 += std::max(offset+xgh, 0);
                  p2 += std::max(1+xgh-offset, 0);
                }
              }
              ssize[0][2+(n1+2*n3)] = std::max(max_off+xgh,0);
              ssize[1][2+(n1+2*n3)] = std::max(1+xgh-min_off,0);
              rsize[0][2+(n1+2*n3)] = p1;
              rsize[1][2+(n1+2*n3)] = p2;
            }
          }
        } else { // 2D
          for (int n1=0; n1<2; n1++) {
            p1 = 0; p2 = 0;
            int max_off = -hnx2; int min_off = hnx2;
            int k = ks;
            for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1+1; i++) {
              int offset = off[0](k,i);
              max_off = std::max(offset+1, max_off);
              min_off = std::min(offset-1, min_off);
              p1 += std::max(offset+xgh, 0);
              p2 += std::max(1+xgh-offset, 0);
            }

            for(int kk = ks; kk<=ks+1; kk++) {
              for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1; i++) {
                int offset = off[1](k,i);
                max_off = std::max(offset+1, max_off);
                min_off = std::min(offset-1, min_off);
                p1 += std::max(offset+xgh, 0);
                p2 += std::max(1+xgh-offset, 0);
              }
            }
            ssize[0][2+n1] = std::max(max_off+xgh,0);
            ssize[1][2+n1] = std::max(1+xgh-min_off,0);
            rsize[0][2+n1] = p1;
            rsize[1][2+n1] = p2;
          }
        }
      } else if(orbital_direction == 2) { // spherical_polar
        for (int n2=0; n2<2; n2++) {
          for (int n1=0; n1<2; n1++) {
            p1 = 0; p2 = 0;
            int max_off = -hnx3; int min_off = hnx3;
            for(int j=js+n2*hnx2; j<=je-(1-n2)*hnx2  ; j++) {
              for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1+1; i++) {
                int offset = off[0](j,i);
                max_off = std::max(offset+1, max_off);
                min_off = std::min(offset-1, min_off);
                p1 += std::max(offset+xgh, 0);
                p2 += std::max(1+xgh-offset, 0);
              }
            }

            for(int j=js+n2*hnx2; j<=je-(1-n2)*hnx2+1; j++) {
              for(int i=is+n1*hnx1; i<=ie-(1-n1)*hnx1  ; i++) {
                int offset = off[1](j,i);
                max_off = std::max(offset+1, max_off);
                min_off = std::min(offset-1, min_off);
                p1 += std::max(offset+xgh, 0);
                p2 += std::max(1+xgh-offset, 0);
              }
            }
            ssize[0][2+(n1+2*n2)] = std::max(max_off+xgh,0);
            ssize[1][2+(n1+2*n2)] = std::max(1+xgh-min_off,0);
            rsize[0][2+(n1+2*n2)] = p1;
            rsize[1][2+(n1+2*n2)] = p2;
          }
        }
      }
    }
  }
  // else { // non-uniform mesh
  // }
  return;
}
