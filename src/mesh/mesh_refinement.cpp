//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh_refinement.cpp
//! \brief implements functions for static/adaptive mesh refinement

// C headers

// C++ headers
#include <algorithm>   // max()
#include <cmath>
#include <cstring>     // strcmp()
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>
#include <tuple>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../parameter_input.hpp"
#include "mesh.hpp"
#include "mesh_refinement.hpp"

//----------------------------------------------------------------------------------------
//! \fn MeshRefinement::MeshRefinement(MeshBlock *pmb, ParameterInput *pin)
//! \brief constructor

MeshRefinement::MeshRefinement(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block_(pmb), deref_count_(0),
    deref_threshold_(pin->GetOrAddInteger("mesh", "derefine_count", 10)),
    AMRFlag_(pmb->pmy_mesh->AMRFlag_) {
  // Create coarse mesh object for parent grid
  pcoarsec = new Coordinates(pmb, pin, true);

  if (NGHOST % 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshRefinement constructor" << std::endl
        << "Selected --nghost=" << NGHOST << " is incompatible with mesh refinement.\n"
        << "Reconfigure with an even number of ghost cells " << std::endl;
    ATHENA_ERROR(msg);
  }
  int qx1;
  int nc1 = pmb->ncells1;

  // In curvilinear grids, prolongation of shared face-centered fields requires face
  // area arrays, these are longer than the arrays used in the prolongation of internal
  // face fields
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0 ||
      std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    fluxinterp_ = true;
  } else {
    fluxinterp_ = false;
  }

  if (fluxinterp_) {
    qx1 = 2;
  } else {
    qx1 = 1;
  }
  fvol_[0][0].NewAthenaArray(nc1);
  fvol_[0][1].NewAthenaArray(nc1);
  fvol_[1][0].NewAthenaArray(nc1);
  fvol_[1][1].NewAthenaArray(nc1);
  sarea_x1_[0][0].NewAthenaArray(nc1+qx1);
  sarea_x1_[0][1].NewAthenaArray(nc1+qx1);
  sarea_x1_[1][0].NewAthenaArray(nc1+qx1);
  sarea_x1_[1][1].NewAthenaArray(nc1+qx1);
  sarea_x2_[0][0].NewAthenaArray(nc1);
  sarea_x2_[0][1].NewAthenaArray(nc1);
  sarea_x2_[0][2].NewAthenaArray(nc1);
  sarea_x2_[1][0].NewAthenaArray(nc1);
  sarea_x2_[1][1].NewAthenaArray(nc1);
  sarea_x2_[1][2].NewAthenaArray(nc1);
  sarea_x3_[0][0].NewAthenaArray(nc1);
  sarea_x3_[0][1].NewAthenaArray(nc1);
  sarea_x3_[1][0].NewAthenaArray(nc1);
  sarea_x3_[1][1].NewAthenaArray(nc1);
  sarea_x3_[2][0].NewAthenaArray(nc1);
  sarea_x3_[2][1].NewAthenaArray(nc1);

  // Create coarse area arrays used in prolongation of shared face-centered fields in
  // curvilinear grids
  if (fluxinterp_) {
    csarea_x1_.NewAthenaArray(nc1);
    csarea_x2_.NewAthenaArray(nc1);
    csarea_x3_.NewAthenaArray(nc1);
  }

  // KGF: probably don't need to preallocate space for pointers in these vectors
  pvars_cc_.reserve(3);
  pvars_fc_.reserve(3);
}


//----------------------------------------------------------------------------------------
//! \fn MeshRefinement::~MeshRefinement()
//! \brief destructor

MeshRefinement::~MeshRefinement() {
  delete pcoarsec;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::RestrictCellCenteredValues(const AthenaArray<Real> &fine,
//!                          AthenaArray<Real> &coarse, int sn, int en,
//!                          int csi, int cei, int csj, int cej, int csk, int cek)
//! \brief restrict cell centered values

void MeshRefinement::RestrictCellCenteredValues(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse, int sn, int en,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int si = (csi - pmb->cis)*2 + pmb->is, ei = (cei - pmb->cis)*2 + pmb->is + 1;

  // store the restricted data in the prolongation buffer for later use
  if (pmb->block_size.nx3>1) { // 3D
    for (int n=sn; n<=en; ++n) {
      for (int ck=csk; ck<=cek; ck++) {
        int k = (ck - pmb->cks)*2 + pmb->ks;
        for (int cj=csj; cj<=cej; cj++) {
          int j = (cj - pmb->cjs)*2 + pmb->js;
          pco->CellVolume(k,j,si,ei,fvol_[0][0]);
          pco->CellVolume(k,j+1,si,ei,fvol_[0][1]);
          pco->CellVolume(k+1,j,si,ei,fvol_[1][0]);
          pco->CellVolume(k+1,j+1,si,ei,fvol_[1][1]);
          for (int ci=csi; ci<=cei; ci++) {
            int i = (ci - pmb->cis)*2 + pmb->is;
            // KGF: add the off-centered quantities first to preserve FP symmetry
            Real tvol = ((fvol_[0][0](i) + fvol_[0][1](i))
                         + (fvol_[0][0](i+1) + fvol_[0][1](i+1)))
                        + ((fvol_[1][0](i) + fvol_[1][1](i))
                           + (fvol_[1][0](i+1) + fvol_[1][1](i+1)));
            // KGF: add the off-centered quantities first to preserve FP symmetry
            coarse(n,ck,cj,ci) =
                (((fine(n,k  ,j  ,i)*fvol_[0][0](i) + fine(n,k  ,j+1,i)*fvol_[0][1](i))
                  + (fine(n,k  ,j  ,i+1)*fvol_[0][0](i+1) +
                     fine(n,k  ,j+1,i+1)*fvol_[0][1](i+1)))
                 + ((fine(n,k+1,j  ,i)*fvol_[1][0](i) + fine(n,k+1,j+1,i)*fvol_[1][1](i))
                    + (fine(n,k+1,j  ,i+1)*fvol_[1][0](i+1) +
                       fine(n,k+1,j+1,i+1)*fvol_[1][1](i+1)))) / tvol;
          }
        }
      }
    }
  } else if (pmb->block_size.nx2>1) { // 2D
    for (int n=sn; n<=en; ++n) {
      for (int cj=csj; cj<=cej; cj++) {
        int j = (cj - pmb->cjs)*2 + pmb->js;
        pco->CellVolume(0,j  ,si,ei,fvol_[0][0]);
        pco->CellVolume(0,j+1,si,ei,fvol_[0][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i = (ci - pmb->cis)*2 + pmb->is;
          // KGF: add the off-centered quantities first to preserve FP symmetry
          Real tvol = (fvol_[0][0](i) + fvol_[0][1](i)) +
                      (fvol_[0][0](i+1) + fvol_[0][1](i+1));

          // KGF: add the off-centered quantities first to preserve FP symmetry
          coarse(n,0,cj,ci) =
              ((fine(n,0,j  ,i)*fvol_[0][0](i) + fine(n,0,j+1,i)*fvol_[0][1](i))
               + (fine(n,0,j ,i+1)*fvol_[0][0](i+1) + fine(n,0,j+1,i+1)*fvol_[0][1](i+1)))
              /tvol;
        }
      }
    }
  } else { // 1D
    int j = pmb->js, cj = pmb->cjs, k = pmb->ks, ck = pmb->cks;
    for (int n=sn; n<=en; ++n) {
      pco->CellVolume(k,j,si,ei,fvol_[0][0]);
      for (int ci=csi; ci<=cei; ci++) {
        int i = (ci - pmb->cis)*2 + pmb->is;
        Real tvol = fvol_[0][0](i) + fvol_[0][0](i+1);
        coarse(n,ck,cj,ci)
            = (fine(n,k,j,i)*fvol_[0][0](i) + fine(n,k,j,i+1)*fvol_[0][0](i+1))/tvol;
      }
    }
  }
}


// over load function  for radiation variables

void MeshRefinement::RestrictCellCenteredValues(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse, int array_order,
    int sn, int en, int csi, int cei, int csj, int cej, int csk, int cek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int si = (csi - pmb->cis)*2 + pmb->is, ei = (cei - pmb->cis)*2 + pmb->is + 1;

  // reverse order for radiation variables
  if (array_order < 0) {
    // store the restricted data in the prolongation buffer for later use
    if (pmb->block_size.nx3>1) { // 3D
      for (int ck=csk; ck<=cek; ck++) {
        int k = (ck - pmb->cks)*2 + pmb->ks;
        for (int cj=csj; cj<=cej; cj++) {
          int j = (cj - pmb->cjs)*2 + pmb->js;
          pco->CellVolume(k,j,si,ei,fvol_[0][0]);
          pco->CellVolume(k,j+1,si,ei,fvol_[0][1]);
          pco->CellVolume(k+1,j,si,ei,fvol_[1][0]);
          pco->CellVolume(k+1,j+1,si,ei,fvol_[1][1]);
          for (int ci=csi; ci<=cei; ci++) {
            int i = (ci - pmb->cis)*2 + pmb->is;
            Real tvol = ((fvol_[0][0](i) + fvol_[0][1](i))
                      + (fvol_[0][0](i+1) + fvol_[0][1](i+1)))
                      + ((fvol_[1][0](i) + fvol_[1][1](i))
                      + (fvol_[1][0](i+1) + fvol_[1][1](i+1)));

            for (int n=sn; n<=en; ++n) {
              // KGF: add the off-centered quantities first to preserve FP symmetry
              // KGF: add the off-centered quantities first to preserve FP symmetry
              coarse(ck,cj,ci,n) =
                  (((fine(k  ,j  ,i,n)*fvol_[0][0](i) + fine(k  ,j+1,i,n)*fvol_[0][1](i))
                    + (fine(k  ,j  ,i+1,n)*fvol_[0][0](i+1) +
                       fine(k  ,j+1,i+1,n)*fvol_[0][1](i+1)))
                   + ((fine(k+1,j  ,i,n)*fvol_[1][0](i)
                       + fine(k+1,j+1,i,n)*fvol_[1][1](i))
                      + (fine(k+1,j  ,i+1,n)*fvol_[1][0](i+1) +
                         fine(k+1,j+1,i+1,n)*fvol_[1][1](i+1)))) / tvol;
            }
          }// end i
        }
      }
    } else if (pmb->block_size.nx2>1) { // 2D
      for (int cj=csj; cj<=cej; cj++) {
        int j = (cj - pmb->cjs)*2 + pmb->js;
        pco->CellVolume(0,j  ,si,ei,fvol_[0][0]);
        pco->CellVolume(0,j+1,si,ei,fvol_[0][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i = (ci - pmb->cis)*2 + pmb->is;
            // KGF: add the off-centered quantities first to preserve FP symmetry
          Real tvol = (fvol_[0][0](i) + fvol_[0][1](i)) +
                        (fvol_[0][0](i+1) + fvol_[0][1](i+1));
          for (int n=sn; n<=en; ++n) {
            // KGF: add the off-centered quantities first to preserve FP symmetry
            coarse(0,cj,ci,n) =
                ((fine(0,j  ,i,n)*fvol_[0][0](i) + fine(0,j+1,i,n)*fvol_[0][1](i))
                 + (fine(0,j ,i+1,n)*fvol_[0][0](i+1)
                    + fine(0,j+1,i+1,n)*fvol_[0][1](i+1)))
                /tvol;
          }
        }
      }
    } else { // 1D
      int j = pmb->js, cj = pmb->cjs, k = pmb->ks, ck = pmb->cks;

      pco->CellVolume(k,j,si,ei,fvol_[0][0]);
      for (int ci=csi; ci<=cei; ci++) {
        int i = (ci - pmb->cis)*2 + pmb->is;
        Real tvol = fvol_[0][0](i) + fvol_[0][0](i+1);
        for (int n=sn; n<=en; ++n) {
          coarse(ck,cj,ci,n)
              = (fine(k,j,i,n)*fvol_[0][0](i) + fine(k,j,i+1,n)*fvol_[0][0](i+1))/tvol;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::RestrictFieldX1(const AthenaArray<Real> &fine
//!     AthenaArray<Real> &coarse, int csi, int cei, int csj, int cej, int csk, int cek)
//! \brief restrict the x1 field data and set them into the coarse buffer

void MeshRefinement::RestrictFieldX1(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int si = (csi - pmb->cis)*2 + pmb->is, ei = (cei - pmb->cis)*2 + pmb->is;

  // store the restricted data in the prolongation buffer for later use
  if (pmb->block_size.nx3>1) { // 3D
    for (int ck=csk; ck<=cek; ck++) {
      int k = (ck - pmb->cks)*2 + pmb->ks;
      for (int cj=csj; cj<=cej; cj++) {
        int j = (cj - pmb->cjs)*2 + pmb->js;
        pco->Face1Area(k,   j,   si, ei, sarea_x1_[0][0]);
        pco->Face1Area(k,   j+1, si, ei, sarea_x1_[0][1]);
        pco->Face1Area(k+1, j,   si, ei, sarea_x1_[1][0]);
        pco->Face1Area(k+1, j+1, si, ei, sarea_x1_[1][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i = (ci - pmb->cis)*2 + pmb->is;
          Real tarea = sarea_x1_[0][0](i) + sarea_x1_[0][1](i) +
                       sarea_x1_[1][0](i) + sarea_x1_[1][1](i);
          coarse(ck,cj,ci) =
              (fine(k  ,j,i)*sarea_x1_[0][0](i) + fine(k  ,j+1,i)*sarea_x1_[0][1](i)
               + fine(k+1,j,i)*sarea_x1_[1][0](i) + fine(k+1,j+1,i)*sarea_x1_[1][1](i)
               )/tarea;
        }
      }
    }
  } else if (pmb->block_size.nx2>1) { // 2D
    int k = pmb->ks;
    for (int cj=csj; cj<=cej; cj++) {
      int j = (cj - pmb->cjs)*2 + pmb->js;
      pco->Face1Area(k,  j,   si, ei, sarea_x1_[0][0]);
      pco->Face1Area(k,  j+1, si, ei, sarea_x1_[0][1]);
      for (int ci=csi; ci<=cei; ci++) {
        int i = (ci - pmb->cis)*2 + pmb->is;
        Real tarea = sarea_x1_[0][0](i) + sarea_x1_[0][1](i);
        coarse(csk,cj,ci) =
            (fine(k,j,i)*sarea_x1_[0][0](i) + fine(k,j+1,i)*sarea_x1_[0][1](i))/tarea;
      }
    }
  } else { // 1D - no restriction, just copy
    for (int ci=csi; ci<=cei; ci++) {
      int i = (ci - pmb->cis)*2 + pmb->is;
      coarse(csk,csj,ci) = fine(pmb->ks,pmb->js,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::RestrictFieldX2(const AthenaArray<Real> &fine
//!     AthenaArray<Real> &coarse, int csi, int cei, int csj, int cej, int csk, int cek)
//! \brief restrict the x2 field data and set them into the coarse buffer

void MeshRefinement::RestrictFieldX2(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int si = (csi - pmb->cis)*2 + pmb->is, ei = (cei - pmb->cis)*2 + pmb->is + 1;

  // store the restricted data in the prolongation buffer for later use
  if (pmb->block_size.nx3>1) { // 3D
    for (int ck=csk; ck<=cek; ck++) {
      int k = (ck - pmb->cks)*2 + pmb->ks;
      for (int cj=csj; cj<=cej; cj++) {
        int j = (cj - pmb->cjs)*2 + pmb->js;
        bool pole = pco->IsPole(j);
        if (!pole) {
          pco->Face2Area(k,   j,  si, ei, sarea_x2_[0][0]);
          pco->Face2Area(k+1, j,  si, ei, sarea_x2_[1][0]);
        } else {
          for (int ci = csi; ci <= cei; ++ci) {
            int i = (ci - pmb->cis) * 2 + pmb->is;
            sarea_x2_[0][0](i) = pco->dx1f(i);
            sarea_x2_[1][0](i) = pco->dx1f(i);
          }
        }
        for (int ci=csi; ci<=cei; ci++) {
          int i = (ci - pmb->cis)*2 + pmb->is;
          Real tarea = sarea_x2_[0][0](i) + sarea_x2_[0][0](i+1) +
                       sarea_x2_[1][0](i) + sarea_x2_[1][0](i+1);
          coarse(ck,cj,ci) =
              (fine(k  ,j,i)*sarea_x2_[0][0](i) + fine(k  ,j,i+1)*sarea_x2_[0][0](i+1)
               +fine(k+1,j,i)*sarea_x2_[1][0](i) + fine(k+1,j,i+1)*sarea_x2_[1][0](i+1))
              /tarea;
        }
      }
    }
  } else if (pmb->block_size.nx2>1) { // 2D
    int k = pmb->ks;
    for (int cj=csj; cj<=cej; cj++) {
      int j = (cj - pmb->cjs)*2 + pmb->js;
      bool pole = pco->IsPole(j);
      if (!pole) {
        pco->Face2Area(k, j, si, ei, sarea_x2_[0][0]);
      } else {
        for (int ci = csi; ci <= cei; ++ci) {
          int i = (ci - pmb->cis) * 2 + pmb->is;
          sarea_x2_[0][0](i) = pco->dx1f(i);
        }
      }
      for (int ci=csi; ci<=cei; ci++) {
        int i = (ci - pmb->cis)*2 + pmb->is;
        Real tarea = sarea_x2_[0][0](i) + sarea_x2_[0][0](i+1);
        coarse(pmb->cks,cj,ci) =
            (fine(k,j,i)*sarea_x2_[0][0](i) + fine(k,j,i+1)*sarea_x2_[0][0](i+1))/tarea;
      }
    }
  } else { // 1D
    int k = pmb->ks, j = pmb->js;
    pco->Face2Area(k, j, si, ei, sarea_x2_[0][0]);
    for (int ci=csi; ci<=cei; ci++) {
      int i = (ci - pmb->cis)*2 + pmb->is;
      Real tarea = sarea_x2_[0][0](i) + sarea_x2_[0][0](i+1);
      coarse(pmb->cks,pmb->cjs,ci) =
          (fine(k,j,i)*sarea_x2_[0][0](i) + fine(k,j,i+1)*sarea_x2_[0][0](i+1))/tarea;
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::RestrictFieldX3(const AthenaArray<Real> &fine
//!     AthenaArray<Real> &coarse, int csi, int cei, int csj, int cej, int csk, int cek)
//! \brief restrict the x3 field data and set them into the coarse buffer

void MeshRefinement::RestrictFieldX3(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int si = (csi - pmb->cis)*2 + pmb->is, ei = (cei - pmb->cis)*2 + pmb->is + 1;

  // store the restricted data in the prolongation buffer for later use
  if (pmb->block_size.nx3>1) { // 3D
    for (int ck=csk; ck<=cek; ck++) {
      int k = (ck - pmb->cks)*2 + pmb->ks;
      for (int cj=csj; cj<=cej; cj++) {
        int j = (cj - pmb->cjs)*2 + pmb->js;
        pco->Face3Area(k,   j,  si, ei, sarea_x3_[0][0]);
        pco->Face3Area(k, j+1,  si, ei, sarea_x3_[0][1]);
        for (int ci=csi; ci<=cei; ci++) {
          int i = (ci - pmb->cis)*2 + pmb->is;
          Real tarea = sarea_x3_[0][0](i) + sarea_x3_[0][0](i+1) +
                       sarea_x3_[0][1](i) + sarea_x3_[0][1](i+1);
          coarse(ck,cj,ci)  =
              (fine(k,j  ,i)*sarea_x3_[0][0](i) + fine(k,j  ,i+1)*sarea_x3_[0][0](i+1)
               + fine(k,j+1,i)*sarea_x3_[0][1](i) + fine(k,j+1,i+1)*sarea_x3_[0][1](i+1)
               ) /tarea;
        }
      }
    }
  } else if (pmb->block_size.nx2>1) { // 2D
    int k = pmb->ks;
    for (int cj=csj; cj<=cej; cj++) {
      int j = (cj - pmb->cjs)*2 + pmb->js;
      pco->Face3Area(k,   j, si, ei, sarea_x3_[0][0]);
      pco->Face3Area(k, j+1, si, ei, sarea_x3_[0][1]);
      for (int ci=csi; ci<=cei; ci++) {
        int i = (ci - pmb->cis)*2 + pmb->is;
        Real tarea = sarea_x3_[0][0](i) + sarea_x3_[0][0](i+1) +
                     sarea_x3_[0][1](i) + sarea_x3_[0][1](i+1);
        coarse(pmb->cks,cj,ci) =
            (fine(k,j  ,i)*sarea_x3_[0][0](i) + fine(k,j  ,i+1)*sarea_x3_[0][0](i+1)
             + fine(k,j+1,i)*sarea_x3_[0][1](i) + fine(k,j+1,i+1)*sarea_x3_[0][1](i+1)
             ) /tarea;
      }
    }
  } else { // 1D
    int k = pmb->ks, j = pmb->js;
    pco->Face3Area(k, j, si, ei, sarea_x3_[0][0]);
    for (int ci=csi; ci<=cei; ci++) {
      int i = (ci - pmb->cis)*2 + pmb->is;
      Real tarea = sarea_x3_[0][0](i) + sarea_x3_[0][0](i+1);
      coarse(pmb->cks,pmb->cjs,ci) =
          (fine(k,j,i)*sarea_x3_[0][0](i) + fine(k,j,i+1)*sarea_x3_[0][0](i+1))/tarea;
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateCellCenteredValues(
//!       const AthenaArray<Real> &coarse,AthenaArray<Real> &fine, int sn, int en,
//!       int si, int ei, int sj, int ej, int sk, int ek)
//! \brief Prolongate cell centered values

void MeshRefinement::ProlongateCellCenteredValues(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int sn, int en, int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  if (pmb->block_size.nx3 > 1) {
    for (int n=sn; n<=en; n++) {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        const Real& x3m = pcoarsec->x3v(k-1);
        const Real& x3c = pcoarsec->x3v(k);
        const Real& x3p = pcoarsec->x3v(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        const Real& fx3m = pco->x3v(fk);
        const Real& fx3p = pco->x3v(fk+1);
        Real dx3fm =  x3c - fx3m;
        Real dx3fp =  fx3p - x3c;
        for (int j = sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          const Real& x2m = pcoarsec->x2v(j-1);
          const Real& x2c = pcoarsec->x2v(j);
          const Real& x2p = pcoarsec->x2v(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          const Real& fx2m = pco->x2v(fj);
          const Real& fx2p = pco->x2v(fj+1);
          Real dx2fm = x2c - fx2m;
          Real dx2fp = fx2p - x2c;
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            const Real& x1m = pcoarsec->x1v(i-1);
            const Real& x1c = pcoarsec->x1v(i);
            const Real& x1p = pcoarsec->x1v(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            const Real& fx1m = pco->x1v(fi);
            const Real& fx1p = pco->x1v(fi+1);
            Real dx1fm = x1c - fx1m;
            Real dx1fp = fx1p - x1c;
            Real ccval = coarse(n,k,j,i);

            // calculate 3D gradients using the minmod limiter
            Real gx1m = (ccval - coarse(n,k,j,i-1))/dx1m;
            Real gx1p = (coarse(n,k,j,i+1) - ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*
                        std::min(std::abs(gx1m), std::abs(gx1p));
            Real gx2m = (ccval - coarse(n,k,j-1,i))/dx2m;
            Real gx2p = (coarse(n,k,j+1,i) - ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*
                        std::min(std::abs(gx2m), std::abs(gx2p));
            Real gx3m = (ccval - coarse(n,k-1,j,i))/dx3m;
            Real gx3p = (coarse(n,k+1,j,i) - ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*
                        std::min(std::abs(gx3m), std::abs(gx3p));

            // KGF: add the off-centered quantities first to preserve FP symmetry
            // interpolate onto the finer grid
            fine(n,fk  ,fj  ,fi  ) = ccval - (gx1c*dx1fm + gx2c*dx2fm + gx3c*dx3fm);
            fine(n,fk  ,fj  ,fi+1) = ccval + (gx1c*dx1fp - gx2c*dx2fm - gx3c*dx3fm);
            fine(n,fk  ,fj+1,fi  ) = ccval - (gx1c*dx1fm - gx2c*dx2fp + gx3c*dx3fm);
            fine(n,fk  ,fj+1,fi+1) = ccval + (gx1c*dx1fp + gx2c*dx2fp - gx3c*dx3fm);
            fine(n,fk+1,fj  ,fi  ) = ccval - (gx1c*dx1fm + gx2c*dx2fm - gx3c*dx3fp);
            fine(n,fk+1,fj  ,fi+1) = ccval + (gx1c*dx1fp - gx2c*dx2fm + gx3c*dx3fp);
            fine(n,fk+1,fj+1,fi  ) = ccval - (gx1c*dx1fm - gx2c*dx2fp - gx3c*dx3fp);
            fine(n,fk+1,fj+1,fi+1) = ccval + (gx1c*dx1fp + gx2c*dx2fp + gx3c*dx3fp);
          }
        }
      }
    }
  } else if (pmb->block_size.nx2 > 1) {
    int k = pmb->cks, fk = pmb->ks;
    for (int n=sn; n<=en; n++) {
      for (int j=sj; j<=ej; j++) {
        int fj = (j - pmb->cjs)*2 + pmb->js;
        const Real& x2m = pcoarsec->x2v(j-1);
        const Real& x2c = pcoarsec->x2v(j);
        const Real& x2p = pcoarsec->x2v(j+1);
        Real dx2m = x2c - x2m;
        Real dx2p = x2p - x2c;
        const Real& fx2m = pco->x2v(fj);
        const Real& fx2p = pco->x2v(fj+1);
        Real dx2fm = x2c - fx2m;
        Real dx2fp = fx2p - x2c;
        for (int i=si; i<=ei; i++) {
          int fi = (i - pmb->cis)*2 + pmb->is;
          const Real& x1m = pcoarsec->x1v(i-1);
          const Real& x1c = pcoarsec->x1v(i);
          const Real& x1p = pcoarsec->x1v(i+1);
          Real dx1m = x1c - x1m;
          Real dx1p = x1p - x1c;
          const Real& fx1m = pco->x1v(fi);
          const Real& fx1p = pco->x1v(fi+1);
          Real dx1fm = x1c - fx1m;
          Real dx1fp = fx1p - x1c;
          Real ccval = coarse(n,k,j,i);

          // calculate 2D gradients using the minmod limiter
          Real gx1m = (ccval - coarse(n,k,j,i-1))/dx1m;
          Real gx1p = (coarse(n,k,j,i+1) - ccval)/dx1p;
          Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*
                      std::min(std::abs(gx1m), std::abs(gx1p));
          Real gx2m = (ccval - coarse(n,k,j-1,i))/dx2m;
          Real gx2p = (coarse(n,k,j+1,i) - ccval)/dx2p;
          Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*
                      std::min(std::abs(gx2m), std::abs(gx2p));

          // KGF: add the off-centered quantities first to preserve FP symmetry
          // interpolate onto the finer grid
          fine(n,fk  ,fj  ,fi  ) = ccval - (gx1c*dx1fm + gx2c*dx2fm);
          fine(n,fk  ,fj  ,fi+1) = ccval + (gx1c*dx1fp - gx2c*dx2fm);
          fine(n,fk  ,fj+1,fi  ) = ccval - (gx1c*dx1fm - gx2c*dx2fp);
          fine(n,fk  ,fj+1,fi+1) = ccval + (gx1c*dx1fp + gx2c*dx2fp);
        }
      }
    }
  } else { // 1D
    int k = pmb->cks, fk = pmb->ks, j = pmb->cjs, fj = pmb->js;
    for (int n=sn; n<=en; n++) {
      for (int i=si; i<=ei; i++) {
        int fi = (i - pmb->cis)*2 + pmb->is;
        const Real& x1m = pcoarsec->x1v(i-1);
        const Real& x1c = pcoarsec->x1v(i);
        const Real& x1p = pcoarsec->x1v(i+1);
        Real dx1m = x1c - x1m;
        Real dx1p = x1p - x1c;
        const Real& fx1m = pco->x1v(fi);
        const Real& fx1p = pco->x1v(fi+1);
        Real dx1fm = x1c - fx1m;
        Real dx1fp = fx1p - x1c;
        Real ccval = coarse(n,k,j,i);

        // calculate 1D gradient using the min-mod limiter
        Real gx1m = (ccval - coarse(n,k,j,i-1))/dx1m;
        Real gx1p = (coarse(n,k,j,i+1) - ccval)/dx1p;
        Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                           std::abs(gx1p));

        // interpolate on to the finer grid
        fine(n,fk  ,fj  ,fi  ) = ccval - gx1c*dx1fm;
        fine(n,fk  ,fj  ,fi+1) = ccval + gx1c*dx1fp;
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateCellCenteredValues(
//        const AthenaArray<Real> &coarse,AthenaArray<Real> &fine,
//      int array_order, int sn, int en,,
//        int si, int ei, int sj, int ej, int sk, int ek)
//  \brief Prolongate cell centered values for radiation variables

void MeshRefinement::ProlongateCellCenteredValues(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int array_order,
    int sn, int en, int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;

  if (array_order < 0) {
    if (pmb->block_size.nx3 > 1) {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        const Real& x3m = pcoarsec->x3v(k-1);
        const Real& x3c = pcoarsec->x3v(k);
        const Real& x3p = pcoarsec->x3v(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        const Real& fx3m = pco->x3v(fk);
        const Real& fx3p = pco->x3v(fk+1);
        Real dx3fm =  x3c - fx3m;
        Real dx3fp =  fx3p - x3c;
        for (int j = sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          const Real& x2m = pcoarsec->x2v(j-1);
          const Real& x2c = pcoarsec->x2v(j);
          const Real& x2p = pcoarsec->x2v(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          const Real& fx2m = pco->x2v(fj);
          const Real& fx2p = pco->x2v(fj+1);
          Real dx2fm = x2c - fx2m;
          Real dx2fp = fx2p - x2c;
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            const Real& x1m = pcoarsec->x1v(i-1);
            const Real& x1c = pcoarsec->x1v(i);
            const Real& x1p = pcoarsec->x1v(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            const Real& fx1m = pco->x1v(fi);
            const Real& fx1p = pco->x1v(fi+1);
            Real dx1fm = x1c - fx1m;
            Real dx1fp = fx1p - x1c;

            for (int n=sn; n<=en; n++) {
              Real ccval = coarse(k,j,i,n);
              // calculate 3D gradients using the minmod limiter
              Real gx1m = (ccval - coarse(k,j,i-1,n))/dx1m;
              Real gx1p = (coarse(k,j,i+1,n) - ccval)/dx1p;
              Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*
                          std::min(std::abs(gx1m), std::abs(gx1p));
              Real gx2m = (ccval - coarse(k,j-1,i,n))/dx2m;
              Real gx2p = (coarse(k,j+1,i,n) - ccval)/dx2p;
              Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*
                          std::min(std::abs(gx2m), std::abs(gx2p));
              Real gx3m = (ccval - coarse(k-1,j,i,n))/dx3m;
              Real gx3p = (coarse(k+1,j,i,n) - ccval)/dx3p;
              Real gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*
                          std::min(std::abs(gx3m), std::abs(gx3p));

              // KGF: add the off-centered quantities first to preserve FP symmetry
              // interpolate onto the finer grid
              fine(fk  ,fj  ,fi  ,n) = ccval - (gx1c*dx1fm + gx2c*dx2fm + gx3c*dx3fm);
              fine(fk  ,fj  ,fi+1,n) = ccval + (gx1c*dx1fp - gx2c*dx2fm - gx3c*dx3fm);
              fine(fk  ,fj+1,fi  ,n) = ccval - (gx1c*dx1fm - gx2c*dx2fp + gx3c*dx3fm);
              fine(fk  ,fj+1,fi+1,n) = ccval + (gx1c*dx1fp + gx2c*dx2fp - gx3c*dx3fm);
              fine(fk+1,fj  ,fi  ,n) = ccval - (gx1c*dx1fm + gx2c*dx2fm - gx3c*dx3fp);
              fine(fk+1,fj  ,fi+1,n) = ccval + (gx1c*dx1fp - gx2c*dx2fm + gx3c*dx3fp);
              fine(fk+1,fj+1,fi  ,n) = ccval - (gx1c*dx1fm - gx2c*dx2fp - gx3c*dx3fp);
              fine(fk+1,fj+1,fi+1,n) = ccval + (gx1c*dx1fp + gx2c*dx2fp + gx3c*dx3fp);
            }
          }
        }
      }
    } else if (pmb->block_size.nx2 > 1) {
      int k = pmb->cks, fk = pmb->ks;
      for (int j=sj; j<=ej; j++) {
        int fj = (j - pmb->cjs)*2 + pmb->js;
        const Real& x2m = pcoarsec->x2v(j-1);
        const Real& x2c = pcoarsec->x2v(j);
        const Real& x2p = pcoarsec->x2v(j+1);
        Real dx2m = x2c - x2m;
        Real dx2p = x2p - x2c;
        const Real& fx2m = pco->x2v(fj);
        const Real& fx2p = pco->x2v(fj+1);
        Real dx2fm = x2c - fx2m;
        Real dx2fp = fx2p - x2c;
        for (int i=si; i<=ei; i++) {
          int fi = (i - pmb->cis)*2 + pmb->is;
          const Real& x1m = pcoarsec->x1v(i-1);
          const Real& x1c = pcoarsec->x1v(i);
          const Real& x1p = pcoarsec->x1v(i+1);
          Real dx1m = x1c - x1m;
          Real dx1p = x1p - x1c;
          const Real& fx1m = pco->x1v(fi);
          const Real& fx1p = pco->x1v(fi+1);
          Real dx1fm = x1c - fx1m;
          Real dx1fp = fx1p - x1c;
          for (int n=sn; n<=en; n++) {
            Real ccval = coarse(k,j,i,n);

            // calculate 2D gradients using the minmod limiter
            Real gx1m = (ccval - coarse(k,j,i-1,n))/dx1m;
            Real gx1p = (coarse(k,j,i+1,n) - ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*
                        std::min(std::abs(gx1m), std::abs(gx1p));
            Real gx2m = (ccval - coarse(k,j-1,i,n))/dx2m;
            Real gx2p = (coarse(k,j+1,i,n) - ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*
                        std::min(std::abs(gx2m), std::abs(gx2p));

            // KGF: add the off-centered quantities first to preserve FP symmetry
            // interpolate onto the finer grid
            fine(fk  ,fj  ,fi  ,n) = ccval - (gx1c*dx1fm + gx2c*dx2fm);
            fine(fk  ,fj  ,fi+1,n) = ccval + (gx1c*dx1fp - gx2c*dx2fm);
            fine(fk  ,fj+1,fi  ,n) = ccval - (gx1c*dx1fm - gx2c*dx2fp);
            fine(fk  ,fj+1,fi+1,n) = ccval + (gx1c*dx1fp + gx2c*dx2fp);
          }
        }
      }
    } else { // 1D
      int k = pmb->cks, fk = pmb->ks, j = pmb->cjs, fj = pmb->js;

      for (int i=si; i<=ei; i++) {
        int fi = (i - pmb->cis)*2 + pmb->is;
        const Real& x1m = pcoarsec->x1v(i-1);
        const Real& x1c = pcoarsec->x1v(i);
        const Real& x1p = pcoarsec->x1v(i+1);
        Real dx1m = x1c - x1m;
        Real dx1p = x1p - x1c;
        const Real& fx1m = pco->x1v(fi);
        const Real& fx1p = pco->x1v(fi+1);
        Real dx1fm = x1c - fx1m;
        Real dx1fp = fx1p - x1c;
        for (int n=sn; n<=en; n++) {
          Real ccval = coarse(k,j,i,n);
          // calculate 1D gradient using the min-mod limiter
          Real gx1m = (ccval - coarse(k,j,i-1,n))/dx1m;
          Real gx1p = (coarse(k,j,i+1,n) - ccval)/dx1p;
          Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                             std::abs(gx1p));
          // interpolate on to the finer grid
          fine(fk  ,fj  ,fi  ,n) = ccval - gx1c*dx1fm;
          fine(fk  ,fj  ,fi+1,n) = ccval + gx1c*dx1fp;
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateSharedFieldX1(const AthenaArray<Real> &coarse,
//!     AthenaArray<Real> &fine, int si, int ei, int sj, int ej, int sk, int ek)
//! \brief prolongate x1 face-centered fields shared between coarse and fine levels

void MeshRefinement::ProlongateSharedFieldX1(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int fsi = (si - pmb->cis)*2 + pmb->is, fei = (ei - pmb->cis)*2 + pmb->is + 1;
  if (pmb->block_size.nx3 > 1) {
    if (fluxinterp_) {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        const Real& x3m = pcoarsec->x3s1(k-1);
        const Real& x3c = pcoarsec->x3s1(k);
        const Real& x3p = pcoarsec->x3s1(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        const Real& fx3m = pco->x3s1(fk);
        const Real& fx3p = pco->x3s1(fk+1);
        Real dfx3 = fx3p - fx3m;
        for (int j=sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          const Real& x2m = pcoarsec->x2s1(j-1);
          const Real& x2c = pcoarsec->x2s1(j);
          const Real& x2p = pcoarsec->x2s1(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          const Real& fx2m = pco->x2s1(fj);
          const Real& fx2p = pco->x2s1(fj+1);
          Real dfx2 = fx2p - fx2m;
          pco->Face1Area(fk,   fj,   fsi, fei, sarea_x1_[0][0]);
          pco->Face1Area(fk,   fj+1, fsi, fei, sarea_x1_[0][1]);
          pco->Face1Area(fk+1, fj,   fsi, fei, sarea_x1_[1][0]);
          pco->Face1Area(fk+1, fj+1, fsi, fei, sarea_x1_[1][1]);
          pcoarsec->Face1Area(k, j, si, ei, csarea_x1_);
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            Real ccval = coarse(k,j,i);
            Real csarea = csarea_x1_(i);
            Real fsa00 = sarea_x1_[0][0](fi);
            Real fsa01 = sarea_x1_[0][1](fi);
            Real fsa10 = sarea_x1_[1][0](fi);
            Real fsa11 = sarea_x1_[1][1](fi);

            Real gx2m = (ccval - coarse(k,j-1,i))/dx2m;
            Real gx2p = (coarse(k,j+1,i) - ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*std::min(std::abs(gx2m),
                                                               std::abs(gx2p));
            Real gx3m = (ccval - coarse(k-1,j,i))/dx3m;
            Real gx3p = (coarse(k+1,j,i) - ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*std::min(std::abs(gx3m),
                                                               std::abs(gx3p));
            fine(fk  ,fj  ,fi) = ccval - gx2c*dfx2*(fsa01 + fsa11)/csarea
                                       - gx3c*dfx3*(fsa10 + fsa11)/csarea;
            fine(fk  ,fj+1,fi) = ccval + gx2c*dfx2*(fsa00 + fsa10)/csarea
                                       - gx3c*dfx3*(fsa10 + fsa11)/csarea;
            fine(fk+1,fj  ,fi) = ccval - gx2c*dfx2*(fsa01 + fsa11)/csarea
                                       + gx3c*dfx3*(fsa00 + fsa01)/csarea;
            fine(fk+1,fj+1,fi) = ccval + gx2c*dfx2*(fsa00 + fsa10)/csarea
                                       + gx3c*dfx3*(fsa00 + fsa01)/csarea;
          }
        }
      }
    } else {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        const Real& x3m = pcoarsec->x3s1(k-1);
        const Real& x3c = pcoarsec->x3s1(k);
        const Real& x3p = pcoarsec->x3s1(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        const Real& fx3m = pco->x3s1(fk);
        const Real& fx3p = pco->x3s1(fk+1);
        for (int j=sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          const Real& x2m = pcoarsec->x2s1(j-1);
          const Real& x2c = pcoarsec->x2s1(j);
          const Real& x2p = pcoarsec->x2s1(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          const Real& fx2m = pco->x2s1(fj);
          const Real& fx2p = pco->x2s1(fj+1);
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            Real ccval = coarse(k,j,i);

            Real gx2m = (ccval - coarse(k,j-1,i))/dx2m;
            Real gx2p = (coarse(k,j+1,i) - ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*std::min(std::abs(gx2m),
                                                             std::abs(gx2p));
            Real gx3m = (ccval - coarse(k-1,j,i))/dx3m;
            Real gx3p = (coarse(k+1,j,i) - ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*std::min(std::abs(gx3m),
                                                             std::abs(gx3p));
            fine(fk  ,fj  ,fi) = ccval - gx2c*(x2c - fx2m) - gx3c*(x3c - fx3m);
            fine(fk  ,fj+1,fi) = ccval + gx2c*(fx2p - x2c) - gx3c*(x3c - fx3m);
            fine(fk+1,fj  ,fi) = ccval - gx2c*(x2c - fx2m) + gx3c*(fx3p - x3c);
            fine(fk+1,fj+1,fi) = ccval + gx2c*(fx2p - x2c) + gx3c*(fx3p - x3c);
          }
        }
      }
    }
  } else if (pmb->block_size.nx2 > 1) {
    int k = pmb->cks, fk = pmb->ks;
    for (int j=sj; j<=ej; j++) {
      int fj = (j - pmb->cjs)*2 + pmb->js;
      const Real& x2m = pcoarsec->x2s1(j-1);
      const Real& x2c = pcoarsec->x2s1(j);
      const Real& x2p = pcoarsec->x2s1(j+1);
      Real dx2m = x2c - x2m;
      Real dx2p = x2p - x2c;
      const Real& fx2m = pco->x2s1(fj);
      const Real& fx2p = pco->x2s1(fj+1);
      for (int i=si; i<=ei; i++) {
        int fi = (i - pmb->cis)*2 + pmb->is;
        Real ccval = coarse(k,j,i);

        Real gx2m = (ccval - coarse(k,j-1,i))/dx2m;
        Real gx2p = (coarse(k,j+1,i) - ccval)/dx2p;
        Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*std::min(std::abs(gx2m),
                                                           std::abs(gx2p));

        fine(fk,fj  ,fi) = ccval - gx2c*(x2c - fx2m);
        fine(fk,fj+1,fi) = ccval + gx2c*(fx2p - x2c);
      }
    }
  } else { // 1D
    for (int i=si; i<=ei; i++) {
      int fi = (i - pmb->cis)*2 + pmb->is;
      fine(0,0,fi) = coarse(0,0,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateSharedFieldX2(const AthenaArray<Real> &coarse,
//!     AthenaArray<Real> &fine, int si, int ei, int sj, int ej, int sk, int ek)
//! \brief prolongate x2 face-centered fields shared between coarse and fine levels

void MeshRefinement::ProlongateSharedFieldX2(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int fsi = (si - pmb->cis)*2 + pmb->is, fei = (ei - pmb->cis)*2 + pmb->is + 1;
  if (pmb->block_size.nx3 > 1) {
    if (fluxinterp_) {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        const Real& x3m = pcoarsec->x3s2(k-1);
        const Real& x3c = pcoarsec->x3s2(k);
        const Real& x3p = pcoarsec->x3s2(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        const Real& fx3m = pco->x3s2(fk);
        const Real& fx3p = pco->x3s2(fk+1);
        Real dfx3 = fx3p - fx3m;
        for (int j=sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          pco->Face2Area(fk,   fj,  fsi, fei, sarea_x2_[0][0]);
          pco->Face2Area(fk+1, fj,  fsi, fei, sarea_x2_[1][0]);
          pcoarsec->Face2Area(k, j, si, ei, csarea_x2_);
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            const Real& x1m = pcoarsec->x1s2(i-1);
            const Real& x1c = pcoarsec->x1s2(i);
            const Real& x1p = pcoarsec->x1s2(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            const Real& fx1m = pco->x1s2(fi);
            const Real& fx1p = pco->x1s2(fi+1);
            Real ccval = coarse(k,j,i);
            Real dfx1 = fx1p - fx1m;
            Real csarea = csarea_x2_(i) + TINY_NUMBER;
            Real fsa00 = sarea_x2_[0][0](fi);
            Real fsa01 = sarea_x2_[0][0](fi+1);
            Real fsa10 = sarea_x2_[1][0](fi);
            Real fsa11 = sarea_x2_[1][0](fi+1);

            Real gx1m = (ccval - coarse(k,j,i-1))/dx1m;
            Real gx1p = (coarse(k,j,i+1) - ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                             std::abs(gx1p));
            Real gx3m = (ccval - coarse(k-1,j,i))/dx3m;
            Real gx3p = (coarse(k+1,j,i) - ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*std::min(std::abs(gx3m),
                                                             std::abs(gx3p));

            fine(fk  ,fj  ,fi) = ccval - gx1c*dfx1*(fsa01 + fsa11)/csarea
                                       - gx3c*dfx3*(fsa10 + fsa11)/csarea;
            fine(fk  ,fj,fi+1) = ccval + gx1c*dfx1*(fsa00 + fsa10)/csarea
                                       - gx3c*dfx3*(fsa10 + fsa11)/csarea;
            fine(fk+1,fj  ,fi) = ccval - gx1c*dfx1*(fsa01 + fsa11)/csarea
                                       + gx3c*dfx3*(fsa00 + fsa01)/csarea;
            fine(fk+1,fj,fi+1) = ccval + gx1c*dfx1*(fsa00 + fsa10)/csarea
                                       + gx3c*dfx3*(fsa00 + fsa01)/csarea;
          }
        }
      }
    } else {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        const Real& x3m = pcoarsec->x3s2(k-1);
        const Real& x3c = pcoarsec->x3s2(k);
        const Real& x3p = pcoarsec->x3s2(k+1);
        Real dx3m = x3c - x3m;
        Real dx3p = x3p - x3c;
        const Real& fx3m = pco->x3s2(fk);
        const Real& fx3p = pco->x3s2(fk+1);
        for (int j=sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            const Real& x1m = pcoarsec->x1s2(i-1);
            const Real& x1c = pcoarsec->x1s2(i);
            const Real& x1p = pcoarsec->x1s2(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            const Real& fx1m = pco->x1s2(fi);
            const Real& fx1p = pco->x1s2(fi+1);
            Real ccval = coarse(k,j,i);

            Real gx1m = (ccval - coarse(k,j,i-1))/dx1m;
            Real gx1p = (coarse(k,j,i+1) - ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                               std::abs(gx1p));
            Real gx3m = (ccval - coarse(k-1,j,i))/dx3m;
            Real gx3p = (coarse(k+1,j,i) - ccval)/dx3p;
            Real gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*std::min(std::abs(gx3m),
                                                               std::abs(gx3p));
            fine(fk  ,fj,fi  ) = ccval - gx1c*(x1c - fx1m) - gx3c*(x3c - fx3m);
            fine(fk  ,fj,fi+1) = ccval + gx1c*(fx1p - x1c) - gx3c*(x3c - fx3m);
            fine(fk+1,fj,fi  ) = ccval - gx1c*(x1c - fx1m) + gx3c*(fx3p - x3c);
            fine(fk+1,fj,fi+1) = ccval + gx1c*(fx1p - x1c) + gx3c*(fx3p - x3c);
          }
        }
      }
    }
  } else if (pmb->block_size.nx2 > 1) {
    int k = pmb->cks, fk = pmb->ks;
    for (int j=sj; j<=ej; j++) {
      int fj = (j - pmb->cjs)*2 + pmb->js;
      for (int i=si; i<=ei; i++) {
        int fi = (i - pmb->cis)*2 + pmb->is;
        const Real& x1m = pcoarsec->x1s2(i-1);
        const Real& x1c = pcoarsec->x1s2(i);
        const Real& x1p = pcoarsec->x1s2(i+1);
        const Real& fx1m = pco->x1s2(fi);
        const Real& fx1p = pco->x1s2(fi+1);
        Real ccval = coarse(k,j,i);

        Real gx1m = (ccval - coarse(k,j,i-1))/(x1c - x1m);
        Real gx1p = (coarse(k,j,i+1) - ccval)/(x1p - x1c);
        Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                           std::abs(gx1p));

        fine(fk,fj,fi  ) = ccval - gx1c*(x1c - fx1m);
        fine(fk,fj,fi+1) = ccval + gx1c*(fx1p - x1c);
      }
    }
  } else {
    for (int i=si; i<=ei; i++) {
      int fi = (i - pmb->cis)*2 + pmb->is;
      Real gxm = (coarse(0,0,i) - coarse(0,0,i-1))
                 /(pcoarsec->x1s2(i) - pcoarsec->x1s2(i-1));
      Real gxp = (coarse(0,0,i+1) - coarse(0,0,i))
                 /(pcoarsec->x1s2(i+1) - pcoarsec->x1s2(i));
      Real gxc = 0.5*(SIGN(gxm) + SIGN(gxp))*std::min(std::abs(gxm),
                                                      std::abs(gxp));
      fine(0,0,fi  ) = fine(0,1,fi  )
                     = coarse(0,0,i) - gxc*(pcoarsec->x1s2(i) - pco->x1s2(fi));
      fine(0,0,fi+1) = fine(0,1,fi+1)
                     = coarse(0,0,i) + gxc*(pco->x1s2(fi+1) - pcoarsec->x1s2(i));
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateSharedFieldX3(const AthenaArray<Real> &coarse,
//!     AthenaArray<Real> &fine, int si, int ei, int sj, int ej, int sk, int ek)
//! \brief prolongate x3 face-centered fields shared between coarse and fine levels

void MeshRefinement::ProlongateSharedFieldX3(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int fsi = (si - pmb->cis)*2 + pmb->is, fei = (ei - pmb->cis)*2 + pmb->is + 1;
  if (pmb->block_size.nx3 > 1) {
    if (fluxinterp_) {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        for (int j=sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          const Real& x2m = pcoarsec->x2s3(j-1);
          const Real& x2c = pcoarsec->x2s3(j);
          const Real& x2p = pcoarsec->x2s3(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          const Real& fx2m = pco->x2s3(fj);
          const Real& fx2p = pco->x2s3(fj+1);
          Real dfx2 = fx2p - fx2m;
          pco->Face3Area(fk,   fj,  fsi, fei, sarea_x3_[0][0]);
          pco->Face3Area(fk, fj+1,  fsi, fei, sarea_x3_[0][1]);
          pcoarsec->Face3Area(k, j, si, ei, csarea_x3_);

          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            const Real& x1m = pcoarsec->x1s3(i-1);
            const Real& x1c = pcoarsec->x1s3(i);
            const Real& x1p = pcoarsec->x1s3(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            const Real& fx1m = pco->x1s3(fi);
            const Real& fx1p = pco->x1s3(fi+1);
            Real ccval = coarse(k,j,i);

            Real gx1m = (ccval - coarse(k,j,i-1))/dx1m;
            Real gx1p = (coarse(k,j,i+1) - ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                             std::abs(gx1p));
            Real gx2m = (ccval - coarse(k,j-1,i))/dx2m;
            Real gx2p = (coarse(k,j+1,i) - ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*std::min(std::abs(gx2m),
                                                             std::abs(gx2p));
            Real dfx1 = fx1p - fx1m;
            Real csarea = csarea_x3_(i);
            Real fsa00 = sarea_x3_[0][0](fi);
            Real fsa01 = sarea_x3_[0][0](fi+1);
            Real fsa10 = sarea_x3_[0][1](fi);
            Real fsa11 = sarea_x3_[0][1](fi+1);
            fine(fk  ,fj  ,fi) = ccval - gx1c*dfx1*(fsa01 + fsa11)/csarea
                                       - gx2c*dfx2*(fsa10 + fsa11)/csarea;
            fine(fk  ,fj,fi+1) = ccval + gx1c*dfx1*(fsa00 + fsa10)/csarea
                                       - gx2c*dfx2*(fsa10 + fsa11)/csarea;
            fine(fk,fj+1  ,fi) = ccval - gx1c*dfx1*(fsa01 + fsa11)/csarea
                                       + gx2c*dfx2*(fsa00 + fsa01)/csarea;
            fine(fk,fj+1,fi+1) = ccval + gx1c*dfx1*(fsa00 + fsa10)/csarea
                                       + gx2c*dfx2*(fsa00 + fsa01)/csarea;
          }
        }
      }
    } else {
      for (int k=sk; k<=ek; k++) {
        int fk = (k - pmb->cks)*2 + pmb->ks;
        for (int j=sj; j<=ej; j++) {
          int fj = (j - pmb->cjs)*2 + pmb->js;
          const Real& x2m = pcoarsec->x2s3(j-1);
          const Real& x2c = pcoarsec->x2s3(j);
          const Real& x2p = pcoarsec->x2s3(j+1);
          Real dx2m = x2c - x2m;
          Real dx2p = x2p - x2c;
          const Real& fx2m = pco->x2s3(fj);
          const Real& fx2p = pco->x2s3(fj+1);
          for (int i=si; i<=ei; i++) {
            int fi = (i - pmb->cis)*2 + pmb->is;
            const Real& x1m = pcoarsec->x1s3(i-1);
            const Real& x1c = pcoarsec->x1s3(i);
            const Real& x1p = pcoarsec->x1s3(i+1);
            Real dx1m = x1c - x1m;
            Real dx1p = x1p - x1c;
            const Real& fx1m = pco->x1s3(fi);
            const Real& fx1p = pco->x1s3(fi+1);
            Real ccval = coarse(k,j,i);

            Real gx1m = (ccval - coarse(k,j,i-1))/dx1m;
            Real gx1p = (coarse(k,j,i+1) - ccval)/dx1p;
            Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                             std::abs(gx1p));
            Real gx2m = (ccval - coarse(k,j-1,i))/dx2m;
            Real gx2p = (coarse(k,j+1,i) - ccval)/dx2p;
            Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*std::min(std::abs(gx2m),
                                                             std::abs(gx2p));
            fine(fk,fj  ,fi  ) = ccval - gx1c*(x1c - fx1m) - gx2c*(x2c - fx2m);
            fine(fk,fj  ,fi+1) = ccval + gx1c*(fx1p - x1c) - gx2c*(x2c - fx2m);
            fine(fk,fj+1,fi  ) = ccval - gx1c*(x1c - fx1m) + gx2c*(fx2p - x2c);
            fine(fk,fj+1,fi+1) = ccval + gx1c*(fx1p - x1c) + gx2c*(fx2p - x2c);
          }
        }
      }
    }
  } else if (pmb->block_size.nx2 > 1) {
    int k = pmb->cks, fk = pmb->ks;
    for (int j=sj; j<=ej; j++) {
      int fj = (j - pmb->cjs)*2 + pmb->js;
      const Real& x2m = pcoarsec->x2s3(j-1);
      const Real& x2c = pcoarsec->x2s3(j);
      const Real& x2p = pcoarsec->x2s3(j+1);
      Real dx2m = x2c - x2m;
      Real dx2p = x2p - x2c;
      const Real& fx2m = pco->x2s3(fj);
      const Real& fx2p = pco->x2s3(fj+1);
      Real dx2fm = x2c - fx2m;
      Real dx2fp = fx2p - x2c;
      for (int i=si; i<=ei; i++) {
        int fi = (i - pmb->cis)*2 + pmb->is;
        const Real& x1m = pcoarsec->x1s3(i-1);
        const Real& x1c = pcoarsec->x1s3(i);
        const Real& x1p = pcoarsec->x1s3(i+1);
        Real dx1m = x1c - x1m;
        Real dx1p = x1p - x1c;
        const Real& fx1m = pco->x1s3(fi);
        const Real& fx1p = pco->x1s3(fi+1);
        Real dx1fm = x1c - fx1m;
        Real dx1fp = fx1p - x1c;
        Real ccval = coarse(k,j,i);

        // calculate 2D gradients using the minmod limiter
        Real gx1m = (ccval - coarse(k,j,i-1))/dx1m;
        Real gx1p = (coarse(k,j,i+1) - ccval)/dx1p;
        Real gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                           std::abs(gx1p));
        Real gx2m = (ccval - coarse(k,j-1,i))/dx2m;
        Real gx2p = (coarse(k,j+1,i) - ccval)/dx2p;
        Real gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*std::min(std::abs(gx2m),
                                                           std::abs(gx2p));

        // interpolate on to the finer grid
        fine(fk,fj  ,fi  ) = fine(fk+1,fj  ,fi  ) = ccval - gx1c*dx1fm-gx2c*dx2fm;
        fine(fk,fj  ,fi+1) = fine(fk+1,fj  ,fi+1) = ccval + gx1c*dx1fp-gx2c*dx2fm;
        fine(fk,fj+1,fi  ) = fine(fk+1,fj+1,fi  ) = ccval - gx1c*dx1fm+gx2c*dx2fp;
        fine(fk,fj+1,fi+1) = fine(fk+1,fj+1,fi+1) = ccval + gx1c*dx1fp+gx2c*dx2fp;
      }
    }
  } else {
    for (int i=si; i<=ei; i++) {
      int fi = (i - pmb->cis)*2 + pmb->is;
      Real gxm = (coarse(0,0,i)   - coarse(0,0,i-1))
                 / (pcoarsec->x1s3(i) - pcoarsec->x1s3(i-1));
      Real gxp = (coarse(0,0,i+1) - coarse(0,0,i))
                 / (pcoarsec->x1s3(i+1) - pcoarsec->x1s3(i));
      Real gxc = 0.5*(SIGN(gxm) + SIGN(gxp))*std::min(std::abs(gxm),
                                                      std::abs(gxp));
      fine(0,0,fi  ) = fine(1,0,fi  )
                     = coarse(0,0,i) - gxc*(pcoarsec->x1s3(i) - pco->x1s3(fi));
      fine(0,0,fi+1) = fine(1,0,fi+1)
                     = coarse(0,0,i) + gxc*(pco->x1s3(fi+1) - pcoarsec->x1s3(i));
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateInternalField(FaceField &fine,
//!                          int si, int ei, int sj, int ej, int sk, int ek)
//! \brief prolongate the internal face-centered fields

void MeshRefinement::ProlongateInternalField(
    FaceField &fine, int si, int ei, int sj, int ej, int sk, int ek) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int fsi = (si - pmb->cis)*2 + pmb->is, fei = (ei - pmb->cis)*2 + pmb->is + 1;
  if (pmb->block_size.nx3 > 1) {
    for (int k=sk; k<=ek; k++) {
      int fk = (k - pmb->cks)*2 + pmb->ks;
      for (int j=sj; j<=ej; j++) {
        int fj = (j - pmb->cjs)*2 + pmb->js;
        pco->Face1Area(fk,   fj,   fsi, fei+1, sarea_x1_[0][0]);
        pco->Face1Area(fk,   fj+1, fsi, fei+1, sarea_x1_[0][1]);
        pco->Face1Area(fk+1, fj,   fsi, fei+1, sarea_x1_[1][0]);
        pco->Face1Area(fk+1, fj+1, fsi, fei+1, sarea_x1_[1][1]);
        pco->Face2Area(fk,   fj,   fsi, fei,   sarea_x2_[0][0]);
        pco->Face2Area(fk,   fj+1, fsi, fei,   sarea_x2_[0][1]);
        pco->Face2Area(fk,   fj+2, fsi, fei,   sarea_x2_[0][2]);
        pco->Face2Area(fk+1, fj,   fsi, fei,   sarea_x2_[1][0]);
        pco->Face2Area(fk+1, fj+1, fsi, fei,   sarea_x2_[1][1]);
        pco->Face2Area(fk+1, fj+2, fsi, fei,   sarea_x2_[1][2]);
        pco->Face3Area(fk,   fj,   fsi, fei,   sarea_x3_[0][0]);
        pco->Face3Area(fk,   fj+1, fsi, fei,   sarea_x3_[0][1]);
        pco->Face3Area(fk+1, fj,   fsi, fei,   sarea_x3_[1][0]);
        pco->Face3Area(fk+1, fj+1, fsi, fei,   sarea_x3_[1][1]);
        pco->Face3Area(fk+2, fj,   fsi, fei,   sarea_x3_[2][0]);
        pco->Face3Area(fk+2, fj+1, fsi, fei,   sarea_x3_[2][1]);
        for (int i=si; i<=ei; i++) {
          int fi = (i - pmb->cis)*2 + pmb->is;
          Real Uxx = 0.0, Vyy = 0.0, Wzz = 0.0;
          Real Uxyz = 0.0, Vxyz = 0.0, Wxyz = 0.0;
#pragma unroll
          for (int jj=0; jj<2; jj++) {
            int js = 2*jj - 1, fjj = fj + jj, fjp = fj + 2*jj;
#pragma unroll
            for (int ii=0; ii<2; ii++) {
              int is = 2*ii - 1, fii = fi + ii, fip = fi + 2*ii;
              Uxx += is*(js*(fine.x2f(fk  ,fjp,fii)*sarea_x2_[0][2*jj](fii) +
                             fine.x2f(fk+1,fjp,fii)*sarea_x2_[1][2*jj](fii))
                         +(fine.x3f(fk+2,fjj,fii)*sarea_x3_[2][  jj](fii) -
                           fine.x3f(fk  ,fjj,fii)*sarea_x3_[0][  jj](fii)));
              Vyy += js*(   (fine.x3f(fk+2,fjj,fii)*sarea_x3_[2][  jj](fii) -
                             fine.x3f(fk  ,fjj,fii)*sarea_x3_[0][  jj](fii))
                            +is*(fine.x1f(fk  ,fjj,fip)*sarea_x1_[0][  jj](fip) +
                                 fine.x1f(fk+1,fjj,fip)*sarea_x1_[1][  jj](fip)));
              Wzz +=     is*(fine.x1f(fk+1,fjj,fip)*sarea_x1_[1][  jj](fip) -
                             fine.x1f(fk  ,fjj,fip)*sarea_x1_[0][  jj](fip))
                         +js*(fine.x2f(fk+1,fjp,fii)*sarea_x2_[1][2*jj](fii) -
                              fine.x2f(fk  ,fjp,fii)*sarea_x2_[0][2*jj](fii));
              Uxyz += is*js*(fine.x1f(fk+1,fjj,fip)*sarea_x1_[1][  jj](fip) -
                             fine.x1f(fk  ,fjj,fip)*sarea_x1_[0][  jj](fip));
              Vxyz += is*js*(fine.x2f(fk+1,fjp,fii)*sarea_x2_[1][2*jj](fii) -
                             fine.x2f(fk  ,fjp,fii)*sarea_x2_[0][2*jj](fii));
              Wxyz += is*js*(fine.x3f(fk+2,fjj,fii)*sarea_x3_[2][  jj](fii) -
                             fine.x3f(fk  ,fjj,fii)*sarea_x3_[0][  jj](fii));
            }
          }
          Real Sdx1 = SQR(pco->dx1f(fi) + pco->dx1f(fi+1));
          Real Sdx2 = SQR(pco->GetEdge2Length(fk+1,fj,fi+1) +
                          pco->GetEdge2Length(fk+1,fj+1,fi+1));
          Real Sdx3 = SQR(pco->GetEdge3Length(fk,fj+1,fi+1) +
                          pco->GetEdge3Length(fk+1,fj+1,fi+1));
          Uxx *= 0.125; Vyy *= 0.125; Wzz *= 0.125;
          Uxyz *= 0.125/(Sdx2 + Sdx3);
          Vxyz *= 0.125/(Sdx1 + Sdx3);
          Wxyz *= 0.125/(Sdx1 + Sdx2);
          fine.x1f(fk  ,fj  ,fi+1) =
              (0.5*(fine.x1f(fk  ,fj  ,fi  )*sarea_x1_[0][0](fi  ) +
                    fine.x1f(fk  ,fj  ,fi+2)*sarea_x1_[0][0](fi+2))
               + Uxx - Sdx3*Vxyz - Sdx2*Wxyz) /sarea_x1_[0][0](fi+1);
          fine.x1f(fk  ,fj+1,fi+1) =
              (0.5*(fine.x1f(fk  ,fj+1,fi  )*sarea_x1_[0][1](fi  ) +
                    fine.x1f(fk  ,fj+1,fi+2)*sarea_x1_[0][1](fi+2))
               + Uxx - Sdx3*Vxyz + Sdx2*Wxyz) /sarea_x1_[0][1](fi+1);
          fine.x1f(fk+1,fj  ,fi+1) =
              (0.5*(fine.x1f(fk+1,fj  ,fi  )*sarea_x1_[1][0](fi  ) +
                    fine.x1f(fk+1,fj  ,fi+2)*sarea_x1_[1][0](fi+2))
               + Uxx + Sdx3*Vxyz - Sdx2*Wxyz) /sarea_x1_[1][0](fi+1);
          fine.x1f(fk+1,fj+1,fi+1) =
              (0.5*(fine.x1f(fk+1,fj+1,fi  )*sarea_x1_[1][1](fi  ) +
                    fine.x1f(fk+1,fj+1,fi+2)*sarea_x1_[1][1](fi+2))
               + Uxx + Sdx3*Vxyz + Sdx2*Wxyz) /sarea_x1_[1][1](fi+1);

          fine.x2f(fk  ,fj+1,fi  ) =
              (0.5*(fine.x2f(fk  ,fj  ,fi  )*sarea_x2_[0][0](fi  ) +
                    fine.x2f(fk  ,fj+2,fi  )*sarea_x2_[0][2](fi  ))
               + Vyy - Sdx3*Uxyz - Sdx1*Wxyz) /sarea_x2_[0][1](fi  );
          fine.x2f(fk  ,fj+1,fi+1) =
              (0.5*(fine.x2f(fk  ,fj  ,fi+1)*sarea_x2_[0][0](fi+1) +
                    fine.x2f(fk  ,fj+2,fi+1)*sarea_x2_[0][2](fi+1))
               + Vyy - Sdx3*Uxyz + Sdx1*Wxyz) /sarea_x2_[0][1](fi+1);
          fine.x2f(fk+1,fj+1,fi  ) =
              (0.5*(fine.x2f(fk+1,fj  ,fi  )*sarea_x2_[1][0](fi  ) +
                    fine.x2f(fk+1,fj+2,fi  )*sarea_x2_[1][2](fi  ))
               + Vyy + Sdx3*Uxyz - Sdx1*Wxyz) /sarea_x2_[1][1](fi  );
          fine.x2f(fk+1,fj+1,fi+1) =
              (0.5*(fine.x2f(fk+1,fj  ,fi+1)*sarea_x2_[1][0](fi+1) +
                    fine.x2f(fk+1,fj+2,fi+1)*sarea_x2_[1][2](fi+1))
               + Vyy + Sdx3*Uxyz + Sdx1*Wxyz) /sarea_x2_[1][1](fi+1);

          fine.x3f(fk+1,fj  ,fi  ) =
              (0.5*(fine.x3f(fk+2,fj  ,fi  )*sarea_x3_[2][0](fi  ) +
                    fine.x3f(fk  ,fj  ,fi  )*sarea_x3_[0][0](fi  ))
               + Wzz - Sdx2*Uxyz - Sdx1*Vxyz) /sarea_x3_[1][0](fi  );
          fine.x3f(fk+1,fj  ,fi+1) =
              (0.5*(fine.x3f(fk+2,fj  ,fi+1)*sarea_x3_[2][0](fi+1) +
                    fine.x3f(fk  ,fj  ,fi+1)*sarea_x3_[0][0](fi+1))
               + Wzz - Sdx2*Uxyz + Sdx1*Vxyz) /sarea_x3_[1][0](fi+1);
          fine.x3f(fk+1,fj+1,fi  ) =
              (0.5*(fine.x3f(fk+2,fj+1,fi  )*sarea_x3_[2][1](fi  ) +
                    fine.x3f(fk  ,fj+1,fi  )*sarea_x3_[0][1](fi  ))
               + Wzz + Sdx2*Uxyz - Sdx1*Vxyz) /sarea_x3_[1][1](fi  );
          fine.x3f(fk+1,fj+1,fi+1) =
              (0.5*(fine.x3f(fk+2,fj+1,fi+1)*sarea_x3_[2][1](fi+1) +
                    fine.x3f(fk  ,fj+1,fi+1)*sarea_x3_[0][1](fi+1))
               + Wzz + Sdx2*Uxyz + Sdx1*Vxyz) /sarea_x3_[1][1](fi+1);
        }
      }
    }
  } else if (pmb->block_size.nx2 > 1) {
    int fk = pmb->ks;
    for (int j=sj; j<=ej; j++) {
      int fj = (j - pmb->cjs)*2 + pmb->js;
      pco->Face1Area(fk,   fj,   fsi, fei+1, sarea_x1_[0][0]);
      pco->Face1Area(fk,   fj+1, fsi, fei+1, sarea_x1_[0][1]);
      pco->Face2Area(fk,   fj,   fsi, fei,   sarea_x2_[0][0]);
      pco->Face2Area(fk,   fj+1, fsi, fei,   sarea_x2_[0][1]);
      pco->Face2Area(fk,   fj+2, fsi, fei,   sarea_x2_[0][2]);
      for (int i=si; i<=ei; i++) {
        int fi = (i - pmb->cis)*2 + pmb->is;
        Real tmp1 = 0.25*(fine.x2f(fk,fj+2,fi+1)*sarea_x2_[0][2](fi+1)
                          - fine.x2f(fk,fj,  fi+1)*sarea_x2_[0][0](fi+1)
                          - fine.x2f(fk,fj+2,fi  )*sarea_x2_[0][2](fi  )
                          + fine.x2f(fk,fj,  fi  )*sarea_x2_[0][0](fi  ));
        Real tmp2 = 0.25*(fine.x1f(fk,fj,  fi  )*sarea_x1_[0][0](fi  )
                          - fine.x1f(fk,fj,  fi+2)*sarea_x1_[0][0](fi+2)
                          - fine.x1f(fk,fj+1,fi  )*sarea_x1_[0][1](fi  )
                          + fine.x1f(fk,fj+1,fi+2)*sarea_x1_[0][1](fi+2));
        fine.x1f(fk,fj  ,fi+1) =
            (0.5*(fine.x1f(fk,fj,  fi  )*sarea_x1_[0][0](fi  )
                  +fine.x1f(fk,fj,  fi+2)*sarea_x1_[0][0](fi+2)) + tmp1)
            /sarea_x1_[0][0](fi+1);
        fine.x1f(fk,fj+1,fi+1) =
            (0.5*(fine.x1f(fk,fj+1,fi  )*sarea_x1_[0][1](fi  )
                  +fine.x1f(fk,fj+1,fi+2)*sarea_x1_[0][1](fi+2)) + tmp1)
            /sarea_x1_[0][1](fi+1);
        fine.x2f(fk,fj+1,fi  ) =
            (0.5*(fine.x2f(fk,fj,  fi  )*sarea_x2_[0][0](fi  )
                  +fine.x2f(fk,fj+2,fi  )*sarea_x2_[0][2](fi  )) + tmp2)
            /sarea_x2_[0][1](fi  );
        fine.x2f(fk,fj+1,fi+1) =
            (0.5*(fine.x2f(fk,fj,  fi+1)*sarea_x2_[0][0](fi+1)
                  +fine.x2f(fk,fj+2,fi+1)*sarea_x2_[0][2](fi+1)) + tmp2)
            /sarea_x2_[0][1](fi+1);
      }
    }
  } else {
    pco->Face1Area(0, 0, fsi, fei+1, sarea_x1_[0][0]);
    for (int i=si; i<=ei; i++) {
      int fi = (i - pmb->cis)*2 + pmb->is;
      Real ph = sarea_x1_[0][0](fi)*fine.x1f(0,0,fi);
      fine.x1f(0,0,fi+1) = ph/sarea_x1_[0][0](fi+1);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::CheckRefinementCondition()
//! \brief Check refinement criteria

void MeshRefinement::CheckRefinementCondition() {
  MeshBlock *pmb = pmy_block_;
  int ret = 0, aret = -1;
  refine_flag_ = 0;

  //! \todo **should be implemented later:**
  //! loop-over refinement criteria
  if (AMRFlag_ != nullptr)
    ret = AMRFlag_(pmb);
  aret = std::max(aret,ret);

  if (aret >= 0)
    deref_count_ = 0;
  if (aret > 0) {
    if (pmb->loc.level == pmb->pmy_mesh->max_level) {
      refine_flag_ = 0;
    } else {
      refine_flag_ = 1;
    }
  } else if (aret < 0) {
    if (pmb->loc.level == pmb->pmy_mesh->root_level) {
      refine_flag_ = 0;
      deref_count_ = 0;
    } else {
      deref_count_++;
      int ec = 0, js, je, ks, ke;
      if (pmb->block_size.nx2 > 1) {
        js = -1;
        je = 1;
      } else {
        js = 0;
        je = 0;
      }
      if (pmb->block_size.nx3 > 1) {
        ks = -1;
        ke = 1;
      } else {
        ks = 0;
        ke = 0;
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=-1; i<=1; i++)
            if (pmb->pbval->nblevel[k+1][j+1][i+1]>pmb->loc.level) ec++;
        }
      }
      if (ec > 0) {
        refine_flag_ = 0;
      } else {
        if (deref_count_ >= deref_threshold_) {
          refine_flag_ = -1;
        } else {
          refine_flag_ = 0;
        }
      }
    }
  }
  return;
}

//! \todo (felker):
//! * consider merging w/ MeshBlock::pvars_cc, etc. See meshblock.cpp

int MeshRefinement::AddToRefinement(AthenaArray<Real> *pvar_cc,
                                     AthenaArray<Real> *pcoarse_cc) {
  pvars_cc_.push_back(std::make_tuple(pvar_cc, pcoarse_cc));
  return static_cast<int>(pvars_cc_.size() - 1);
}

int MeshRefinement::AddToRefinement(FaceField *pvar_fc, FaceField *pcoarse_fc) {
  pvars_fc_.push_back(std::make_tuple(pvar_fc, pcoarse_fc));
  return static_cast<int>(pvars_fc_.size() - 1);
}

//! Currently, only called in 2x functions in bvals_refine.cpp:
//! __________
//! - BoundaryValues::RestrictGhostCellsOnSameLevel()--- to perform additional
//! restriction on primitive Hydro standard/coarse arrays (only for GR) without changing
//! the var_cc/coarse_buf pointer members of the HydroBoundaryVariable.
//!
//! - BoundaryValues::ProlongateGhostCells()--- to ensure prolongation occurs on conserved
//! (not primitive) variable standard/coarse arrays for Hydro, PassiveScalars
//!
//! Should probably consolidate this function and std::vector of tuples with
//! BoundaryVariable interface ptr members. Too much independent switching of ptrs!
//! __________
//! Even though we currently do not have special GR functionality planned for
//! PassiveScalars::coarse_r_ like Hydro::coarse_prim_
//! (it is never transferred in Mesh::LoadBalancingAndAdaptiveMeshRefinement)
//! the physical (non-periodic) boundary functions will still apply only to the PRIMITIVE
//! scalar variable arrays, thus S/AMR demand
//! 1. AthenaArray<Real> PassiveScalars::coarse_r
//! 2. ability to switch (s, coarse_s) and (r, coarse_r) ptrs in MeshRefinement::bvals_cc_

void MeshRefinement::SetHydroRefinement(HydroBoundaryQuantity hydro_type) {
  //! \todo (felker):
  //! * make more general so it can be used as SetPassiveScalarsRefinement()
  //! * e.g. refer to "int Hydro::refinement_idx" instead of assuming that
  //!   the correct tuple is in the first vector entry
  Hydro *ph = pmy_block_->phydro;
  //! hard-coded assumption that, if multilevel, then Hydro is always present
  //! and enrolled in mesh refinement in the first pvars_cc_ vector entry
  switch (hydro_type) {
    case (HydroBoundaryQuantity::cons): {
      pvars_cc_.front() = std::make_tuple(&ph->u, &ph->coarse_cons_);
      break;
    }
    case (HydroBoundaryQuantity::prim): {
      pvars_cc_.front() = std::make_tuple(&ph->w, &ph->coarse_prim_);
      break;
    }
  }
  return;
}
