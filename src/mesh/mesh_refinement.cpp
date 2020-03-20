//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh_refinement.cpp
//  \brief implements functions for static/adaptive mesh refinement

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

// BD: debug
#include "../lagrange_interp.hpp"
//-

//----------------------------------------------------------------------------------------
//! \fn MeshRefinement::MeshRefinement(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructor

MeshRefinement::MeshRefinement(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block_(pmb), deref_count_(0),
    deref_threshold_(pin->GetOrAddInteger("mesh", "derefine_count", 10)),
    AMRFlag_(pmb->pmy_mesh->AMRFlag_) {

  coutCyan("MeshRefinement::MeshRefinement\n");

  // Create coarse mesh object for parent grid
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    pcoarsec = new Cartesian(pmb, pin, true);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    pcoarsec = new Cylindrical(pmb, pin, true);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    pcoarsec = new SphericalPolar(pmb, pin, true);
  } else if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    pcoarsec = new Minkowski(pmb, pin, true);
  } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
    pcoarsec = new Schwarzschild(pmb, pin, true);
  } else if (std::strcmp(COORDINATE_SYSTEM, "kerr-schild") == 0) {
    pcoarsec = new KerrSchild(pmb, pin, true);
  } else if (std::strcmp(COORDINATE_SYSTEM, "gr_user") == 0) {
    pcoarsec = new GRUser(pmb, pin, true);
  }

  if (NGHOST % 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshRefinement constructor" << std::endl
        << "Selected --nghost=" << NGHOST << " is incompatible with mesh refinement.\n"
        << "Reconfigure with an even number of ghost cells " << std::endl;
    ATHENA_ERROR(msg);
  }

  int nc1 = pmb->ncells1;
  fvol_[0][0].NewAthenaArray(nc1);
  fvol_[0][1].NewAthenaArray(nc1);
  fvol_[1][0].NewAthenaArray(nc1);
  fvol_[1][1].NewAthenaArray(nc1);
  sarea_x1_[0][0].NewAthenaArray(nc1+1);
  sarea_x1_[0][1].NewAthenaArray(nc1+1);
  sarea_x1_[1][0].NewAthenaArray(nc1+1);
  sarea_x1_[1][1].NewAthenaArray(nc1+1);
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

  // KGF: probably don't need to preallocate space for pointers in these vectors
  pvars_cc_.reserve(3);
  pvars_fc_.reserve(3);
  pvars_vc_.reserve(3);
}


//----------------------------------------------------------------------------------------
//! \fn MeshRefinement::~MeshRefinement()
//  \brief destructor

MeshRefinement::~MeshRefinement() {
  coutCyan("MeshRefinement::~MeshRefinement\n");
  delete pcoarsec;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::RestrictCellCenteredValues(const AthenaArray<Real> &fine,
//                           AthenaArray<Real> &coarse, int sn, int en,
//                           int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict cell centered values

void MeshRefinement::RestrictCellCenteredValues(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse, int sn, int en,
    int csi, int cei, int csj, int cej, int csk, int cek) {

  coutCyan("MeshRefinement::RestrictCellCenteredValues\n");

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  int si = (csi - pmb->cis)*2 + pmb->is, ei = (cei - pmb->cis)*2 + pmb->is + 1;

  // BD: interpolation test
  //
  // note:
  // during boundary communication this will fail as restriction is performed on
  // nodes that are unpopulated


  AthenaArray<Real>& fine_ = const_cast<AthenaArray<Real>&>(fine);
  if (false) {
    // BD debug
    bool debug = true;
    int debug_block_id = 0;
    //---


    int j = pmb->js, cj = pmb->cjs, k = pmb->ks, ck = pmb->cks;

    // settings for interpolator
    const int interp_order = 2;
    const int interp_dim = 1;

    Real origin[1] = {pco->x1v(0)};
    Real delta[1] = {pco->dx1v(0)};
    int size[1] = {pmb->ncells1};

    for (int ci=csi; ci<=cei; ci++) {
      int i = (ci - pmb->cis)*2 + pmb->is;

      // target coord
      Real coord[1] = {pcoarsec->x1v(i)};

      // set up n-dimensional interpolator
      LagrangeInterpND<interp_order, interp_dim> * pinterp = nullptr;
      pinterp = new LagrangeInterpND<interp_order, interp_dim>(origin, delta,
                                                               size, coord);

      // test interpolation to point on first component
      AthenaArray<Real> src;
      for (int n=sn; n<=en; n++) {
        src.InitWithShallowSlice(fine_, n, 1);
        coarse(n, ck, cj, ci) = pinterp->eval(src.data());
      }

      // clean-up
      delete pinterp;

    }

    // kill on MeshBlock id
    if (debug)
      if (pmb->gid == debug_block_id) {
        coutBoldRed("\npcoarsec->x1v\n");
        pcoarsec->x1v.print_all();

        coutBoldRed("pco->x1v\n");
        pco->x1v.print_all();


        coutBoldRed("[post]fine\n");
        fine_.print_all();

        coutBoldRed("[post]coarse\n");
        coarse.print_all();

        coutBoldRed("[exact]coarse\n");
        pmb->DebugWaveMeshBlock(coarse, 0, pmb->ncc1, 0, 0, 0, 0, false, true);
        coarse.print_all();

        // Q();
      }

    return;
  }
  //---

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
    // printf("1d_restr");
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

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::RestrictFieldX1(const AthenaArray<Real> &fine
//      AthenaArray<Real> &coarse, int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the x1 field data and set them into the coarse buffer

void MeshRefinement::RestrictFieldX1(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  coutCyan("MeshRefinement::RestrictFieldX1\n");

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
//      AthenaArray<Real> &coarse, int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the x2 field data and set them into the coarse buffer

void MeshRefinement::RestrictFieldX2(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  coutCyan("MeshRefinement::RestrictFieldX2\n");

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
//      AthenaArray<Real> &coarse, int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict the x3 field data and set them into the coarse buffer

void MeshRefinement::RestrictFieldX3(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse,
    int csi, int cei, int csj, int cej, int csk, int cek) {
  coutCyan("MeshRefinement::RestrictFieldX3\n");

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
//! \fn void MeshRefinement::RestrictVertexCenteredValues(const AthenaArray<Real> &fine,
//                           AthenaArray<Real> &coarse, int sn, int en,
//                           int csi, int cei, int csj, int cej, int csk, int cek)
//  \brief restrict vertex centered values

void MeshRefinement::RestrictVertexCenteredValues(
    const AthenaArray<Real> &fine, AthenaArray<Real> &coarse, int sn, int en,
    int csi, int cei, int csj, int cej, int csk, int cek) {

  // For vertex centered values this becomes an injection / identity operator

  coutCyan("MeshRefinement::RestrictVertexCenteredValues\n");

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;


  // store the restricted data in the prolongation buffer for later use
  if (pmb->block_size.nx3 > 1) { // 3D
    for (int n=sn; n<=en; ++n){
      for (int ck=csk; ck<=cek; ck++) {
        int k = (ck - pmb->cks)*2 + pmb->ks;
        for (int cj=csj; cj<=cej; cj++) {
          int j = (cj - pmb->cjs)*2 + pmb->js;
          for (int ci=csi; ci<=cei; ci++) {
            int i = (ci - pmb->cis)*2 + pmb->is;
            coarse(n, ck, cj, ci) = fine(n, k, j, i);
          }
        }
      }
    }

  } else if (pmb->block_size.nx2 > 1) { // 2D
    for (int n=sn; n<=en; ++n){
      for (int cj=csj; cj<=cej; cj++) {
        int j = (cj - pmb->cjs)*2 + pmb->js;
        for (int ci=csi; ci<=cei; ci++) {
          int i = (ci - pmb->cis)*2 + pmb->is;
          coarse(n, 0, cj, ci) = fine(n, 0, j, i);
        }
      }
    }

    // debug
    // for (int cj=csj; cj<=cej; cj++) {
    //   int j = (cj - pmb->cjs)*2 + pmb->js;
    //   for (int ci=csi; ci<=cei; ci++) {
    //     int i = (ci - pmb->cis)*2 + pmb->is;
    //     printf("(ci, cj, i, j)=(%d, %d, %d, %d)\n", ci, cj, i, j);

    //     printf("fine(0, 0, j:%d, i:%d)=%1.3f\n", j, i, fine(0, 0, j, i));
    //   }
    // }
    //-

  } else { // 1D
    int j = pmb->js, cj = pmb->cjs, k = pmb->ks, ck = pmb->cks;
    for (int n=sn; n<=en; ++n) {
      for (int ci=csi; ci<=cei; ci++) {
        int i = (ci - pmb->cis)*2 + pmb->is;
        coarse(n, ck, cj, ci) = fine(n, k, j, i);
      }
    }


    // debug
    // for (int ci=csi; ci<=cei; ci++) {
    //   int i = (ci - pmb->cis)*2 + pmb->is;
    //   printf("(ci, i)=(%d, %d)\n", ci, i);
    //   printf("fine(0, k:%d, j:%d, i:%d)=%1.3f\n", k, j, i, fine(0, k, j, i));
    // }
    //-
  }



}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateCellCenteredValues(
//        const AthenaArray<Real> &coarse,AthenaArray<Real> &fine, int sn, int en,,
//        int si, int ei, int sj, int ej, int sk, int ek)
//  \brief Prolongate cell centered values

void MeshRefinement::ProlongateCellCenteredValues(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int sn, int en, int si, int ei, int sj, int ej, int sk, int ek) {
  coutCyan("MeshRefinement::ProlongateCellCenteredValues\n");

  coutBoldRed("[pre]coarse\n");
  coarse.print_all();

  coutBoldRed("[pre]fine\n");
  fine.print_all();

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;

  // BD debug: re-populate coarse grid with exact solution
  // need mutable
  AthenaArray<Real>& coarse_ = const_cast<AthenaArray<Real>&>(coarse);

  if (FILL_WAVE_COARSE_P) {
    pmb->DebugWaveMeshBlock(coarse_,
                            pmb->cims, pmb->cipe,
                            pmb->cjme, pmb->cjpe,
                            pmb->ckme, pmb->ckpe,
                            false, true);

    coutBoldRed("Warning: coarse buffer overridden..\n");
  }


  // BD: interpolation test
  if (false) {
    int k = pmb->cks, fk = pmb->ks, j = pmb->cjs, fj = pmb->js;

    // settings for interpolator
    const int interp_order = 2;
    const int interp_dim = 1;

    Real origin[1] = {pcoarsec->x1v(0)};
    Real delta[1] = {pcoarsec->dx1v(0)};
    int size[1] = {pmb->ncc1};

    for (int i=si; i<=ei; i++) {
      for (int ii=0; ii<=1; ii++) {
        // fine idx
        int fi = (i - pmb->cis)*2 + pmb->is + ii;

        // target coord
        Real coord[1] = {pco->x1v(fi)};

        // set up n-dimensional interpolator
        LagrangeInterpND<interp_order, interp_dim> * pinterp = nullptr;
        pinterp = new LagrangeInterpND<interp_order, interp_dim>(origin, delta,
                                                                 size, coord);


        // test interpolation to point on first component
        AthenaArray<Real> src;
        for (int n=sn; n<=en; n++) {
          src.InitWithShallowSlice(coarse_, n, 1);
          fine(n, fk, fj, fi) = pinterp->eval(src.data());
        }

        // clean-up
        delete pinterp;

        // BD debug: overwrite with soln
        // pmb->DebugWaveMeshBlock(fine, fi, fi, fj, fj, fk, fk, false);
        //---
      }
    }

    return;
  }
  //---


  // BD: debug- disable slope-limiter for the moment
  // this brings behaviour into line with 'master'
  bool slope_limit = false;

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
            Real ccval = coarse_(n,k,j,i);

            // calculate 3D gradients using the minmod limiter
            Real gx1c, gx2c, gx3c;
            if (slope_limit) {

              Real gx1m = (ccval - coarse_(n,k,j,i-1))/dx1m;
              Real gx1p = (coarse_(n,k,j,i+1) - ccval)/dx1p;
              gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*
                std::min(std::abs(gx1m), std::abs(gx1p));

              Real gx2m = (ccval - coarse_(n,k,j-1,i))/dx2m;
              Real gx2p = (coarse_(n,k,j+1,i) - ccval)/dx2p;
              gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*
                std::min(std::abs(gx2m), std::abs(gx2p));

              Real gx3m = (ccval - coarse_(n,k-1,j,i))/dx3m;
              Real gx3p = (coarse_(n,k+1,j,i) - ccval)/dx3p;
              gx3c = 0.5*(SIGN(gx3m) + SIGN(gx3p))*
                std::min(std::abs(gx3m), std::abs(gx3p));
            } else { // use 2nd ordered centred
                gx1c = (coarse_(n,k,j,i+1) - coarse_(n,k,j,i-1))/(2*dx1m);
                gx2c = (coarse_(n,k,j+1,i) - coarse_(n,k,j-1,i))/(2*dx2m);
                gx3c = (coarse_(n,k+1,j,i) - coarse_(n,k-1,j,i))/(2*dx3m);
            }
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
          Real ccval = coarse_(n,k,j,i);

          Real gx1c, gx2c;
          // calculate 2D gradients using the minmod limiter
          if (slope_limit) {
            Real gx1m = (ccval - coarse_(n,k,j,i-1))/dx1m;
            Real gx1p = (coarse_(n,k,j,i+1) - ccval)/dx1p;
            gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*
              std::min(std::abs(gx1m), std::abs(gx1p));
            Real gx2m = (ccval - coarse_(n,k,j-1,i))/dx2m;
            Real gx2p = (coarse_(n,k,j+1,i) - ccval)/dx2p;
            gx2c = 0.5*(SIGN(gx2m) + SIGN(gx2p))*
              std::min(std::abs(gx2m), std::abs(gx2p));
          } else { // use 2nd ordered centered
              gx1c = (coarse_(n,k,j,i+1) - coarse_(n,k,j,i-1))/(2*dx1m);
              gx2c = (coarse_(n,k,j+1,i) - coarse_(n,k,j+1,i))/(2*dx2m);
          }
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
        Real ccval = coarse_(n,k,j,i);

        Real gx1c;

        // calculate 1D gradient using the min-mod limiter
        if (slope_limit) {
          Real gx1m = (ccval - coarse_(n,k,j,i-1))/dx1m;
          Real gx1p = (coarse_(n,k,j,i+1) - ccval)/dx1p;
          gx1c = 0.5*(SIGN(gx1m) + SIGN(gx1p))*std::min(std::abs(gx1m),
                                                        std::abs(gx1p));
        } else { // use 2nd order centered
          gx1c = (coarse_(n,k,j,i+1) - coarse_(n,k,j,i-1))/(2*dx1m);
        }
        // interpolate on to the finer grid
        fine(n,fk  ,fj  ,fi  ) = ccval - gx1c*dx1fm;
        fine(n,fk  ,fj  ,fi+1) = ccval + gx1c*dx1fp;
      }
    }
  }

  for (int i=si; i<=ei; i++) {
    int fi = (i - pmb->cis)*2 + pmb->is;
    printf("(i, fi) = (%d, %d)\n", i, fi);
  }
  printf("\n");

  coutBoldRed("[post]coarse\n");
  coarse_.print_all();

  coutBoldRed("[post]fine\n");
  fine.print_all();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateSharedFieldX1(const AthenaArray<Real> &coarse,
//      AthenaArray<Real> &fine, int si, int ei, int sj, int ej, int sk, int ek)
//  \brief prolongate x1 face-centered fields shared between coarse and fine levels

void MeshRefinement::ProlongateSharedFieldX1(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int si, int ei, int sj, int ej, int sk, int ek) {
  coutCyan("MeshRefinement::ProlongateSharedFieldX1\n");

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  if (pmb->block_size.nx3 > 1) {
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
//      AthenaArray<Real> &fine, int si, int ei, int sj, int ej, int sk, int ek)
//  \brief prolongate x2 face-centered fields shared between coarse and fine levels

void MeshRefinement::ProlongateSharedFieldX2(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int si, int ei, int sj, int ej, int sk, int ek) {
  coutCyan("MeshRefinement::ProlongateSharedFieldX2\n");

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  if (pmb->block_size.nx3 > 1) {
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
//      AthenaArray<Real> &fine, int si, int ei, int sj, int ej, int sk, int ek)
//  \brief prolongate x3 face-centered fields shared between coarse and fine levels

void MeshRefinement::ProlongateSharedFieldX3(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int si, int ei, int sj, int ej, int sk, int ek) {
  coutCyan("MeshRefinement::ProlongateSharedFieldX3\n");

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  if (pmb->block_size.nx3 > 1) {
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
//                           int si, int ei, int sj, int ej, int sk, int ek)
//  \brief prolongate the internal face-centered fields

void MeshRefinement::ProlongateInternalField(
    FaceField &fine, int si, int ei, int sj, int ej, int sk, int ek) {
  coutCyan("MeshRefinement::ProlongateInternalField\n");

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

// Temporary for poly interpolation [mma]
// InterpolatingPolynomial[{{x0, x1, x2}, {f0, f1, f2}} // Transpose, x] /. x -> fn
// InputForm[%]
//
// TODO: Use David's branch / rewrite in Barycentric form
Real interp1_2(Real fn,
               Real x0, Real x1, Real x2,
               Real f0, Real f1, Real f2) {

  return (f0 + (fn - x0)*((-f0 + f1)/(-x0 + x1) +
   ((fn - x1)*(-((-f0 + f1)/(-x0 + x1)) + (-f1 + f2)/(-x1 + x2)))/(-x0 + x2)));
}

Real interp1_3(Real fn,
               Real x0, Real x1, Real x2, Real x3,
               Real f0, Real f1, Real f2, Real f3) {
  return (f0 + (fn - x0)*((-f0 + f1)/(-x0 + x1) +
   (fn - x1)*((-((-f0 + f1)/(-x0 + x1)) + (-f1 + f2)/(-x1 + x2))/(-x0 + x2) +
     ((fn - x2)*(-((-((-f0 + f1)/(-x0 + x1)) + (-f1 + f2)/(-x1 + x2))/(-x0 + x2)) +
        (-((-f1 + f2)/(-x1 + x2)) + (-f2 + f3)/(-x2 + x3))/(-x1 + x3)))/(-x0 + x3))));
}

Real interp1_4(Real fn,
               Real x0, Real x1, Real x2, Real x3, Real x4,
               Real f0, Real f1, Real f2, Real f3, Real f4) {

  return (f0 + (fn - x0)*((-f0 + f1)/(-x0 + x1) +
   (fn - x1)*((-((-f0 + f1)/(-x0 + x1)) + (-f1 + f2)/(-x1 + x2))/(-x0 + x2) +
     (fn - x2)*((-((-((-f0 + f1)/(-x0 + x1)) + (-f1 + f2)/(-x1 + x2))/(-x0 + x2)) +
         (-((-f1 + f2)/(-x1 + x2)) + (-f2 + f3)/(-x2 + x3))/(-x1 + x3))/(-x0 + x3) +
       ((fn - x3)*(-((-((-((-f0 + f1)/(-x0 + x1)) + (-f1 + f2)/(-x1 + x2))/(-x0 + x2)) +
             (-((-f1 + f2)/(-x1 + x2)) + (-f2 + f3)/(-x2 + x3))/(-x1 + x3))/(-x0 + x3)) +
          (-((-((-f1 + f2)/(-x1 + x2)) + (-f2 + f3)/(-x2 + x3))/(-x1 + x3)) +
            (-((-f2 + f3)/(-x2 + x3)) + (-f3 + f4)/(-x3 + x4))/(-x2 + x4))/(-x1 + x4)))/
        (-x0 + x4)))));
}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::ProlongateVertexCenteredValues(
//        const AthenaArray<Real> &coarse,AthenaArray<Real> &fine, int sn, int en,,
//        int si, int ei, int sj, int ej, int sk, int ek)
//  \brief Prolongate vertex centered values

void MeshRefinement::ProlongateVertexCenteredValues(
    const AthenaArray<Real> &coarse, AthenaArray<Real> &fine,
    int sn, int en, int si, int ei, int sj, int ej, int sk, int ek) {
  coutCyan("MeshRefinement::ProlongateVertexCenteredValues\n");

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;


  // BD debug: re-populate coarse grid with exact solution
  // need mutable
  AthenaArray<Real>& coarse_ = const_cast<AthenaArray<Real>&>(coarse);

  if (FILL_WAVE_COARSE_P) {
    pmb->DebugWaveMeshBlock(coarse_,
                            pmb->cims, pmb->cipe,
                            pmb->cjme, pmb->cjpe,
                            pmb->ckme, pmb->ckpe,
                            false, true);

    coutBoldRed("Warning: coarse buffer overridden..\n");
  }
  // pmb->pmr->RestrictVertexCenteredValues(fine, coarse_, sn, en,
  //                                        pmb->civs, pmb->cive, 0, 0, 0, 0);

  //--



  coutBoldRed("(block_size.nx1, pmb->is, pmb->cis, si, ei) = ");
  printf("(%d, %d, %d, %d, %d)\n\n",
         pmb->block_size.nx1, pmb->is, pmb->cis, si, ei);

  coutBoldRed("(block_size.nx2, pmb->js, pmb->cjs, sj, ej) = ");
  printf("(%d, %d, %d, %d, %d)\n\n",
         pmb->block_size.nx2, pmb->js, pmb->cjs, sj, ej);

  // note:
  // once done, branch with
  // pmb->block_size.nx3, nx2, n1 > 1 condition


  // 1d- hard-coded operators
  if (false) {

    int k = pmb->cks, fk = pmb->ks, j = pmb->cjs, fj = pmb->js;

    coarse_.print_all();

    // for interpolation bias switching
    int bi = pmb->cnghost + pmb->block_size.nx1 / 4; // block_size is fine no.

    printf("  bi = %d\n", bi);

    for (int n=sn; n<=en; n++) {
      for (int i=si; i<=ei; i++) {
        int fi;
        Real gx1c;

        if (i < bi) {
          // subtract coarse ghosts, double, add fine ghosts
          fi = (i - pmb->cis) * 2 + pmb->is;
          coutBoldYellow("bias to right: (i, fi)=");
          printf("(%d, %d)\n", i, fi);

          // bias towards right (note how populated)
          const Real& x1m = pcoarsec->x1f(i-1);
          const Real& x1c = pcoarsec->x1f(i);
          const Real& x1p = pcoarsec->x1f(i+1);
          // Real dx1m = x1c - x1m;
          // Real dx1p = x1p - x1c;

          const Real& fx1m = pco->x1f(fi);
          const Real& fx1p = pco->x1f(fi+1);
          // Real dx1fm = x1c - fx1m;
          // Real dx1fp = fx1p - x1c;

          // Real ccval = coarse_(n,k,j,i);

          // .. interp.
          printf("(x1m, x1c, x1p, fx1m, fx1p) ="
                 "(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)\n",
                 x1m, x1c, x1p, fx1m, fx1p);

          //gx1c = (coarse_(n,k,j,i+1) - coarse_(n,k,j,i-1))/(2*dx1m);
          //-

          // fine(n,fk  ,fj  ,fi  ) = ccval - gx1c * dx1fm;
          // fine(n,fk  ,fj  ,fi+1) = ccval + gx1c * dx1fp;
          // fine(n,fk  ,fj  ,fi  ) = fine(n,fk  ,fj  ,fi + 2 * NGHOST);
          // fine(n,fk  ,fj  ,fi+1) = fine(n,fk  ,fj  ,fi + 2 * NGHOST - 1);

          // polynomial interpolation
          Real cxm = pcoarsec->x1f(i-1);
          Real cxc = pcoarsec->x1f(i);
          Real cxp = pcoarsec->x1f(i+1);

          Real cfm = coarse_(n,k,j,i-1);
          Real cfc = coarse_(n,k,j,i);
          Real cfp = coarse_(n,k,j,i+1);

          Real fxc = pco->x1f(fi);
          Real fxp = pco->x1f(fi+1);

          // ord: 2
          // fine(n,fk  ,fj  ,fi  ) = interp1_2(fxc,
          //                                    cxm, cxc, cxp,
          //                                    cfm, cfc, cfp);
          // fine(n,fk  ,fj  ,fi+1) = interp1_2(fxp,
          //                                    cxm, cxc, cxp,
          //                                    cfm, cfc, cfp);

          // // ord: 3
          Real cxp1 = pcoarsec->x1f(i+2);
          Real cfp1 = coarse_(n,k,j,i+2);

          fine(n,fk  ,fj  ,fi  ) = interp1_3(fxc,
                                             cxm, cxc, cxp, cxp1,
                                             cfm, cfc, cfp, cfp1);
          fine(n,fk  ,fj  ,fi+1) = interp1_3(fxp,
                                             cxm, cxc, cxp, cxp1,
                                             cfm, cfc, cfp, cfp1);

          // ord: 4
          // Real cxp1 = pcoarsec->x1f(i+2);
          // Real cfp1 = coarse_(n,k,j,i+2);
          // Real cxp2 = pcoarsec->x1f(i+3);
          // Real cfp2 = coarse_(n,k,j,i+3);

          // fine(n,fk  ,fj  ,fi  ) = interp1_4(fxc,
          //                                    cxm, cxc, cxp, cxp1, cxp2,
          //                                    cfm, cfc, cfp, cfp1, cfp2);
          // fine(n,fk  ,fj  ,fi+1) = interp1_4(fxp,
          //                                    cxm, cxc, cxp, cxp1, cxp2,
          //                                    cfm, cfc, cfp, cfp1, cfp2);

          // debug with replacement
          // pmb->DebugWaveMeshBlock(fine, fi, fi, fj, fj, fk, fk, false);
          // pmb->DebugWaveMeshBlock(fine, fi+1, fi+1, fj, fj, fk, fk, false);

        } else if (i > bi) {
          // subtract coarse ghosts, double, add fine ghosts, compensate for vert.
          fi = (i - pmb->cis - 1) * 2 + pmb->is + 1;
          coutBoldYellow("bias to left: (i, fi)=");
          printf("(%d, %d)\n", i, fi);

          // bias towards left (note how populated)
          const Real& x1m = pcoarsec->x1f(i-1);
          const Real& x1c = pcoarsec->x1f(i);
          const Real& x1p = pcoarsec->x1f(i+1);
          Real dx1m = x1c - x1m;
          Real dx1p = x1p - x1c;

          const Real& fx1m = pco->x1f(fi);
          const Real& fx1p = pco->x1f(fi+1);
          Real dx1fm = x1c - fx1m;
          Real dx1fp = fx1p - x1c;

          Real ccval = coarse_(n,k,j,i);

          // .. interp.
          printf("(x1m, x1c, x1p, fx1m, fx1p) ="
                 "(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)\n",
                 x1m, x1c, x1p, fx1m, fx1p);

          // gx1c = (coarse_(n,k,j,i+1) - coarse_(n,k,j,i-1))/(2*dx1m);
          //-

          // fine(n,fk  ,fj  ,fi  ) = ccval - gx1c * dx1fm;
          // fine(n,fk  ,fj  ,fi+1) = ccval + gx1c * dx1fp;
          // fine(n,fk  ,fj  ,fi) = fine(n,fk  ,fj  ,fi - 2);
          // fine(n,fk  ,fj  ,fi+1) = fine(n,fk  ,fj  ,fi - 3);

          // polynomial interpolation
          Real cxm = pcoarsec->x1f(i-1);
          Real cxc = pcoarsec->x1f(i);
          Real cxp = pcoarsec->x1f(i+1);

          Real cfm = coarse_(n,k,j,i-1);
          Real cfc = coarse_(n,k,j,i);
          Real cfp = coarse_(n,k,j,i+1);

          Real fxc = pco->x1f(fi);
          Real fxp = pco->x1f(fi+1);

          // // ord: 2
          // fine(n,fk  ,fj  ,fi  ) = interp1_2(fxc,
          //                                    cxm, cxc, cxp,
          //                                    cfm, cfc, cfp);
          // fine(n,fk  ,fj  ,fi+1) = interp1_2(fxp,
          //                                    cxm, cxc, cxp,
          //                                    cfm, cfc, cfp);

          // // ord: 3
          Real cxm1 = pcoarsec->x1f(i-2);
          Real cfm1 = coarse_(n,k,j,i-2);
          fine(n,fk  ,fj  ,fi  ) = interp1_3(fxc,
                                             cxm1, cxm, cxc, cxp,
                                             cfm1, cfm, cfc, cfp);
          fine(n,fk  ,fj  ,fi+1) = interp1_3(fxp,
                                             cxm1, cxm, cxc, cxp,
                                             cfm1, cfm, cfc, cfp);

          // ord: 4
          // Real cxm1 = pcoarsec->x1f(i-2);
          // Real cfm1 = coarse_(n,k,j,i-2);
          // Real cxm2 = pcoarsec->x1f(i-3);
          // Real cfm2 = coarse_(n,k,j,i-3);
          // fine(n,fk  ,fj  ,fi  ) = interp1_4(fxc,
          //                                    cxm2, cxm1, cxm, cxc, cxp,
          //                                    cfm2, cfm1, cfm, cfc, cfp);
          // fine(n,fk  ,fj  ,fi+1) = interp1_4(fxp,
          //                                    cxm2, cxm1, cxm, cxc, cxp,
          //                                    cfm2, cfm1, cfm, cfc, cfp);

          // debug with replacement
          // pmb->DebugWaveMeshBlock(fine, fi, fi, fj, fj, fk, fk, false);
          // pmb->DebugWaveMeshBlock(fine, fi+1, fi+1, fj, fj, fk, fk, false);


        } else { // center node
          // central node location inferred from fine parameters
          fi = pmb->is + pmb->block_size.nx1 / 2;
          coutBoldYellow("centered: (i, fi)=");
          printf("(%d, %d)\n", i, fi);

          // centered
          coutBoldRed("ADD CENTRAL NODE (inject)\n");
          // .. interp.

          //-
          Q();
        }

      }
    }

    fine.print_all();
    // if (si > 4)
    //   Q();

    // to avoid lower fcn overwrites
    return;
  }


  // 1d- LagrangeInterpND
  if (pmb->pmy_mesh->ndim == 1) {
    // BD debug
    bool debug = false;
    int debug_block_id = 0;
    //---

    int k = pmb->cks, fk = pmb->ks, j = pmb->cjs, fj = pmb->js;

    // BD debug: print additional info
    if (debug) {
      coutBoldRed("\npcoarsec->x1f\n");
      pcoarsec->x1f.print_all();

      coutBoldRed("pco->x1f\n");
      pco->x1f.print_all();

      coutBoldRed("[pre]coarse\n");
      coarse_.print_all();

      coutBoldRed("[pre]fine\n");
      fine.print_all();
    }
    //---

    // settings for interpolator
    const int interp_order = 2 * NGHOST;
    const int interp_dim = 1;

    Real origin[1] = {pcoarsec->x1f(0)};
    Real delta[1] = {pcoarsec->dx1f(0)};
    int size[1] = {pmb->ncv1};

    for (int i=si; i<=ei; i++) {
      for (int ii=0; ii<=1; ii++) {
        // fine idx
        int fi = (i - pmb->civs) * 2 + pmb->ivs + ii;

        // target coord
        Real coord[1] = {pco->x1f(fi)};

        // set up n-dimensional interpolator
        LagrangeInterpND<interp_order, interp_dim> * pinterp = nullptr;
        pinterp = new LagrangeInterpND<interp_order, interp_dim>(origin, delta,
                                                                 size, coord);


        // test interpolation to point on first component
        AthenaArray<Real> src;
        for (int n=sn; n<=en; n++) {
          src.InitWithShallowSlice(coarse_, n, 1);
          fine(n, fk, fj, fi) = pinterp->eval(src.data());
        }

        // clean-up
        delete pinterp;

        // BD debug: overwrite with soln
        // pmb->DebugWaveMeshBlock(fine, fi, fi, fj, fj, fk, fk, false);
        //---
      }
    }

    // kill on MeshBlock id
    if (debug)
      if (pmb->gid == debug_block_id) {
        coutBoldRed("\npcoarsec->x1f\n");
        pcoarsec->x1f.print_all();

        coutBoldRed("pco->x1f\n");
        pco->x1f.print_all();


        coutBoldRed("[post]coarse\n");
        coarse_.print_all();
        coutBoldRed("[post]fine\n");
        fine.print_all();

        coutBoldRed("[exact]fine\n");
        pmb->DebugWaveMeshBlock(fine, pmb->ims, pmb->ipe, 0, 0, 0, 0, false);
        fine.print_all();

        Q();
      }

  }

  // 2d LagrangeInterpND
  if (pmb->pmy_mesh->ndim == 2) {
    // BD debug
    bool debug = true;
    int debug_block_id = 3;
    //---

    int k = pmb->cks, fk = pmb->ks;

    // BD debug: print additional info
    if (debug) {
      coutBoldRed("\npcoarsec->x1f\n");
      pcoarsec->x1f.print_all();

      coutBoldRed("\npcoarsec->x2f\n");
      pcoarsec->x2f.print_all();

      coutBoldRed("pco->x2f\n");
      pco->x2f.print_all();

      coutBoldRed("[pre]coarse\n");
      coarse_.print_all();

      coutBoldRed("[pre]fine\n");
      fine.print_all();
    }
    //---

    AthenaArray<Real>& coarse_ = const_cast<AthenaArray<Real>&>(coarse);

    // settings for interpolator
    const int interp_order = 2 * NGHOST;
    const int interp_dim = 2;

    Real origin[2] = {pcoarsec->x1f(0), pcoarsec->x2f(0)};
    Real delta[2] = {pcoarsec->dx1f(0), pcoarsec->dx2f(0)};
    int size[2] = {pmb->ncv1, pmb->ncv2};

    // for interpolation bias switching
    // int bi = pmb->cnghost + pmb->block_size.nx1 / 4; // block_size is fine no.
    // int bj = pmb->cnghost + pmb->block_size.nx2 / 4; // block_size is fine no.

    // printf("  bi, bj = %d, %d\n", bi, bj);

    for (int n=sn; n<=en; n++) {
      for (int j=sj; j<=ej; j++) {
        for (int jj=0; jj<=1; jj++) {
          // current fine index
          int fj = (j - pmb->cjvs)*2 + pmb->jvs + jj;

          for (int i=si; i<=ei; i++) {

            for (int ii=0; ii<=1; ii++) {
              // current coarse index
              int fi = (i - pmb->civs)*2 + pmb->ivs + ii;

              // target coord
              Real coord[2] = {pco->x1f(fi), pco->x2f(fj)};

              // set up n-dimensional interpolator
              LagrangeInterpND<interp_order, interp_dim> * pinterp = nullptr;
              pinterp = new LagrangeInterpND<interp_order, interp_dim>(origin, delta,
                                                                       size, coord);

              // test interpolation to point on variable first component
              AthenaArray<Real> src;
              for (int n=sn; n<=en; n++) {
                src.InitWithShallowSlice(coarse_, n, 1);
                fine(n, fk, fj, fi) = pinterp->eval(src.data());
              }

              // clean-up
              delete pinterp;

              // BD debug: overwrite with soln
              // pmb->DebugWaveMeshBlock(fine, fi, fi, fj, fj, fk, fk, false);
              // pmb->DebugWaveMeshBlock(fine, fi+1, fi+1, fj, fj, fk, fk, false);
              // pmb->DebugWaveMeshBlock(fine, fi, fi, fj+1, fj+1, fk, fk, false);
              // pmb->DebugWaveMeshBlock(fine, fi+1, fi+1, fj+1, fj+1, fk, fk, false);
              //---

            }

          }


        }
      }
    }

    // kill on MeshBlock id
    if (debug)
      if (pmb->gid == debug_block_id) {
        coutBoldRed("\npcoarsec->x1f\n");
        pcoarsec->x1f.print_all();

        coutBoldRed("pco->x1f\n");
        pco->x1f.print_all();

        coutBoldRed("\npcoarsec->x2f\n");
        pcoarsec->x2f.print_all();

        coutBoldRed("pco->x2f\n");
        pco->x2f.print_all();

        coutBoldRed("[post]coarse\n");
        coarse_.print_all();
        coutBoldRed("[post]fine\n");
        fine.print_all();

        coutBoldRed("[exact]fine\n");
        pmb->DebugWaveMeshBlock(fine, pmb->ims, pmb->ipe, pmb->jms, pmb->jpe,
                                0, 0, false);
        fine.print_all();

        Q();
      }


  }


  // BD debug: dump regardless of break
  coutBoldRed("[post]coarse\n");
  coarse_.print_all();

  coutBoldRed("[post]fine\n");
  fine.print_all();
  //---

}

//----------------------------------------------------------------------------------------
//! \fn void MeshRefinement::CheckRefinementCondition()
//  \brief Check refinement criteria

void MeshRefinement::CheckRefinementCondition() {
  coutCyan("MeshRefinement::CheckRefinementCondition\n");

  MeshBlock *pmb = pmy_block_;
  int ret = 0, aret = -1;
  refine_flag_ = 0;

  // *** should be implemented later ***
  // loop-over refinement criteria
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

// TODO(felker): consider merging w/ MeshBlock::pvars_cc, etc. See meshblock.cpp

int MeshRefinement::AddToRefinement(AthenaArray<Real> *pvar_in,
                                    AthenaArray<Real> *pcoarse_in) {
  coutCyan("MeshRefinement::AddToRefinement\n");

  // BD: TODO: one should do this differently...
  if (PREFER_VC) {
    pvars_vc_.push_back(std::make_tuple(pvar_in, pcoarse_in));
    return static_cast<int>(pvars_vc_.size() - 1);
  } else {
    pvars_cc_.push_back(std::make_tuple(pvar_in, pcoarse_in));
    return static_cast<int>(pvars_cc_.size() - 1);
  }
}

int MeshRefinement::AddToRefinement(FaceField *pvar_fc, FaceField *pcoarse_fc) {
  coutCyan("MeshRefinement::AddToRefinement\n");

  pvars_fc_.push_back(std::make_tuple(pvar_fc, pcoarse_fc));
  return static_cast<int>(pvars_fc_.size() - 1);
}

// Currently, only called in 2x functions in bvals_refine.cpp:
// ----------
// - BoundaryValues::RestrictGhostCellsOnSameLevel()--- to perform additional
// restriction on primitive Hydro standard/coarse arrays (only for GR) without changing
// the var_cc/coarse_buf pointer members of the HydroBoundaryVariable.

// - BoundaryValues::ProlongateGhostCells()--- to ensure prolongation occurs on conserved
// (not primitive) variable standard/coarse arrays for Hydro, PassiveScalars

// Should probably consolidate this function and std::vector of tuples with
// BoundaryVariable interface ptr members. Too much independent switching of ptrs!
// ----------
// Even though we currently do not have special GR functionality planned for
// PassiveScalars::coarse_r_ like Hydro::coarse_prim_
// (it is never transferred in Mesh::LoadBalancingAndAdaptiveMeshRefinement)
// the physical (non-periodic) boundary functions will still apply only to the PRIMITIVE
// scalar variable arrays, thus S/AMR demand 1) AthenaArray<Real> PassiveScalars::coarse_r
// 2) ability to switch (s, coarse_s) and (r, coarse_r) ptrs in MeshRefinement::bvals_cc_

void MeshRefinement::SetHydroRefinement(HydroBoundaryQuantity hydro_type) {
  coutCyan("MeshRefinement::SetHydroRefinement\n");

  // TODO(felker): make more general so it can be used as SetPassiveScalarsRefinement()
  // e.g. refer to "int Hydro::refinement_idx" instead of assuming that the correct tuple
  // is in the first vector entry
  Hydro *ph = pmy_block_->phydro;
  // hard-coded assumption that, if multilevel, then Hydro is always present
  // and enrolled in mesh refinement in the first pvars_cc_ vector entry
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
