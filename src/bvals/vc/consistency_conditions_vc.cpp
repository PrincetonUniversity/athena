//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file consistency_conditions_vc.cpp
//  \brief additive unpack of vertex-centered data needs averaging

// C headers

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "bvals_vc.hpp"


//////////////////////////////////
// debug functions
void print_int_arr3(int *arr, int dim3=3, int dim2=3, int dim1=3){
  for (int k=0; k<dim3; k++) {
    printf("\nk=%d:\n", k); // slab for readability
    for (int j=0; j<dim2; j++) {
      for (int i=0; i<dim1; i++) {
        printf("%d, ", arr[i + dim1 * (j + dim2 * (k))]);
      }
      printf("\n");
    }
  }

}


void VertexCenteredBoundaryVariable::ErrorUnknownMultiplicity() {
  std::stringstream msg;
  msg << "### FATAL ERROR" << std::endl
      << "Node multiplicty for vertex-centered unknown." << std::endl;
  ATHENA_ERROR(msg);
  return;

}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ZeroVertexGhosts()
//  \brief zero out ghost zones

void VertexCenteredBoundaryVariable::ZeroVertexGhosts() {
  coutYellow("VertexCenteredBoundaryVariable::ZeroVertexGhosts\n");
  // additive unpack used to populate ghosts entails that old ghost data needs
  // to be cleaned

  // BD: TODO coarse_buf needs to be zerod also...
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  if (pmb->block_size.nx3 > 1) {

    // fill bands
    for (int n=nl_; n<=nu_; ++n) {
      for (int k=0; k<NGHOST; k++)
        for (int j=0; j<=pmb->jpe; j++)
#pragma omp simd
          for (int i=0; i<=pmb->ipe; i++) {
            var(n, k, j, i) = 0;
            var(n, pmb->kps + k, j, i) = 0;
          }

      for (int k=pmb->kvs; k<=pmb->kve; k++)
        for (int j=0; j<=pmb->jpe; j++)
#pragma omp simd
          for (int i=0; i<NGHOST; i++) {
            var(n, k, j, i) = 0;
            var(n, k, j, pmb->ips + i) = 0;
          }

      for (int k=pmb->kvs; k<=pmb->kve; k++)
        for (int j=0; j<NGHOST; j++)
#pragma omp simd
          for (int i=pmb->ivs; i<=pmb->ive; i++) {
            var(n, k, j, i) = 0;
            var(n, k, pmb->jps + j, i) = 0;
          }

    }

  } else if (pmb->block_size.nx2 > 1) {
    for (int n=nl_; n<=nu_; ++n) {
      // top and bottom bands
      for (int j=0; j<NGHOST; j++) {
#pragma omp simd
        for (int i=0; i<=pmb->ipe; i++) {
          var(n, 0, j, i) = 0;
          var(n, 0, pmb->jps + j, i) = 0;
        }
      }
      // east and west truncated bands
      for (int j=pmb->jvs; j<=pmb->jve; j++) {
#pragma omp simd
        for (int i=0; i<NGHOST; i++) {
          var(n, 0, j, i) = 0;
          var(n, 0, j, pmb->jps + i) = 0;
        }
      }
    }
  } else {
    for (int n=nl_; n<=nu_; ++n) {
#pragma omp simd
      for (int i=pmb->ims; i<=pmb->ime; ++i) {
        var(n, 0, 0, i) = 0;
      }
#pragma omp simd
      for (int i=pmb->ips; i<=pmb->ipe; ++i) {
        var(n, 0, 0, i) = 0;
      }
    }
  }

  return;
  }


void VertexCenteredBoundaryVariable::_FinalizeVert3(){
  // During refinement, blocks are split and refined blocks are in contact.
  // This function modifies division factors to that of a MeshBlock that isn't
  // refined.

  MeshBlock *pmb = pmy_block_;
  Mesh *pm = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_vc;

  int ng = NGHOST;

  // level of current mesh-block
  int &lev = pmb->loc.level;

  // dimensions
  int bh_nx[3] = {pmb->block_size.nx1 / 2,
                  pmb->block_size.nx2 / 2,
                  pmb->block_size.nx3 / 2};

  int x_cnr[3] = {pmb->ivs + bh_nx[0],
                  pmb->jvs + bh_nx[1],
                  pmb->kvs + bh_nx[2]};


  // int ox1_cnr[2] = {0, pmb->ive};
  // int ox2_cnr[2] = {0, pmb->jve};
  // int ox3_cnr[2] = {0, pmb->kve};

  int fac;

  // store multiplicative factors
  // int nmul[5][5][5];

  // for(int k=0; k<5; k++)
  //   for(int j=0; j<5; j++)
  //     for(int i=0; i<5; i++)
  //       nmul[k][j][i] = 0;

  // iterator over neighbors
  // for (int n=0; n<pbval_->nneighbor; n++) {
  //   NeighborBlock& nb = pbval_->neighbor[n];
  //   NeighborConnect nb_t = nb.ni.type;
  //   if (nb_t == NeighborConnect::face) {
  //     // two of {ox1, ox2, ox3} == 0
  //   } else if (nb_t == NeighborConnect::edge) {
  //     // one of ox1 == 0 | ox2 == 0 | ox3 == 0
  //     // remainder non-zero
  //     // nmul[2 * (nb.ni.ox3 + 1)][2 * (nb.ni.ox2 + 1)][2 * (nb.ni.ox1 + 1)] += 2;
  //   } else if (nb_t == NeighborConnect::corner) {
  //     //nmul[2 * (nb.ni.ox3 + 1)][2 * (nb.ni.ox2 + 1)][2 * (nb.ni.ox1 + 1)] += 4;
  //     // nmul[2 * (nb.ni.ox3 + 1)][2 * (nb.ni.ox2 + 1)][2 * (nb.ni.ox1 + 1)] += 1;
  //     if (nb.snb.level>=lev) {
  //       const int iox3 = (nb.ni.ox3+1)/2;
  //       const int iox2 = (nb.ni.ox2+1)/2;
  //       const int iox1 = (nb.ni.ox1+1)/2;
  //       var(0,
  //           ng * (1 - iox3) + ox3_cnr[iox3],
  //           ng * (1 - iox2) + ox2_cnr[iox2],
  //           ng * (1 - iox1) + ox1_cnr[iox1]) /= 8.;
  //     }
  //   }

  // }

  // printf("nmul\n");
  // for (int k=0; k<5; k++) {
  //   printf("\nk=%d:\n", k); // slab for readability
  //   for (int j=0; j<5; j++) {
  //     for (int i=0; i<5; i++) {
  //       printf("%d, ", nmul[k][j][i]);
  //     }
  //     printf("\n");
  //   }
  // }

  return;
  /*
  //-----ox3:

  // ox3 < 0, west edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][1][0] == lev)    // W lvl
         + 2 * (pmb->pbval->nblevel[0][1][0] > lev) // [for ref.] W lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      var(0, k, x_cnr[1], pmb->ivs) /= fac / 2.;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][1][0]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, x_cnr[1], pmb->ime-ix) /= 2.;
        }

  }

  // ox3 < 0, east edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][1][2] == lev)    // E lvl
         + 2 * (pmb->pbval->nblevel[0][1][2] > lev) // [for ref.] E lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
        var(0, k, x_cnr[1], pmb->ive) /= fac / 2.;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][1][2]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, x_cnr[1], pmb->ips+ix) /= 2.;
        }

  }

  // ox3 < 0, north edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][0][1] == lev)    // N lvl
         + 2 * (pmb->pbval->nblevel[0][0][1] > lev) // [for ref.] N lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      var(0, k, pmb->jvs, x_cnr[0]) /= fac / 2.;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][0][1]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jme-ix, x_cnr[0]) /= 2.;
        }

  }

  // ox3 < 0, south edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][2][1] == lev)    // S lvl
         + 2 * (pmb->pbval->nblevel[0][2][1] > lev) // [for ref.] S lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      var(0, k, pmb->jve, x_cnr[0]) /= fac / 2.;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][2][1]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jps+ix, x_cnr[0]) /= 2.;
        }

  }

  // ox3 < 0, internal cross
  if (pmb->pbval->nblevel[0][1][1] > lev)
    for (int k=pmb->kms; k<pmb->kvs; ++k) {
      for (int ix=pmb->iis; ix<=pmb->iie; ++ix) {
        var(0, k, x_cnr[1], ix) /= 2.;
      }
      for (int ix=pmb->jis; ix<=pmb->jie; ++ix) {
        var(0, k, ix, x_cnr[0]) /= 2.;
      }
    }
  //-------------------------------------------

  return;
  */

  //int fac;


  //-----top-cap ghost-zones:

  // ox3 < 0, west edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][1][0] == lev)    // W lvl
         + 2 * (pmb->pbval->nblevel[0][1][0] > lev) // [for ref.] W lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      var(0, k, x_cnr[1], pmb->ivs) /= fac;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][1][0]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, x_cnr[1], pmb->ime-ix) /= 2.;
        }

    // correct remaining part of edge
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      for (int ix=0; ix<bh_nx[1]-1; ++ix) {
        var(0, k, pmb->jis+ix, pmb->ivs) /= 2.;
        var(0, k, pmb->jie-ix, pmb->ivs) /= 2.;
      }
  }

  // ox3 < 0, east edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][1][2] == lev)    // E lvl
         + 2 * (pmb->pbval->nblevel[0][1][2] > lev) // [for ref.] E lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
        var(0, k, x_cnr[1], pmb->ive) /= fac;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][1][2]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, x_cnr[1], pmb->ips+ix) /= 2.;
        }

    // correct remaining part of edge
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      for (int ix=0; ix<bh_nx[1]-1; ++ix) {
        var(0, k, pmb->jis+ix, pmb->ive) /= 2.;
        var(0, k, pmb->jie-ix, pmb->ive) /= 2.;
      }

  }


  // ox3 < 0, north edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][0][1] == lev)    // N lvl
         + 2 * (pmb->pbval->nblevel[0][0][1] > lev) // [for ref.] N lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      var(0, k, pmb->jvs, x_cnr[0]) /= fac;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][0][1]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jme-ix, x_cnr[0]) /= 2.;
        }

    // correct remaining part of edge
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      for (int ix=0; ix<bh_nx[0]-1; ++ix) {
        var(0, k, pmb->jvs, pmb->iis+ix) /= 2.;
        var(0, k, pmb->jvs, pmb->iie-ix) /= 2.;
      }

  }

  // ox3 < 0, south edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] > lev)     // face lvl
         + (pmb->pbval->nblevel[0][2][1] == lev)    // S lvl
         + 2 * (pmb->pbval->nblevel[0][2][1] > lev) // [for ref.] S lvl
         );

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      var(0, k, pmb->jve, x_cnr[0]) /= fac;

    // correct perp. outward strip
    if (pmb->pbval->nblevel[0][2][1]>lev)
      for (int k=pmb->kms; k<pmb->kvs; ++k)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jps+ix, x_cnr[0]) /= 2.;
        }

    // correct remaining part of edge
    for (int k=pmb->kms; k<pmb->kvs; ++k)
      for (int ix=0; ix<bh_nx[0]-1; ++ix) {
        var(0, k, pmb->jve, pmb->iis+ix) /= 2.;
        var(0, k, pmb->jve, pmb->iie-ix) /= 2.;
      }

  }

  // ox3 < 0, NW corner
  fac = (1
         + (pmb->pbval->nblevel[0][0][0] >= lev)
         + (pmb->pbval->nblevel[0][0][1] >= lev)
         + (pmb->pbval->nblevel[0][1][0] >= lev));

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k) {
      if ((pmb->pbval->nblevel[0][0][0] >= lev))
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jvs, ix) /= 2.;
          var(0, k, ix, pmb->ivs) /= 2.;
        }
      var(0, k, pmb->jvs, pmb->ivs) /= fac;
    }
  }

  // ox3 < 0, NE corner
  fac = (1
         + (pmb->pbval->nblevel[0][0][2] >= lev)
         + (pmb->pbval->nblevel[0][0][1] >= lev)
         + (pmb->pbval->nblevel[0][1][2] >= lev));

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k) {
      if ((pmb->pbval->nblevel[0][0][2] >= lev))
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jvs, pmb->ipe-ix) /= 2.;
          var(0, k, ix, pmb->ive) /= 2.;
        }
      var(0, k, pmb->jvs, pmb->ive) /= fac;
    }
  }


  // ox3 < 0, SW corner
  fac = (1
         + (pmb->pbval->nblevel[0][2][0] >= lev)
         + (pmb->pbval->nblevel[0][1][0] >= lev)
         + (pmb->pbval->nblevel[0][2][1] >= lev));

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k) {
      if (pmb->pbval->nblevel[0][2][0] >= lev)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jve, ix) /= 2.;
          var(0, k, pmb->jpe-ix, pmb->ivs) /= 2.;
        }
      var(0, k, pmb->jve, pmb->ivs) /= fac;
    }
  }

  // ox3 < 0, SE corner
  fac = (1
         + (pmb->pbval->nblevel[0][2][2] >= lev)
         + (pmb->pbval->nblevel[0][2][1] >= lev)
         + (pmb->pbval->nblevel[0][1][2] >= lev));

  if (fac > 1) {
    for (int k=pmb->kms; k<pmb->kvs; ++k) {
      if (pmb->pbval->nblevel[0][2][2] >= lev)
        for (int ix=0; ix<ng; ++ix) {
          var(0, k, pmb->jve, pmb->ipe-ix) /= 2.;
          var(0, k, pmb->jpe-ix, pmb->ive) /= 2.;
        }

      var(0, k, pmb->jve, pmb->ive) /= fac;
    }
  }

  // ox3 < 0, internal cross
  if (pmb->pbval->nblevel[0][1][1] > lev)
    for (int k=pmb->kms; k<pmb->kvs; ++k) {
      for (int ix=pmb->iis; ix<=pmb->iie; ++ix) {
        var(0, k, x_cnr[1], ix) /= 2.;
      }
      for (int ix=pmb->jis; ix<=pmb->jie; ++ix) {
        var(0, k, ix, x_cnr[0]) /= 2.;
      }
    }

  //-------------------------------------------

  //-----top-face:

  // ox3 < 0, west edge split
  fac = (1
         + (pmb->pbval->nblevel[0][1][1] == lev)
         + (pmb->pbval->nblevel[0][1][0] == lev)
         + (pmb->pbval->nblevel[1][1][0] == lev)
         + 2 * (pmb->pbval->nblevel[0][1][1] > lev)
         + 2 * (pmb->pbval->nblevel[0][1][0] > lev)
         + 2 * (pmb->pbval->nblevel[1][1][0] > lev)
         );

  if (fac > 1) {
    var(0, pmb->kvs, x_cnr[1], pmb->ivs) /= fac;
  }
  //-------------------------------------------

}

void VertexCenteredBoundaryVariable::_FinalizeVertexConsistency3(
  int ox1, int ox2, int ox3){

  int ox[3] = {ox1, ox2, ox3};

  MeshBlock *pmb = pmy_block_;
  Mesh *pm = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_vc;

  int ng = NGHOST;
  // int b_nx1 = pmb->block_size.nx1;
  // int b_nx2 = pmb->block_size.nx2;
  // int b_nx3 = pmb->block_size.nx3;

  int &lev = pmb->loc.level;

  // overall approach:
  // initially assume no refinement
  // correct afterwards if false


  // individual logic for faces hard-coded
  //
  // only one oxi is non-zero, thus with (e.g.)
  // ox3 = ox2 = 0, ox1 = 1
  // we can inspect the face-connected neighbor level with:
  // pmb->pbval->nblevel[ox3 + 1][ox2 + 1][ox1 + 1]

  // face - connecting neighbour level
  int nb_lev = pmb->pbval->nblevel[ox[2] + 1][ox[1] + 1][ox[0] + 1];

  printf("(lev, nb_lev; ox1, ox2, ox3) = (%d, %d; %d, %d, %d)\n",
         lev, nb_lev, ox[0], ox[1], ox[2]);

  // don't apply if neighbour is coarser

  // if (nb_lev < lev) {
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();
  // }


  if (nb_lev < lev)
    return;

  // if (nb_lev > lev) {
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();
  // }

  // if ((ox3 == 1)
  //     && (pmb->pbval->nblevel[2][0][1] >= lev)
  //     && (pmb->pbval->nblevel[2][1][0] >= lev)
  //     && (pmb->pbval->nblevel[2][1][2] >= lev)
  //     && (pmb->pbval->nblevel[2][2][0] >= lev)
  //     && (pmb->pbval->nblevel[2][1][1] == lev + 1)) {
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();
  // }

  printf("->cont.\n");

  int s[3] = {pmb->ims, pmb->jms, pmb->kms};
  int e[3] = {pmb->ipe, pmb->jpe, pmb->kpe};

  // build non-ghost ranges
  int svs[3] = {pmb->ivs, pmb->jvs, pmb->kvs};
  int evs[3] = {pmb->ive, pmb->jve, pmb->kve};

  // combine ghost ranging
  int vs[2][3] = {{pmb->ivs, pmb->jvs, pmb->kvs},
                  {pmb->ive, pmb->jve, pmb->kve}};

  // If refinement is enabled then the consistency condition must be further
  // modified. This is due to new interfaces that can appear if contacting a
  // neighbor on a finer level. Moreover, coarser neighbors do not additively
  // unpack and must be ignored.
  //
  // Types of finer connection and treatment:
  //   corner: can be ignored
  //   face: restrict ranges as indicated below

  // [BC / refinement] take care of coarser neighbors not contributing
  if (pmb->pbval->nblevel[1][1][0] < lev) // ox1==-1 face is coarser
    s[0] += ng;
  if (pmb->pbval->nblevel[1][1][2] < lev) // ox1==+1 face is coarser
    e[0] -= ng;

  if (pmb->pbval->nblevel[1][0][1] < lev) // ox2==-1 face is coarser
    s[1] += ng;
  if (pmb->pbval->nblevel[1][2][1] < lev) // ox2==+1 face is coarser
    e[1] -= ng;

  if (pmb->pbval->nblevel[0][1][1] < lev) // ox3==-1 face is coarser
    s[2] += ng;
  if (pmb->pbval->nblevel[2][1][1] < lev) // ox3==+1 face is coarser
    e[2] -= ng;


  // select face by reducing range for salient axis
  int ix_dimf = (ox[2] != 0) ? 2 : (ox[1] != 0) ? 1 : 0;

  // start / end indices confined for face under consideration
  if (ox[ix_dimf] < 0)
    s[ix_dimf] = e[ix_dimf] = svs[ix_dimf];
  else if (ox[ix_dimf] > 0)
    s[ix_dimf] = e[ix_dimf] = evs[ix_dimf];


  // apply averaging condition
  for (int n_=nl_; n_<=nu_; ++n_) {
    for (int k=s[2]; k<=e[2]; ++k)
      for (int j=s[1]; j<=e[1]; ++j)
        for (int i=s[0]; i<=e[0]; ++i) {
          var(n_, k, j, i) /= 2.;
        }
  }


  return;


  // compensate for finer levels touching current face
  if (pm->multilevel) {
    int bh_nx[3] = {pmb->block_size.nx1 / 2,
                    pmb->block_size.nx2 / 2,
                    pmb->block_size.nx3 / 2};
    // internal corner coordinates for face
    int x_cnr[3] = {ng + (ox[0] + 1) * bh_nx[0],
                    ng + (ox[1] + 1) * bh_nx[1],
                    ng + (ox[2] + 1) * bh_nx[2]};

    printf("bh_nx, x_cnr = (%d, %d, %d), (%d, %d, %d)\n",
           bh_nx[0], bh_nx[1], bh_nx[2], x_cnr[0], x_cnr[1], x_cnr[2]);

    // Relative to the current face we additionally need:
    // north: [0][1]
    // east:  [1][2]
    // south: [2][1]
    // west:  [1][0]
    // neighbor levels

    //int nb_N_lev = []

    // Assumption: face is at same level as block

    // n_ iter!

    // treat ox3 axis N/S face edges [off-the-face]
    // int cfix = 2;

    // cfix = 0;
    // if ((nb_lev > lev) && (ox3 == -1)
    //     // && (pmb->pbval->nblevel[cfix][2][1] == lev)
    //     // && (pmb->pbval->nblevel[cfix][2][1] == lev)
    //     // && (pmb->pbval->nblevel[cfix][1][0] == -1)
    //     // && (pmb->pbval->nblevel[cfix][1][2] == -1)
    //     ) {
    //   coutBoldRed("MB::UWIL gid = ");
    //   printf("%d\n", pmb->gid);
    //   Q();
    // }

    // if ((nb_lev > lev) && (ox1 == 1)
    //     // && (pmb->pbval->nblevel[cfix][2][1] == lev)
    //     // && (pmb->pbval->nblevel[cfix][2][1] == lev)
    //     // && (pmb->pbval->nblevel[cfix][1][0] == -1)
    //     // && (pmb->pbval->nblevel[cfix][1][2] == -1)
    //     ) {
    //   coutBoldRed("MB::UWIL gid = ");
    //   printf("%d\n", pmb->gid);
    //   Q();
    // }

    // if ((nb_lev > lev) && (ox2 == -1)
    //     // && (pmb->pbval->nblevel[cfix][2][1] == lev)
    //     // && (pmb->pbval->nblevel[cfix][2][1] == lev)
    //     // && (pmb->pbval->nblevel[cfix][1][0] == -1)
    //     // && (pmb->pbval->nblevel[cfix][1][2] == -1)
    //     ) {
    //   coutBoldRed("MB::UWIL gid = ");
    //   printf("%d\n", pmb->gid);
    //   Q();
    // }


    if (nb_lev > lev) {
      // face is touching partitioned

      if (ox1 != 0) {
        int ix_ox1 = (ox1 + 1) / 2;
        int sgn_ox1 = (ox1 < 1) ? -1 : 1;

        int dI_ns[2] = {
          (pmb->pbval->nblevel[0][1][2*ix_ox1] < nb_lev) ? ng + 1: 0,
          (pmb->pbval->nblevel[2][1][2*ix_ox1] < nb_lev) ? ng + 1: 0
        };

        int dI_we[2] = {
          (pmb->pbval->nblevel[1][0][2*ix_ox1] < nb_lev) ? ng + 1: 0,
          (pmb->pbval->nblevel[1][2][2*ix_ox1] < nb_lev) ? ng + 1: 0
        };

        // loop over ghost layer
        for (int J=0; J<=ng; J++) {
          Real fac = (J==0) ? 1.5 : 2.;
          Real fac_coarser = (J==0) ? 1.25 : 1.5;

          int ix_J = vs[ix_ox1][0] + sgn_ox1 * J;

          // north-south
          for (int I=svs[2] - ng + dI_ns[0]; I<=evs[2] + ng - dI_ns[1]; I++) {
            var(0, I, x_cnr[1], ix_J) /= fac;
          }

          // east-west
          for (int I=svs[1] - ng + dI_we[0]; I<=evs[1] + ng - dI_we[1]; I++) {
            var(0, x_cnr[2], I, ix_J) /= fac;
          }

          for (int I=0; I<2; I++) {
            if (dI_ns[I] > 0) // neighbor not refined [edge-cleanup]
              var(0, vs[I][2], x_cnr[1], ix_J) /= fac_coarser;

            if (dI_we[I] > 0)
              var(0, x_cnr[2], vs[I][1], ix_J) /= fac_coarser;
          }

        }
        // finally correct the centered ghost
        var(0, x_cnr[2], x_cnr[1], vs[ix_ox1][0]) /= 20./18.;

      }

      if (ox2 != 0) {
        int ix_ox2 = (ox2 + 1) / 2;
        int sgn_ox2 = (ox2 < 1) ? -1 : 1;

        int dI_ns[2] = {
          (pmb->pbval->nblevel[0][2*ix_ox2][1] < nb_lev) ? ng + 1: 0,
          (pmb->pbval->nblevel[2][2*ix_ox2][1] < nb_lev) ? ng + 1: 0
        };

        int dI_we[2] = {
          (pmb->pbval->nblevel[1][2*ix_ox2][0] < nb_lev) ? ng + 1: 0,
          (pmb->pbval->nblevel[1][2*ix_ox2][2] < nb_lev) ? ng + 1: 0
        };

        // loop over ghost layer
        for (int J=0; J<=ng; J++) {
          Real fac = (J==0) ? 1.5 : 2.;
          Real fac_coarser = (J==0) ? 1.25 : 1.5;

          int ix_J = vs[ix_ox2][0] + sgn_ox2 * J;

          // north-south
          for (int I=svs[0] - ng + dI_ns[0]; I<=evs[0] + ng - dI_ns[1]; I++) {
            var(0, x_cnr[2], ix_J, I) /= fac;
          }

          // east-west
          for (int I=svs[2] - ng + dI_we[0]; I<=evs[2] + ng - dI_we[1]; I++) {
            var(0, I, ix_J, x_cnr[0]) /= fac;
          }

          for (int I=0; I<2; I++) {
            if (dI_ns[I] > 0) // neighbor not refined [edge-cleanup]
              var(0, x_cnr[2], ix_J, vs[I][0]) /= fac_coarser;

            if (dI_we[I] > 0)
              var(0, vs[I][2], ix_J, x_cnr[0]) /= fac_coarser;
          }

        }
        // finally correct the centered ghost
        var(0, x_cnr[2], vs[ix_ox2][1], x_cnr[0]) /= 20./18.;

      }

      if (ox3 != 0) {
        int ix_ox3 = (ox3 + 1) / 2;
        int sgn_ox3 = (ox3 < 1) ? -1 : 1;

        int dI_ns[2] = {
          (pmb->pbval->nblevel[2*ix_ox3][0][1] < nb_lev) ? ng + 1: 0,
          (pmb->pbval->nblevel[2*ix_ox3][2][1] < nb_lev) ? ng + 1: 0
        };

        int dI_we[2] = {
          (pmb->pbval->nblevel[2*ix_ox3][1][0] < nb_lev) ? ng + 1: 0,
          (pmb->pbval->nblevel[2*ix_ox3][1][2] < nb_lev) ? ng + 1: 0
        };

        // loop over ghost layer
        for (int J=0; J<=ng; J++) {
          Real fac = (J==0) ? 1.5 : 2.;
          Real fac_coarser = (J==0) ? 1.25 : 1.5;

          int ix_J = vs[ix_ox3][0] + sgn_ox3 * J;

          // north-south
          for (int I=svs[1] - ng + dI_ns[0]; I<=evs[1] + ng - dI_ns[1]; I++) {
            var(0, ix_J, I, x_cnr[0]) /= fac;
          }

          // east-west
          for (int I=svs[0] - ng + dI_we[0]; I<=evs[0] + ng - dI_we[1]; I++) {
            var(0, ix_J, x_cnr[1], I) /= fac;
          }

          for (int I=0; I<2; I++) {
            if (dI_ns[I] > 0) // neighbor not refined [edge-cleanup]
              var(0, ix_J, vs[I][1], x_cnr[0]) /= fac_coarser;

            if (dI_we[I] > 0)
              var(0, ix_J, x_cnr[1], vs[I][0]) /= fac_coarser;
          }

        }
        // finally correct the centered ghost
        var(0, vs[ix_ox3][2], x_cnr[1], x_cnr[0]) /= 20./18.;
      }



    } else if (nb_lev == lev) {

      // correct edges with unrefined neighbours
      coutBoldRed("nb_lev == lev\n");

      if (ox1 != 0) {
        int ix_ox1 = (ox1 + 1) / 2;
        int sgn_ox1 = (ox1 < 1) ? -1 : 1;

        // loop over ghost layer
        for (int J=0; J<=ng; J++) {
          Real fac = (J==0) ? 1.5 : 2.;
          Real fac_coarser = (J==0) ? 1.25 : 1.5;

          int ix_J = vs[ix_ox1][0] + sgn_ox1 * J;

          // fac = std::sqrt(fac);
          // fac_coarser = std::sqrt(fac_coarser);

          // edge clean-up required
          for (int dir=0; dir<2; dir++) {
            // north-south
            if (pmb->pbval->nblevel[2*dir][1][2*ix_ox1] > lev) {
              for (int I=1; I<=ng; I++) {
                var(0, vs[dir][2] + (2*dir-1)*I, x_cnr[1], ix_J) /= fac;
              }
              var(0, vs[dir][2], x_cnr[1], ix_J) /= fac_coarser;
            }

            // west-east
            if (pmb->pbval->nblevel[1][2*dir][2*ix_ox1] > lev) {
              for (int I=1; I<=ng; I++) {
                var(0, x_cnr[2], vs[dir][1] + (2*dir-1)*I, ix_J) /= fac;
              }
              var(0, x_cnr[2], vs[dir][1], ix_J) /= fac_coarser;
            }

          }
        }

      }

      if (ox2 != 0) {
        int ix_ox2 = (ox2 + 1) / 2;
        int sgn_ox2 = (ox2 < 1) ? -1 : 1;

        // loop over ghost layer
        for (int J=0; J<=ng; J++) {
          Real fac = (J==0) ? 1.5 : 2.;
          Real fac_coarser = (J==0) ? 1.25 : 1.5;

          int ix_J = vs[ix_ox2][0] + sgn_ox2 * J;

          // fac = std::sqrt(fac);
          // fac_coarser = std::sqrt(fac_coarser);

          // edge clean-up required
          for (int dir=0; dir<2; dir++) {
            // north-south
            if (pmb->pbval->nblevel[2*dir][2*ix_ox2][1] > lev) {
              for (int I=1; I<=ng; I++) {
                var(0, x_cnr[2], ix_J, vs[dir][0] + (2*dir-1)*I) /= fac;
              }
              var(0, x_cnr[2], ix_J, vs[dir][0]) /= fac_coarser;
            }

            // west-east
            if (pmb->pbval->nblevel[1][2*ix_ox2][2*dir] > lev) {
              for (int I=1; I<=ng; I++) {
                var(0, vs[dir][2] + (2*dir-1)*I, ix_J, x_cnr[0]) /= fac;
              }
              var(0, vs[dir][2], ix_J, x_cnr[0]) /= fac_coarser;
            }

          }
        }

      }

      if (ox3 != 0) {
        int ix_ox3 = (ox3 + 1) / 2;
        int sgn_ox3 = (ox3 < 1) ? -1 : 1;

        // loop over ghost layer
        for (int J=0; J<=ng; J++) {
          Real fac = (J==0) ? 1.5 : 2.;
          Real fac_coarser = (J==0) ? 1.25 : 1.5;

          int ix_J = vs[ix_ox3][0] + sgn_ox3 * J;

          // fac = std::sqrt(fac);
          // fac_coarser = std::sqrt(fac_coarser);

          // edge clean-up required
          for (int dir=0; dir<2; dir++) {
            // north-south
            if (pmb->pbval->nblevel[2*ix_ox3][2*dir][1] > lev) {
              for (int I=1; I<=ng; I++) {
                var(0, ix_J, vs[dir][1] + (2*dir-1)*I, x_cnr[0]) /= fac;
              }
              var(0, ix_J, vs[dir][1], x_cnr[0]) /= fac_coarser;
            }

            // west-east
            if (pmb->pbval->nblevel[2*ix_ox3][1][2*dir] > lev) {
              for (int I=1; I<=ng; I++) {
                var(0, ix_J, x_cnr[1], vs[dir][0] + (2*dir-1)*I) /= fac;
              }
              var(0, ix_J, x_cnr[1], vs[dir][0]) /= fac_coarser;
            }
          }

        }
      }

    }


  }

  // int si = pmb->ims, ei = pmb->ipe;
  // int sj = pmb->jms, ej = pmb->jpe;
  // int sk = pmb->kms, ek = pmb->kpe;

  // // for ghost-zones
  // // int gsi = pmb->ims, gei = pmb->ime;
  // // int gsj = pmb->jms, gej = pmb->jme;
  // // int gsk = pmb->kms, gek = pmb->kme;

  // // for treating corner interiors
  // int sI = -1, eI = 1, sJ = -1, eJ = 1, sK = -1, eK = 1;


  // if (ox1 < 0)
  //   si = ei = pmb->ivs, sI = eI = 0;
  // else if (ox1 > 0)
  //   si = ei = pmb->ive, sI = eI = 0;

  // if (ox2 < 0)
  //   sj = ej = pmb->jvs, sJ = eJ = 0;
  // else if (ox2 > 0)
  //   sj = ej = pmb->jve, sJ = eJ = 0;

  // if (ox3 < 0)
  //   sk = ek = pmb->kvs, sK = eK = 0;
  // else if (ox3 > 0)
  //   sk = ek = pmb->kve, sK = eK = 0;

  // for (int n_=nl_; n_<=nu_; ++n_) {
  //   // interior of face
  //   for (int k=sk; k<=ek; ++k)
  //     for (int j=sj; j<=ej; ++j)
  //       for (int i=si; i<=ei; ++i) {
  //         var(n_, k, j, i) /= 2.;
  //       }
  // }



  // // iterate over different corners
  // for (int K=sK; K<=eK; K+=2)
  //   for (int J=sJ; J<=eJ; J+=2)
  //     for (int I=sI; I<=eI; I+=2) {
  //       int cnb_lev = pmb->pbval->nblevel[K + 1][J + 1][I + 1];
  //       printf("cnb_lev = %d\n", cnb_lev);
  //     }

}

void VertexCenteredBoundaryVariable::_FinalizeVert3a(){

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  int &lev = pmb->loc.level;


  const int ng = NGHOST;

  // idx components for addressing block vertices and midpoints
  const int bh_nx[3] = {pmb->block_size.nx3/2,
                        pmb->block_size.nx2/2,
                        pmb->block_size.nx1/2};
  const int k_c[3] = {pmb->kvs, pmb->kvs+bh_nx[2], pmb->kve};
  const int j_c[3] = {pmb->jvs, pmb->jvs+bh_nx[1], pmb->jve};
  const int i_c[3] = {pmb->ivs, pmb->ivs+bh_nx[0], pmb->ive};


  coutBoldRed("pmb->pbval->nblevel\n");
  print_int_arr3(&(pmb->pbval->nblevel[0][0][0]));
  coutBoldRed("lev =");
  printf("%d\n", lev);

  int fac = 1;

  // nb finer or same level
  int nb_rel[3][3][3];
  int nb_ifl[3][3][3];
  //int sl[3][3][3];
  for (int k=0; k<=2; ++k)
    for (int j=0; j<=2; ++j)
      for (int i=0; i<=2; ++i) {
        nb_rel[k][j][i] = pmb->pbval->nblevel[k][j][i] >= lev;
        nb_ifl[k][j][i] = pmb->pbval->nblevel[k][j][i] > lev;
        //sl[k][j][i] = pmb->pbval->nblevel[k][j][i] == lev;
      }

  // extract node multiplicity
  int node_mult[5][5][5] = {{{0}}};  // BD: conv. to short-int

  for (int n=0; n<pbval_->nneighbor; n++) {
    coutMagenta("iterating over neighbors; (n, pbval_->nneighbor)=");
    printf("(%d, %d)\n", n, pbval_->nneighbor);
    NeighborBlock& nb = pbval_->neighbor[n];

    // only refined and equal-level blocks contribute
    if (nb.snb.level >= lev) {
      // retain the type of neighbor connect
      NeighborConnect nb_t = nb.ni.type;

      int ox3 = nb.ni.ox3, ox2 = nb.ni.ox2, ox1 = nb.ni.ox1;

      if (nb_t == NeighborConnect::face) {
        // two of {ox1, ox2, ox3} == 0
        // remainder non-zero
        // printf("f: (ox1, ox2, ox3, l)=(%d, %d, %d, %d)\n",
        //        nb.ni.ox1, nb.ni.ox2, nb.ni.ox3,
        //        pmb->pbval->nblevel[ix_ox3][ix_ox2][ix_ox1]);

      } else if (nb_t == NeighborConnect::edge) {
        // one of ox1 == 0 | ox2 == 0 | ox3 == 0
        // remainder non-zero

        // printf("e: (ox1, ox2, ox3, l)=(%d, %d, %d, %d)\n",
        //        nb.ni.ox1, nb.ni.ox2, nb.ni.ox3,
        //        pmb->pbval->nblevel[ox3+1][iox2][iox1]);


      } else if (nb_t == NeighborConnect::corner) {
        // ox1, ox2, ox3 all non-zero
        int K = 2*(ox3 + 1), J = 2 * (ox2 + 1), I = 2 * (ox1 + 1);
        for (int k=0; k<2; ++k)
          for (int j=0; j<2; ++j)
            for (int i=0; i<2; ++i)
              node_mult[K - k * ox3][J - j * ox2][I - i * ox1]++;

      }

    }

  }


  // coutBoldRed("node_mult\n");
  // print_int_arr3(&node_mult[0][0][0], 5, 5, 5);

  // BD: debug
  // temporarily override multiplicites just using block
  const int tk_c[5] = {pmb->kms, pmb->kvs, pmb->kvs+pmb->block_size.nx3/2,
                       pmb->kve, pmb->kpe};
  const int tj_c[5] = {pmb->jms, pmb->jvs, pmb->jvs+pmb->block_size.nx2/2,
                       pmb->jve, pmb->jpe};
  const int ti_c[5] = {pmb->ims, pmb->ivs, pmb->ivs+pmb->block_size.nx1/2,
                       pmb->ive, pmb->ipe};



  for (int k=0; k<5; ++k)
    for (int j=0; j<5; ++j)
      for (int i=0; i<5; ++i)
        node_mult[k][j][i] = (int) (var(0, tk_c[k], tj_c[j], ti_c[i]) + 0.5);

  //-- print
  for (int k=0; k<5; k++) {
    printf("\nk=%d:\n", k); // slab for readability
    for (int j=0; j<5; j++) {
      for (int i=0; i<5; i++) {
        printf("%d, ", node_mult[k][j][i]);
      }
      printf("\n");
    }
  }


  //----------------------------------------------------------------------------
  // apply node multiplicity factors
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // cube vertices
  //----------------------------------------------------------------------------

  for (int k=0; k<2; ++k)
    for (int j=0; j<2; ++j)
      for (int i=0; i<2; ++i) {
        const int K = 2*(k+1)-1, J = 2*(j+1)-1, I = 2*(i+1)-1;
        var(0, k_c[2*k], j_c[2*j], i_c[2*i]) /= node_mult[K][J][I];
      }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // external ghost corner edges (relative to faces)
  //----------------------------------------------------------------------------

  for (int K=0; K<2; ++K) {

    for (int J=0; J<2; ++J)
      for (int I=0; I<2; ++I) {
        int node_fac = node_mult[2*(K+1)-1][4*J][2*(I+1)-1];
        for (int i=0; i<ng; ++i) {
          const int off = (J-1) * ng + J;
          var(0, k_c[2*K], j_c[2*J]+off+i, i_c[2*I]) /= node_fac;
        }

        node_fac = node_mult[2*(K+1)-1][2*(J+1)-1][4*I];
        for (int i=0; i<ng; ++i) {
          const int off = (I-1) * ng + I;
          var(0, k_c[2*K], j_c[2*J], i_c[2*I]+off+i) /= node_fac;
        }

        node_fac = node_mult[4*K][2*(J+1)-1][2*(I+1)-1];
        for (int i=0; i<ng; ++i) {
          const int off = (K-1) * ng + K;
          var(0, k_c[2*K]+off+i, j_c[2*J], i_c[2*I]) /= node_fac;
        }


      }

  }


  /*
  if (node_mult[1][0][1] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0]-ng+ix, i_c[0]) /= node_mult[1][0][1];

  if (node_mult[1][1][0] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0], i_c[0]-ng+ix) /= node_mult[1][1][0];



  if (node_mult[1][1][4] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0], i_c[2]+1+ix) /= node_mult[1][1][4];

  if (node_mult[1][0][3] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0]-ng+ix, i_c[2]) /= node_mult[1][0][3];



  if (node_mult[1][3][0] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2], i_c[0]-ng+ix) /= node_mult[1][3][0];

  if (node_mult[1][4][1] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2]+1+ix, i_c[0]) /= node_mult[1][4][1];


  if (node_mult[1][3][4] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2], i_c[2]+1+ix) /= node_mult[1][3][4];


  if (node_mult[1][4][3] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2]+1+ix, i_c[2]) /= node_mult[1][4][3];
  */

  //


  /*
  if (node_mult[1][1][0] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0], i_c[0]-ng+ix) /= node_mult[1][1][0];



  if (node_mult[1][1][4] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0], i_c[2]+1+ix) /= node_mult[1][1][4];


  if (node_mult[1][3][0] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2], i_c[0]-ng+ix) /= node_mult[1][3][0];

  if (node_mult[1][3][4] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2], i_c[2]+1+ix) /= node_mult[1][3][4];
  */


  /*
  if (node_mult[1][0][3] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0]-ng+ix, i_c[2]) /= node_mult[1][0][3];

  if (node_mult[1][0][1] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[0]-ng+ix, i_c[0]) /= node_mult[1][0][1];

  if (node_mult[1][4][1] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2]+1+ix, i_c[0]) /= node_mult[1][4][1];

  if (node_mult[1][4][3] > 1)
    for (int ix=0; ix<ng; ++ix)
      var(0, k_c[0], j_c[2]+1+ix, i_c[2]) /= node_mult[1][4][3];
  */

  /*
  for (int K=0; K<2; ++K) {

    for (int J=0; J<2; ++J)
      for (int I=0; I<2; ++I) {
        int node_fac = node_mult[2*(K+1)-1][4*J][2*(I+1)-1];
        for (int i=0; i<ng; ++i) {
          const int off = (J-1) * ng + J;
          var(0, k_c[2*K], j_c[2*J]+off+i, i_c[2*I]) /= node_fac;
        }
      }

    for (int I=0; I<2; ++I)
      for (int J=0; J<2; ++J) {
        int node_fac = node_mult[2*(K+1)-1][2*(J+1)-1][4*I];
        for (int i=0; i<ng; ++i) {
          const int off = (I-1) * ng + I;
          var(0, k_c[2*K], j_c[2*J], i_c[2*I]+off+i) /= node_fac;
        }
      }


  }


  for (int K=0; K<2; ++K) {

    for (int J=0; J<2; ++J)
      for (int I=0; I<2; ++I) {
        int node_fac = node_mult[4*J][2*(K+1)-1][2*(I+1)-1];
        for (int i=0; i<ng; ++i) {
          const int off = (K-1) * ng + K;
          var(0, k_c[2*K]+off+i, j_c[2*J], i_c[2*I]) /= node_fac;
        }
      }

  }
  */



  //----------------------------------------------------------------------------
  // external ghost corners relative to faces
  //----------------------------------------------------------------------------

  for (int K=0; K<2; ++K)
    for (int J=0; J<2; ++J)
      for (int I=0; I<2; ++I) {
        int node_fac = node_mult[2*(K+1)-1][4*J][4*I];
        if (node_fac > 1) {
          for (int j=0; j<ng; ++j)
            for (int i=0; i<ng; ++i) {
              const int jo = (J-1) * ng + J, io = (I-1) * ng + I;
              var(0, k_c[2*K], j_c[2*J]+jo+j, i_c[2*I]+io+i) /= node_fac;
            }
        }
      }

  for (int K=0; K<2; ++K)
    for (int J=0; J<2; ++J)
      for (int I=0; I<2; ++I) {
        int node_fac = node_mult[4*J][2*(K+1)-1][4*I];
        if (node_fac > 1) {
          for (int j=0; j<ng; ++j)
            for (int i=0; i<ng; ++i) {
              const int jo = (J-1) * ng + J, io = (I-1) * ng + I;
              var(0, j_c[2*J]+jo+j, k_c[2*K], i_c[2*I]+io+i) /= node_fac;
            }
        }
      }

  for (int K=0; K<2; ++K)
    for (int J=0; J<2; ++J)
      for (int I=0; I<2; ++I) {
        int node_fac = node_mult[4*J][4*I][2*(K+1)-1];
        if (node_fac > 1) {
          for (int j=0; j<ng; ++j)
            for (int i=0; i<ng; ++i) {
              const int jo = (J-1) * ng + J, io = (I-1) * ng + I;
              var(0, j_c[2*J]+jo+j, i_c[2*I]+io+i, k_c[2*K]) /= node_fac;
            }
        }
      }
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  // external ghost edges relative to faces
  //----------------------------------------------------------------------------

  // north edge

  /*
  fac = node_mult[1][1][2];

  if (fac > 1) {
    int fac_ = (fac % 2 == 0) ? fac : fac - 1;
    for (int ix=0; ix<ng; ++ix) {
      var(0, k_c[0], j_c[0]-ng+ix, i_c[1]) /= (fac - 2);
    }

    for (int ix=0; ix<ng; ++ix) {
      for (int i=1; i<bh_nx[0]; ++i) {
        var(0, k_c[0], j_c[0]-ng+ix, i_c[0]+i) /= (fac_ - 2);
        var(0, k_c[0], j_c[0]-ng+ix, i_c[1]+i) /= (fac_ - 2);
      }
    }

  }

  fac = node_mult[1][3][2];

  if (fac > 1) {
    int fac_ = (fac % 2 == 0) ? fac : fac - 1;
    for (int ix=0; ix<ng; ++ix) {
      var(0, k_c[0], j_c[2]+1+ix, i_c[1]) /= (fac - 2);
    }

    for (int ix=0; ix<ng; ++ix) {
      for (int i=1; i<bh_nx[0]; ++i) {
        var(0, k_c[0], j_c[2]+1+ix, i_c[0]+i) /= (fac_ - 2);
        var(0, k_c[0], j_c[2]+1+ix, i_c[1]+i) /= (fac_ - 2);
      }
    }

  }
  */
  /*

  fac = node_mult[1][2][1];

  if (fac > 1) {
    int fac_ = (fac % 2 == 0) ? fac : fac - 1;
    for (int ix=0; ix<ng; ++ix) {
      var(0, k_c[0], j_c[1], i_c[0]-ng+ix) /= (fac - 2);
    }

    for (int ix=0; ix<ng; ++ix) {
      for (int i=1; i<bh_nx[1]; ++i) {
        var(0, k_c[0], j_c[0]+i, i_c[0]-ng+ix) /= (fac_ - 2);
        var(0, k_c[0], j_c[1]+i, i_c[0]-ng+ix) /= (fac_ - 2);
      }
    }

  }


  fac = node_mult[1][2][3];

  if (fac > 1) {
    int fac_ = (fac % 2 == 0) ? fac : fac - 1;
    for (int ix=0; ix<ng; ++ix) {
      var(0, k_c[0], j_c[1], i_c[2]+1+ix) /= (fac - 2);
    }

    for (int ix=0; ix<ng; ++ix) {
      for (int i=1; i<bh_nx[1]; ++i) {
        var(0, k_c[0], j_c[0]+i, i_c[2]+1+ix) /= (fac_ - 2);
        var(0, k_c[0], j_c[1]+i, i_c[2]+1+ix) /= (fac_ - 2);
      }
    }

  }
  */

  for (int K=0; K<2; ++K)
    for (int I=0; I<2; ++I) {
      // ox3 faces; north-south
      int node_fac = node_mult[2*(K+1)-1][2*(I+1)-1][2];
      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[2*K], j_c[2*I]+off+ix, i_c[1]) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[2]; ++i) {
            var(0, k_c[2*K], j_c[2*I]+off+ix, i_c[0]+i) /= (node_fac_ - 2);
            var(0, k_c[2*K], j_c[2*I]+off+ix, i_c[1]+i) /= (node_fac_ - 2);
          }
        }

      }

      // ox3 faces; east-west
      node_fac = node_mult[2*(K+1)-1][2][2*(I+1)-1];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[2*K], j_c[1], i_c[2*I]+off+ix) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[1]; ++i) {
            var(0, k_c[2*K], j_c[0]+i, i_c[2*I]+off+ix) /= (node_fac_ - 2);
            var(0, k_c[2*K], j_c[1]+i, i_c[2*I]+off+ix) /= (node_fac_ - 2);
          }
        }

      }

      // ox2 faces; north-south
      node_fac = node_mult[2][2*(K+1)-1][2*(I+1)-1];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[1], j_c[2*K], i_c[2*I]+off+ix) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[0]; ++i) {
            var(0, k_c[0]+i, j_c[2*K], i_c[2*I]+off+ix) /= (node_fac_ - 2);
            var(0, k_c[1]+i, j_c[2*K], i_c[2*I]+off+ix) /= (node_fac_ - 2);
          }
        }

      }

      // ox2 faces; east-west
      node_fac = node_mult[2*(I+1)-1][2*(K+1)-1][2];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[2*I]+off+ix, j_c[2*K], i_c[1]) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[2]; ++i) {
            var(0, k_c[2*I]+off+ix, j_c[2*K], i_c[0]+i) /= (node_fac_ - 2);
            var(0, k_c[2*I]+off+ix, j_c[2*K], i_c[1]+i) /= (node_fac_ - 2);
          }
        }

      }

      // ox1 faces; north-south
      node_fac = node_mult[2][2*(I+1)-1][2*(K+1)-1];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[2*I]+off+ix, j_c[1], i_c[2*K]) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[1]; ++i) {
            var(0, k_c[2*I]+off+ix, j_c[0]+i, i_c[2*K]) /= (node_fac_ - 2);
            var(0, k_c[2*I]+off+ix, j_c[1]+i, i_c[2*K]) /= (node_fac_ - 2);
          }
        }

      }

      // ox1 faces; east-west
      node_fac = node_mult[2*(I+1)-1][2][2*(K+1)-1];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[1], j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[0]; ++i) {
            var(0, k_c[0]+i, j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac_ - 2);
            var(0, k_c[1]+i, j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac_ - 2);
          }
        }

      }

      // below are permuted
      // node_fac = node_mult[2][2*(K+1)-1][2*(I+1)-1];

      // if (node_fac > 1) {
      //   int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
      //   const int off = (I-1) * ng + I;

      //   for (int ix=0; ix<ng; ++ix) {
      //     var(0, k_c[1], j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac - 2);
      //   }

      //   for (int ix=0; ix<ng; ++ix) {
      //     for (int i=1; i<bh_nx[1]; ++i) {
      //       var(0, k_c[0]+i, j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac_ - 2);
      //       var(0, k_c[1]+i, j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac_ - 2);
      //     }
      //   }

      // }

      /*
      node_fac = node_mult[2][2*(K+1)-1][2*(I+1)-1];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[1], j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[1]; ++i) {
            var(0, k_c[0]+i, j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac_ - 2);
            var(0, k_c[1]+i, j_c[2*I]+off+ix, i_c[2*K]) /= (node_fac_ - 2);
          }
        }

      }
      */

    }

  if (false)
  for (int K=0; K<2; ++K)
    for (int I=0; I<2; ++I) {
      int node_fac = node_mult[2*(I+1)-1][2*(K+1)-1][2];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[2*I]+off+ix, j_c[2*K], i_c[1]) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[2]; ++i) {
            var(0, k_c[2*I]+off+ix, j_c[2*K], i_c[0]+i) /= (node_fac_ - 2);
            var(0, k_c[2*I]+off+ix, j_c[2*K], i_c[1]+i) /= (node_fac_ - 2);
          }
        }

      }

      node_fac = node_mult[2][2*(K+1)-1][2*(I+1)-1];

      if (node_fac > 1) {
        int node_fac_ = (node_fac % 2 == 0) ? node_fac : node_fac - 1;
        const int off = (I-1) * ng + I;

        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[1], j_c[2*K], i_c[2*I]+off+ix) /= (node_fac - 2);
        }

        for (int ix=0; ix<ng; ++ix) {
          for (int i=1; i<bh_nx[2]; ++i) {
            var(0, k_c[0]+i, j_c[2*K], i_c[2*I]+off+ix) /= (node_fac_ - 2);
            var(0, k_c[1]+i, j_c[2*K], i_c[2*I]+off+ix) /= (node_fac_ - 2);
          }
        }

      }

    }

  //----------------------------------------------------------------------------

  // if (node_mult[1][0][0] > 1) {
  //   for (int jx=0; jx<ng; ++jx)
  //     for (int ix=0; ix<ng; ++ix)
  //       const jo = (J-1) * ng + jx;
  //       var(0, k_c[0], j_c[0]-ng+jx, i_c[0]-ng+ix) /= node_mult[1][0][0];
  // }

  /*
  if (node_mult[1][0][4] > 1) {
    for (int jx=0; jx<ng; ++jx)
      for (int ix=0; ix<ng; ++ix)
        var(0, k_c[0], tj_c[0]+jx, ti_c[3]+1+ix) /= node_mult[1][0][4];
  }


  if (node_mult[1][4][0] > 1) {
    for (int jx=0; jx<ng; ++jx)
      for (int ix=0; ix<ng; ++ix)
        var(0, k_c[0], tj_c[3]+1+jx, ti_c[0]+ix) /= node_mult[1][4][0];
  }

  if (node_mult[1][4][4] > 1) {
    for (int jx=0; jx<ng; ++jx)
      for (int ix=0; ix<ng; ++ix)
        var(0, k_c[0], tj_c[3]+1+jx, ti_c[3]+1+ix) /= node_mult[1][4][4];
  }
  */

  /*
  //----------------------------
  // take care of block vertices

  // k_c[0], j_c[0], i_c[0]

  fac = (nb_rel[0][0][0] +
         nb_rel[0][1][0] +
         nb_rel[0][0][1] +
         nb_rel[0][1][1] +
         nb_rel[1][0][0] +
         nb_rel[1][1][0] +
         nb_rel[1][0][1] +
         nb_rel[1][1][1]);

  if (fac > 1)
    var(0, k_c[0], j_c[0], i_c[0]) /= fac;

  for (int J=0; J<=2; J+=2)
    for (int I=0; I<=2; I+=2) {
      fac = (nb_rel[0][0][0] +
             nb_rel[0][1][0] +
             nb_rel[0][0][1] +
             nb_rel[0][1][1] +
             nb_rel[1][0][0] +
             nb_rel[1][1][0] +
             nb_rel[1][0][1] +
             nb_rel[1][1][1]);

      if (fac > 1)
        var(0, k_c[0], j_c[J], i_c[I]) /= fac;

    }

  // ox3 < 0, NW corner
  fac = (1
         + (pmb->pbval->nblevel[0][0][0] >= lev)
         + (pmb->pbval->nblevel[0][0][1] >= lev)
         + (pmb->pbval->nblevel[0][1][0] >= lev));

  if (fac > 1) {
    if ((pmb->pbval->nblevel[0][0][0] >= lev))
      for (int k=k_c[0]-ng; k<k_c[0]; ++k) {
        for (int ix=0; ix<ng; ++ix) {
          var(0, k_c[0], j_c[0], ix) /= 2.;
          var(0, k_c[0], ix, i_c[0]) /= 2.;
        }
        var(0, k, pmb->jvs, pmb->ivs) /= fac;
    }
  }
  */

  // if (!((std::abs(var(0, k_c[0], j_c[0], i_c[0]) - 1.) < 0.001) &&
  //       (std::abs(var(0, k_c[0], j_c[0], i_c[2]) - 1.) < 0.001) )) {
  //   var.print_all("%1.1f");
  //   coutBoldRed("\n\nerror found:");
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();
  // }

  for (int k=0; k<5; ++k)
    for (int j=0; j<5; ++j)
      for (int i=0; i<5; ++i) {
        if (node_mult[k][j][i] == 3) {

          // var.print_all("%1.1f");

          // [k,j,i]
          // for (int k=pmb->kms; k<=pmb->kpe; ++k) {
          //   printf("\nk=%d:\n", k); // slab for readability
          //   for (int j=pmb->jms; j<=pmb->jpe; ++j) {
          //     for (int i=pmb->ims; i<=pmb->ipe; ++i) {
          //       printf("%1.0f, ", var(0, k, j, i));
          //     }
          //     printf("\n");
          //   }
          // }

          // [j,i,k]
          // for (int j=pmb->jms; j<=pmb->jpe; ++j) {
          //   printf("\nj=%d:\n", j); // slab for readability
          //   for (int i=pmb->ims; i<=pmb->ipe; ++i) {
          //     for (int k=pmb->kms; k<=pmb->kpe; ++k) {
          //       printf("%1.0f, ", var(0, k, j, i));
          //     }
          //     printf("\n");
          //   }
          // }

          // [i,k,j]
          for (int i=pmb->ims; i<=pmb->ipe; ++i) {
            printf("\ni=%d:\n", i); // slab for readability
            for (int k=pmb->kms; k<=pmb->kpe; ++k) {
              for (int j=pmb->jms; j<=pmb->jpe; ++j) {
                printf("%1.0f, ", var(0, k, j, i));
              }
              printf("\n");
            }
          }


          coutBoldRed("\n\nerror found:");
          coutBoldRed("MB::UWIL gid = ");
          printf("%d\n", pmb->gid);
          Q();

        }
      }

  if (false)
  for (int k=0; k<=NGHOST; ++k)
    for (int j=0; j<=pmb->jpe; ++j)
      for (int i=0; i<=pmb->ipe; ++i)
        if (!((std::abs(var(0, k, j, i) - 1.) < 0.001)
              || (std::abs(var(0, k, j, i)) < 0.001))) {
          var.print_all("%1.1f");
          coutBoldRed("\n\nerror found:");
          coutBoldRed("MB::UWIL gid = ");
          printf("%d\n", pmb->gid);
          Q();
        }

  /*
  for (int k=0; k<=NGHOST; ++k)
    for (int j=pmb->jve; j<=pmb->jpe; ++j)
      for (int i=0; i<=pmb->ivs; ++i)
        if (!((std::abs(var(0, k, j, i) - 1.) < 0.001)
              || (std::abs(var(0, k, j, i)) < 0.001))) {
          var.print_all("%1.1f");
          coutBoldRed("\n\nerror found:");
          coutBoldRed("MB::UWIL gid = ");
          printf("%d\n", pmb->gid);
          Q();
        }

  for (int k=0; k<=NGHOST; ++k)
    for (int j=0; j<=pmb->jvs; ++j)
      for (int i=pmb->ive; i<=pmb->ipe; ++i)
        if (!((std::abs(var(0, k, j, i) - 1.) < 0.001)
              || (std::abs(var(0, k, j, i)) < 0.001))) {
          var.print_all("%1.1f");
          coutBoldRed("\n\nerror found:");
          coutBoldRed("MB::UWIL gid = ");
          printf("%d\n", pmb->gid);
          Q();
        }

  for (int k=0; k<=NGHOST; ++k)
    for (int j=pmb->jve; j<=pmb->jpe; ++j)
      for (int i=pmb->ive; i<=pmb->ipe; ++i)
        if (!((std::abs(var(0, k, j, i) - 1.) < 0.001)
              || (std::abs(var(0, k, j, i)) < 0.001))) {
          var.print_all("%1.1f");
          coutBoldRed("\n\nerror found:");
          coutBoldRed("MB::UWIL gid = ");
          printf("%d\n", pmb->gid);
          Q();
        }
  */

}

void VertexCenteredBoundaryVariable::_FinalizeVert3noref(){
    _FinalizeVertexConsistency3(0, 0, -1);
    _FinalizeVertexConsistency3(0, 0, 1);

    _FinalizeVertexConsistency3(0, -1, 0);
    _FinalizeVertexConsistency3(0, 1, 0);

    _FinalizeVertexConsistency3(-1, 0, 0);
    _FinalizeVertexConsistency3(1, 0, 0);
}

void VertexCenteredBoundaryVariable::_FinalizeVert2(){

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  int &mylevel = pmb->loc.level;

  int ng = NGHOST;
  int bh_nx1 = pmb->block_size.nx1 / 2;
  int bh_nx2 = pmb->block_size.nx2 / 2;

  int node_mult[3][3][3];

  for(int i=0; i<=2; i++)
    for(int j=0; j<=2; j++)
      for(int k=0; k<=2; k++)
        node_mult[k][j][i] = 1;

  NeighborConnect nb_type[3][3][3];


  for (int n=0; n<pbval_->nneighbor; n++) {
    // BD: debug
    coutMagenta("iterating over neighbors; (n, pbval_->nneighbor)=");
    printf("(%d, %d)\n", n, pbval_->nneighbor);
    //--

    NeighborBlock& nb = pbval_->neighbor[n];
    int ox1 = nb.ni.ox1, ox2 = nb.ni.ox2;

    // infer information about shared vertex multiplicity

    // common to same and finer
    if (mylevel <= nb.snb.level) {
      // neighbor at the same or finer level

      // retain the type of neighbor connect
      nb_type[0][ox2 + 1][ox1 + 1] = nb.ni.type;

      // retain refinement information
      // nb_fi1[0][ox2 + 1][ox1 + 1] = nb.ni.fi1;

      // central, shared node of neighbor [edge or face]
      node_mult[0][ox2 + 1][ox1 + 1] += 1;

      // north-south contributions to 2d corner
      if (ox1 == 0) {
        node_mult[0][ox2 + 1][0] += 1;
        node_mult[0][ox2 + 1][2] += 1;

        // finer level -> only one "sub-block"
        if (mylevel < nb.snb.level) {
          if (nb.ni.fi1 == 1) {
            node_mult[0][ox2 + 1][0] -= 1;
          } else {
            node_mult[0][ox2 + 1][2] -= 1;
          }
        }
      }

      // east-west
      if (ox2 == 0) {
        node_mult[0][0][ox1 + 1] += 1;
        node_mult[0][2][ox1 + 1] += 1;

        // finer level -> only one "sub-block"
        if (mylevel < nb.snb.level) {
          if (nb.ni.fi1 == 1) {
            node_mult[0][0][ox1 + 1] -= 1;
          } else {
            node_mult[0][2][ox1 + 1] -= 1;
          }
        }
      }

    }
  }

  // apply conditions

  for (int j=0; j<=2; j++) {
    for (int i=0; i<=2; i++) {
      int nblev = pmb->pbval->nblevel[1][j][i];

      int ix_cnr = ng + i * bh_nx1;
      int jx_cnr = ng + j * bh_nx2;

      // apply conditions on internal shared centers
      if (node_mult[0][j][i] > 1) {
        for (int n_=nl_; n_<=nu_; ++n_) {
          var(n_, 0, jx_cnr, ix_cnr) /= node_mult[0][j][i];
        }
      }

      if (mylevel <= nblev) {
        int ox1 = i-1, ox2 = j-1;
        int sgn_ox1 = (ox1 < 0) ? -1 : 1;
        int sgn_ox2 = (ox2 < 0) ? -1 : 1;

        if (nb_type[0][j][i] == NeighborConnect::edge) {
          for (int n_=nl_; n_<=nu_; ++n_) {
            for (int ix=1; ix<=ng; ix++) {
              var(n_, 0, jx_cnr, ix_cnr + sgn_ox1 * ix) /= 2.;
              var(n_, 0, jx_cnr + sgn_ox2 * ix, ix_cnr) /= 2.;
            }
          }
        } else if (nb_type[0][j][i] == NeighborConnect::face) {

          if (ox1 == 0) {
            for (int n_=nl_; n_<=nu_; ++n_) {
              for (int ix=1; ix<bh_nx1; ix++) {
                var(n_, 0, jx_cnr, ix_cnr - ix) /= 2.;
                var(n_, 0, jx_cnr, ix_cnr + ix) /= 2.;
              }
            }

          } else { // ox1 != 0 and ox2 = 0
            for (int n_=nl_; n_<=nu_; ++n_) {
              for (int ix=1; ix<bh_nx2; ix++) {
                var(n_, 0, jx_cnr - ix, ix_cnr) /= 2.;
                var(n_, 0, jx_cnr + ix, ix_cnr) /= 2.;
              }
            }

          }

          if (mylevel < nblev) {
            if (ox1 == 0) {
              for (int n_=nl_; n_<=nu_; ++n_) {
                for (int ix=1; ix<=ng; ix++) {
                  var(n_, 0, jx_cnr + sgn_ox2 * ix, ix_cnr) /= 2.;
                }
              }

            } else { // ox1 != 0 and ox2 = 0
              for (int n_=nl_; n_<=nu_; ++n_) {
                for (int ix=1; ix<=ng; ix++) {
                  var(n_, 0, jx_cnr, ix_cnr + sgn_ox1 * ix) /= 2.;
                }
              }

            }


          }

        }

      }

    }
  }



}

void VertexCenteredBoundaryVariable::_FinalizeVert1(){

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  int &lev = pmb->loc.level;

  // idx components for addressing block vertices and midpoints
  const int i_c[3] = {pmb->ivs, pmb->ivs+pmb->block_size.nx1/2, pmb->ive};

  // 1d
  for (int n=0; n<pbval_->nneighbor; n++) {
    // Same and finer level neighbours additively unpack;
    // this leads to nodes with multiplicity

    NeighborBlock& nb = pbval_->neighbor[n];
    if (lev <= nb.snb.level) {
      // nb.ni.ox1 < 0 => pmb->ivs; nb.ni.ox1 > 0 => pmb->ive
      for (int n_=nl_; n_<=nu_; ++n_)
        var(n_, 0, 0, i_c[nb.ni.ox1+1]) /= 2.;
    }

  }

}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::FinalizeVertexConsistency()
//  \brief apply division factors to shared vertices

void VertexCenteredBoundaryVariable::FinalizeVertexConsistency() {
  coutYellow("VertexCenteredBoundaryVariable::FinalizeVertexConsistency\n");

  // boundary set; impose (dimensionally dependent) consistency condition
  // NOTE: these only make sense for periodic BC
  // see bvals_fc for approach to counting connectivity
  MeshBlock *pmb = pmy_block_;
  Mesh *pm = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_vc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  int &mylevel = pmb->loc.level;

  // DEBUG
  coutMagenta("mylevel: ");
  printf("%d\n", mylevel);

  printf("pre-application:\n");
  var.print_all("%1.0f");

  if(pmb->gid == 3)
    Q();

  // different conditions in different dimensions
  if (pmb->block_size.nx3 > 1)
    _FinalizeVert3a();
    // _FinalizeVert3noref();
  else if (pmb->block_size.nx2 > 1)
    _FinalizeVert2();
  else
    _FinalizeVert1();

  return;



  //
  printf("{{{\n");
  //var.print_all("%1.1f");

  // [k,j,i]
  // for (int k=pmb->kms; k<=pmb->kpe; ++k) {
  //   printf("\nk=%d:\n", k); // slab for readability
  //   for (int j=pmb->jms; j<=pmb->jpe; ++j) {
  //     for (int i=pmb->ims; i<=pmb->ipe; ++i) {
  //       printf("%1.1f, ", var(0, k, j, i));
  //     }
  //     printf("\n");
  //   }
  // }

  // [j,i,k]
  // for (int j=pmb->jms; j<=pmb->jpe; ++j) {
  //   printf("\nj=%d:\n", j); // slab for readability
  //   for (int i=pmb->ims; i<=pmb->ipe; ++i) {
  //     for (int k=pmb->kms; k<=pmb->kpe; ++k) {
  //       printf("%1.1f, ", var(0, k, j, i));
  //     }
  //     printf("\n");
  //   }
  // }

  // [i,k,j]
  // for (int i=pmb->ims; i<=pmb->ipe; ++i) {
  //   printf("\ni=%d:\n", i); // slab for readability
  //   for (int k=pmb->kms; k<=pmb->kpe; ++k) {
  //     for (int j=pmb->jms; j<=pmb->jpe; ++j) {
  //       printf("%1.1f, ", var(0, k, j, i));
  //     }
  //     printf("\n");
  //   }
  // }

  // flip idx for readability
  // var.print_all("%1.1f", false, true);
  // coarse_var.print_all();
  // if (mylevel == 2)
  // Q();

  coutBoldRed("MB::UWIL gid = ");
  printf("%d\n", pmb->gid);

  // printf("x1f: ");
  // pmb->pcoord->x1f.print_data("%1.2f");
  // printf("x2f: ");
  // pmb->pcoord->x2f.print_data("%1.2f");

  /////////////////////////////
  // printf("node_mult: \n");

  int kmax = (pmb->block_size.nx3 > 1) ? 2 : 0;
  // for (int k=0; k<=kmax; k++) {
  //   int kk = 0;
  //   if (pmb->block_size.nx3 > 1) {
  //     printf("\nk=%d:\n", k); // slab for readability
  //     kk = k;
  //   }

  //   for (int j=0; j<=2; j++) {
  //     for (int i=0; i<=2; i++) {
  //       printf("%d, ", node_mult[kk][j][i]);
  //     }
  //     printf("\n");
  //   }
  // }
  /////////////////////////////


  /////////////////////////////
  printf("pbvals->nblevel\n");

  for (int k=0; k<=kmax; k++) {
    int kk = 1;
    if (pmb->block_size.nx3 > 1) {
      printf("\nk=%d:\n", k); // slab for readability
      kk = k;
    }

    for (int j=0; j<=2; j++) {
      for (int i=0; i<=2; i++) {
        printf("%d, ", pmb->pbval->nblevel[kk][j][i]);
        // if (pmb->pbval->nblevel[1][j][i] == -1) {
        //   Q();
        // }
      }
      printf("\n");
    }
  }

  // reorder
  // for (int j=0; j<=2; j++) {
  //   printf("\nj=%d:\n", j); // slab for readability
  //   for (int k=0; k<=2; k++) {
  //     for (int i=0; i<=2; i++) {
  //       printf("%d, ", pmb->pbval->nblevel[k][j][i]);
  //     }
  //     printf("\n");
  //   }
  // }

  // // reorder
  // for (int i=0; i<=2; i++) {
  //   printf("\ni=%d:\n", i); // slab for readability
  //   for (int k=0; k<=2; k++) {
  //     for (int j=0; j<=2; j++) {
  //       printf("%d, ", pmb->pbval->nblevel[k][j][i]);
  //     }
  //     printf("\n");
  //   }
  // }


  /////////////////////////////
  coutMagenta("mylevel: ");
  printf("%d\n", mylevel);

  // filter error on subsets of meshblock

  // for (int k=0; k<pmb->kvs; k++) {
  //   printf("\nk=%d:\n", k); // slab for readability
  //   for (int j=pmb->jms; j<=pmb->jpe; j++) {
  //     for (int i=pmb->ims; i<=pmb->ipe; i++) {
  //       printf("%1.1f, ", var(0, k, j, i));
  //       if (std::abs(var(0, k, j, i) - 1.) > 0.001) {
  //         coutBoldRed("\n\nerror found:");
  //         printf(" @ (%d,%d,%d)\n\n",
  //                k, j, i);
  //         coutBoldRed("MB::UWIL gid = ");
  //         printf("%d\n", pmb->gid);
  //         Q();

  //       }
  //     }
  //     printf("\n");
  //   }
  // }

  // if (std::abs(var(0, 0, 4, 6) - 1.) > 0.001) {
  //   coutBoldRed("\n\nerror found:");
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();

  // }

  // for (int k=0; k<NGHOST; ++k)
  //   for (int ix=0; ix<=NGHOST; ++ix)
  //     if (std::abs(var(0, k, pmb->jvs + pmb->block_size.nx2/2, pmb->ivs-ix) - 1.) > 0.001) {
  //       coutBoldRed("\n\nerror found:");
  //       coutBoldRed("MB::UWIL gid = ");
  //       printf("%d\n", pmb->gid);
  //       Q();
  //     }


  // for (int k=0; k<=NGHOST; ++k)
  //   for (int j=pmb->jms; j<=pmb->jpe; ++j)
  //     for (int i=pmb->ims; i<=pmb->ipe; ++i)
  //       if ((std::abs(var(0, k, j, i) - 1.) > 0.001)
  //           && (std::abs(var(0, k, j, i)) > 0.001)) {
  //         coutBoldRed("\n\nerror found:");
  //         coutBoldRed("MB::UWIL gid = ");
  //         printf("%d\n", pmb->gid);
  //         Q();
  //       }

  // if ((std::abs(var(0, pmb->kvs,
  //                   pmb->jvs+pmb->block_size.nx2 / 2,
  //                   pmb->ivs) - 1.) > 0.001)) {
  //   coutBoldRed("\n\nerror found:");
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();
  // }

  // if ((std::abs(var(0, pmb->kvs,
  //                   pmb->jvs+pmb->block_size.nx2 / 2,
  //                   pmb->ivs) - 1.) > 0.001)) {
  //   coutBoldRed("\n\nerror found:");
  //   coutBoldRed("MB::UWIL gid = ");
  //   printf("%d\n", pmb->gid);
  //   Q();
  // }

  // if(pmb->gid == 5)
  //   Q();
  // if(pmb->gid == 515)
  //   Q();

  // if(pmb->gid == 96)
  //   Q();

  // if(pmb->gid == 351)
  //   Q();
  // if(pmb->gid == 3)
  //   Q();
  // if(pmb->gid == 349)
  //   Q();

  // analyse top-cap
  // if(pmb->gid == 494)
  //   Q();

  // ox3 == -1
  // if(pmb->gid == 494)
  //   Q();

  // ox3 == +1
  // if(pmb->gid == 3)
  //   Q();

  // ox2 == +1
  // if(pmb->gid == 5)
  //   Q();

  // ox2 == -1
  // if(pmb->gid == 362)
  //   Q();

  // ox1 == -1
  // if(pmb->gid == 267)
  //   Q();

  // ox1 == +1
  // if(pmb->gid == 6)
  //   Q();

  // if(pmb->gid == 4)
  //   Q();
  printf("}}}\n");
  return;
}
