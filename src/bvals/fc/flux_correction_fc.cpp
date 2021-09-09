//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_fc.cpp
//! \brief functions that perform flux correction for face-centered variables

// C headers

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../globals.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/buffer_utils.hpp"
#include "bvals_fc.hpp"

// this is not added in flux_correction_cc.cpp:
// #include "../bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//! \todo (felker):
//! - break-up the long functions in this file

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadFluxBoundaryBufferSameLevel(Real *buf,
//!                                                   const NeighborBlock& nb)
//! \brief Set EMF correction buffers for sending to a block on the same level

int FaceCenteredBoundaryVariable::LoadFluxBoundaryBufferSameLevel(
    Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  OrbitalAdvection *porb = pmb->porb;

  // KGF: shearing box:
  AthenaArray<Real> &bx1 = pmb->pfield->b.x1f;
  Real qomL = pbval_->qomL_;

  int p = 0;
  if (nb.ni.type == NeighborConnect::face) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        // pack e2
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++] = e2(k,j,i);
        }
        // pack e3
        // KGF: shearing box
        // shift azmuthal velocity if shearing boundary blocks
        if (nb.shear && pmb->pmy_mesh->sts_loc == TaskType::main_int) {
          if(porb->orbital_advection_defined) {
            for (int k=pmb->ks; k<=pmb->ke; k++) {
              for (int j=pmb->js; j<=pmb->je+1; j++)
                buf[p++] = e3(k,j,i);
            }
          } else {
            if (nb.fid == BoundaryFace::inner_x1) {
              for (int k=pmb->ks; k<=pmb->ke; k++) {
                for (int j=pmb->js; j<=pmb->je+1; j++)
                  buf[p++] = e3(k,j,i) - 0.5*qomL*(bx1(k,j,i) + bx1(k,j-1,i));
              }
            } else if (nb.fid == BoundaryFace::outer_x1) {
              for (int k=pmb->ks; k<=pmb->ke; k++) {
                for (int j=pmb->js; j<=pmb->je+1; j++)
                  buf[p++] = e3(k,j,i) + 0.5*qomL*(bx1(k,j,i) + bx1(k,j-1,i));
              }
            }
          }
        } else {
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              buf[p++] = e3(k,j,i);
          }
        } // KGF: shearing box
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
          int j;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        // pack e1
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            buf[p++] = e1(k,j,i);
        }
        // pack e3
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++] = e3(k,j,i);
        }
        // x3 direction
      } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
        int k;
        if (nb.fid == BoundaryFace::inner_x3) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        // pack e1
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            buf[p++] = e1(k,j,i);
        }
        // pack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            buf[p++] = e2(k,j,i);
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k = pmb->ks;
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        // pack e2
        if (pbval_->shearing_box == 2 && nb.shear
            && pmb->pmy_mesh->sts_loc == TaskType::main_int) {
          if (nb.fid == BoundaryFace::inner_x1) {
            for (int j=pmb->js; j<=pmb->je; j++)
              buf[p++] = e2(k,j,i) + qomL*bx1(k,j,i);
          } else if (nb.fid == BoundaryFace::outer_x1) {
            for (int j=pmb->js; j<=pmb->je; j++)
              buf[p++] = e2(k,j,i) - qomL*bx1(k,j,i);
          }
        } else {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++] = e2(k,j,i);
        }
        // pack e3
        // KGF: shearing box
        // shift azmuthal velocity if shearing boundary blocks
        if (pbval_->shearing_box == 1 && nb.shear
            && pmb->pmy_mesh->sts_loc == TaskType::main_int) {
          if(porb->orbital_advection_defined) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              buf[p++] = e3(k,j,i);
          } else {
            if (nb.fid == BoundaryFace::inner_x1) {
              for (int j=pmb->js; j<=pmb->je+1; j++)
                buf[p++] = e3(k,j,i) - 0.5*qomL*(bx1(k,j,i) + bx1(k,j-1,i));
            } else if (nb.fid == BoundaryFace::outer_x1) {
              for (int j=pmb->js; j<=pmb->je+1; j++)
                buf[p++] = e3(k,j,i) + 0.5*qomL*(bx1(k,j,i) + bx1(k,j-1,i));
            }
          }
        } else {
          for (int j=pmb->js; j<=pmb->je+1; j++)
            buf[p++] = e3(k,j,i);
        } // KGF: shearing box
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        // pack e1
        for (int i=pmb->is; i<=pmb->ie; i++)
          buf[p++] = e1(k,j,i);
        // pack e3
        for (int i=pmb->is; i<=pmb->ie+1; i++)
          buf[p++] = e3(k,j,i);
      }
    } else { // 1D
      int i, j = pmb->js, k = pmb->ks;
      if (nb.fid == BoundaryFace::inner_x1) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }

      // pack e2
      if (pbval_->shearing_box == 2 && nb.shear
          && pmb->pmy_mesh->sts_loc == TaskType::main_int) {
        if (nb.fid == BoundaryFace::inner_x1) {
          buf[p++] = e2(k,j,i) + qomL*bx1(k,j,i);
        } else if (nb.fid == BoundaryFace::outer_x1) {
          buf[p++] = e2(k,j,i) - qomL*bx1(k,j,i);
        }
      } else {
        buf[p++] = e2(k,j,i);
      }

      // pack e3
      buf[p++] = e3(k,j,i);
    }
  } else if (nb.ni.type == NeighborConnect::edge) {
    // x1x2 edge (both 2D and 3D)
    if (nb.eid >= 0 && nb.eid < 4) {
      int i, j;
      if ((nb.eid & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((nb.eid & 2) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      // pack e3
      // KGF: shearing box
      // shift azmuthal velocity if shearing boundary blocks
      if (pbval_->shearing_box == 1
          && nb.shear && nb.ni.ox1 != 0
          && pmb->pmy_mesh->sts_loc == TaskType::main_int) {
        if(porb->orbital_advection_defined) {
          for (int k=pmb->ks; k<=pmb->ke; k++)
            buf[p++] = e3(k,j,i);
        } else {
          if (nb.ni.ox1 == -1) {
            for (int k=pmb->ks; k<=pmb->ke; k++)
              buf[p++] = e3(k,j,i)  - 0.5*qomL*(bx1(k,j,i) + bx1(k,j-1,i));
          } else if (nb.ni.ox1 == 1) {
            for (int k=pmb->ks; k<=pmb->ke; k++)
              buf[p++] = e3(k,j,i)  + 0.5*qomL*(bx1(k,j,i) + bx1(k,j-1,i));
          }
        }
      } else {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          buf[p++] = e3(k,j,i);
      }         // KGF: shearing box
      // x1x3 edge
    } else if (nb.eid >= 4 && nb.eid < 8) {
      int i, k;
      if ((nb.eid & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((nb.eid & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      // pack e2
      if (pbval_->shearing_box == 2
          && nb.shear && nb.ni.ox1 != 0
          && pmb->pmy_mesh->sts_loc == TaskType::main_int) {
        if (nb.ni.ox1 == -1) {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++] = e2(k,j,i) + qomL*bx1(k,j,i);
        } else if (nb.ni.ox1 == 1) {
          for (int j=pmb->js; j<=pmb->je; j++)
            buf[p++] = e2(k,j,i) - qomL*bx1(k,j,i);
        }
      } else {
        for (int j=pmb->js; j<=pmb->je; j++)
          buf[p++] = e2(k,j,i);
      }
      // x2x3 edge
    } else if (nb.eid >= 8 && nb.eid < 12) {
      int j, k;
      if ((nb.eid & 1) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      if ((nb.eid & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      // pack e1
      for (int i=pmb->is; i<=pmb->ie; i++)
        buf[p++] = e1(k,j,i);
    }
  }
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadFluxBoundaryBufferToCoarser(Real *buf,
//!                                                        const NeighborBlock& nb)
//! \brief Set EMF correction buffers for sending to a block on the coarser level

int FaceCenteredBoundaryVariable::LoadFluxBoundaryBufferToCoarser(
    Real *buf, const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  // use the surface area aray as the edge length array
  AthenaArray<Real> &le1 = pbval_->sarea_[0];
  AthenaArray<Real> &le2 = pbval_->sarea_[1];
  int p = 0;
  if (nb.ni.type == NeighborConnect::face) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        // restrict and pack e2
        for (int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          for (int j=pmb->js; j<=pmb->je; j+=2) {
            Real el1 = pco->GetEdge2Length(k,j,i);
            Real el2 = pco->GetEdge2Length(k,j+1,i);
            buf[p++] = (e2(k,j,i)*el1 + e2(k,j+1,i)*el2)/(el1 + el2);
          }
        }
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          for (int j=pmb->js; j<=pmb->je+1; j+=2) {
            bool pole = pco->IsPole(j);
            Real el1, el2;
            if (!pole) {
              el1 = pco->GetEdge3Length(k,j,i);
              el2 = pco->GetEdge3Length(k+1,j,i);
            } else {
              el1 = pco->dx3f(k);
              el2 = pco->dx3f(k+1);
            }
            buf[p++] = (e3(k,j,i)*el1 + e3(k+1,j,i)*el2)/(el1 + el2);
          }
        }
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        bool pole = pco->IsPole(j);
        // restrict and pack e1
        for (int k=pmb->ks; k<=pmb->ke+1; k+=2) {
          if (!pole || !GENERAL_RELATIVITY) {
            pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          } else {
            for (int i = pmb->is; i <= pmb->ie+1; i+=2) {
              le1(i) = pco->dx1f(i);
              le1(i+1) = pco->dx1f(i+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++] = (e1(k,j,i)*le1(i) + e1(k,j,i+1)*le1(i+1))/(le1(i) + le1(i+1));
        }
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          if (!pole) {
            pco->Edge3Length(k,   j, pmb->is, pmb->ie+1, le1);
            pco->Edge3Length(k+1, j, pmb->is, pmb->ie+1, le2);
          } else {
            for (int i = pmb->is; i <= pmb->ie+1; i+=2) {
              le1(i) = pco->dx3f(k);
              le2(i) = pco->dx3f(k+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++] = (e3(k,j,i)*le1(i) + e3(k+1,j,i)*le2(i))/(le1(i) + le2(i));
        }
        // x3 direction
      } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
        int k;
        if (nb.fid == BoundaryFace::inner_x3) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        // restrict and pack e1
        for (int j=pmb->js; j<=pmb->je+1; j+=2) {
          bool pole = pco->IsPole(j);
          if (!pole || !GENERAL_RELATIVITY) {
            pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
          } else {
            for (int i = pmb->is; i <= pmb->ie; i+=2) {
              le1(i) = pco->dx1f(i);
              le1(i+1) = pco->dx1f(i+1);
            }
          }
          for (int i=pmb->is; i<=pmb->ie; i+=2)
            buf[p++] = (e1(k,j,i)*le1(i) + e1(k,j,i+1)*le1(i+1))/(le1(i) + le1(i+1));
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          pco->Edge2Length(k,   j, pmb->is, pmb->ie+1, le1);
          pco->Edge2Length(k, j+1, pmb->is, pmb->ie+1, le2);
          for (int i=pmb->is; i<=pmb->ie+1; i+=2)
            buf[p++] = (e2(k,j,i)*le1(i) + e2(k,j+1,i)*le2(i))/(le1(i) + le2(i));
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k = pmb->ks;
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1 = pco->GetEdge2Length(k,j,i);
          Real el2 = pco->GetEdge2Length(k,j+1,i);
          buf[p++] = (e2(k,j,i)*el1 + e2(k,j+1,i)*el2)/(el1 + el2);
        }
        // pack e3
        for (int j=pmb->js; j<=pmb->je+1; j+=2)
          buf[p++] = e3(k,j,i);
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        bool pole = pco->IsPole(j);
        // restrict and pack e1
        if (!pole || !GENERAL_RELATIVITY) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        } else {
          for (int i = pmb->is; i <= pmb->ie; i+=2) {
            le1(i) = pco->dx1f(i);
            le1(i+1) = pco->dx1f(i+1);
          }
        }
        for (int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++] = (e1(k,j,i)*le1(i) + e1(k,j,i+1)*le1(i+1))/(le1(i) + le1(i+1));
        // pack e3
        for (int i=pmb->is; i<=pmb->ie+1; i+=2)
          buf[p++] = e3(k,j,i);
      }
    } else { // 1D
      int i, j=pmb->js, k=pmb->ks;
      if (nb.fid == BoundaryFace::inner_x1) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      // pack e2 and e3
      buf[p++] = e2(k,j,i);
      buf[p++] = e3(k,j,i);
    }
  } else if (nb.ni.type == NeighborConnect::edge) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if (nb.eid >= 0 && nb.eid < 4) {
        int i, j;
        if ((nb.eid & 1) == 0) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        if ((nb.eid & 2) == 0) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        bool pole = pco->IsPole(j);
        // restrict and pack e3
        for (int k=pmb->ks; k<=pmb->ke; k+=2) {
          Real el1, el2;
          if (!pole) {
            el1 = pco->GetEdge3Length(k,j,i);
            el2 = pco->GetEdge3Length(k+1,j,i);
          } else {
            el1 = pco->dx3f(k);
            el2 = pco->dx3f(k+1);
          }
          buf[p++] = (e3(k,j,i)*el1 + e3(k+1,j,i)*el2)/(el1 + el2);
        }
        // x1x3 edge
      } else if (nb.eid >= 4 && nb.eid < 8) {
        int i, k;
        if ((nb.eid & 1) == 0) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        if ((nb.eid & 2) == 0) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        // restrict and pack e2
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          Real el1 = pco->GetEdge2Length(k,j,i);
          Real el2 = pco->GetEdge2Length(k,j+1,i);
          buf[p++] = (e2(k,j,i)*el1 + e2(k,j+1,i)*el2)/(el1 + el2);
        }
        // x2x3 edge
      } else if (nb.eid >= 8 && nb.eid < 12) {
        int j, k;
        if ((nb.eid & 1) == 0) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        bool pole = pco->IsPole(j);
        if ((nb.eid & 2) == 0) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        // restrict and pack e1
        if (!pole || !GENERAL_RELATIVITY) {
          pco->Edge1Length(k, j, pmb->is, pmb->ie, le1);
        } else {
          for (int i = pmb->is; i <= pmb->ie; i+=2) {
            le1(i) = pco->dx1f(i);
            le1(i+1) = pco->dx1f(i+1);
          }
        }
        for (int i=pmb->is; i<=pmb->ie; i+=2)
          buf[p++] = (e1(k,j,i)*le1(i) + e1(k,j,i+1)*le1(i+1))/(le1(i) + le1(i+1));
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      // x1x2 edge
      int i, j;
      if ((nb.eid & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((nb.eid & 2) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      // pack e3
      buf[p++] = e3(pmb->ks,j,i);
    }
  }
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int FaceCenteredBoundaryVariable::LoadFluxBoundaryBufferToPolar(Real *buf,
//!                                                        const SimpleNeighborBlock &nb,
//!                                                          bool is_north)
//! \brief Load EMF values along polar axis into send buffers

int FaceCenteredBoundaryVariable::LoadFluxBoundaryBufferToPolar(
    Real *buf, const SimpleNeighborBlock &nb, bool is_north) {
  MeshBlock *pmb = pmy_block_;
  int count = 0;
  int j = is_north ? pmb->js : pmb->je + 1;
  for (int i = pmb->is; i <= pmb->ie; ++i) {
    Real val = 0.0;
    for (int k = pmb->ks; k <= pmb->ke; ++k) {  // avoid double counting right ends
      val += pmb->pfield->e.x1e(k, j, i);
    }
    buf[count++] = val / (pmb->ke - pmb->ks + 1);
  }
  return count;
}

//----------------------------------------------------------------------------------------
//! \brief helper function for below SendFluxCorrection()

void FaceCenteredBoundaryVariable::CopyPolarBufferSameProcess(
    const SimpleNeighborBlock& snb, int ssize, int polar_block_index, bool is_north) {
  // Locate target buffer
  // 1) which MeshBlock?
  MeshBlock *ptarget_block = pmy_mesh_->FindMeshBlock(snb.gid);
  // 2) which element in vector of BoundaryVariable *?
  FaceCenteredBoundaryVariable *ptarget_pfbval =
      static_cast<FaceCenteredBoundaryVariable *>(
          ptarget_block->pbval->bvars[bvar_index]);
  Real *target_buf, *send_buf;

  if (is_north) {
    target_buf= ptarget_pfbval->flux_north_recv_[pmy_block_->loc.lx3];
    send_buf = flux_north_send_[polar_block_index];
    BoundaryStatus &target_flag = ptarget_pfbval->flux_north_flag_[pmy_block_->loc.lx3];
    target_flag = BoundaryStatus::arrived;
  } else {
    target_buf= ptarget_pfbval->flux_south_recv_[pmy_block_->loc.lx3];
    send_buf = flux_south_send_[polar_block_index];
    BoundaryStatus &target_flag = ptarget_pfbval->flux_south_flag_[pmy_block_->loc.lx3];
    target_flag = BoundaryStatus::arrived;
  }
  std::memcpy(target_buf, send_buf, ssize*sizeof(Real));
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SendFluxCorrection()
//! \brief Restrict, pack and send the surface EMF to the coarse neighbor(s) if needed

void FaceCenteredBoundaryVariable::SendFluxCorrection() {
  MeshBlock *pmb=pmy_block_;

  // Send non-polar EMF values
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if ((nb.ni.type != NeighborConnect::face) && (nb.ni.type != NeighborConnect::edge))
      break;
    if (bd_var_flcor_.sflag[nb.bufid] == BoundaryStatus::completed) continue;
    int p = 0;
    if (nb.snb.level == pmb->loc.level) {
      if ((nb.ni.type == NeighborConnect::face)
          || ((nb.ni.type == NeighborConnect::edge)
              && (edge_flag_[nb.eid]))) {
        p = LoadFluxBoundaryBufferSameLevel(bd_var_flcor_.send[nb.bufid], nb);
      } else {
        continue;
      }
    } else if (nb.snb.level == pmb->loc.level-1) {
      p = LoadFluxBoundaryBufferToCoarser(bd_var_flcor_.send[nb.bufid], nb);
    } else {
      continue;
    }
    if (nb.snb.rank == Globals::my_rank) { // on the same MPI rank
      CopyFluxCorrectionBufferSameProcess(nb, p);
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&(bd_var_flcor_.req_send[nb.bufid]));
#endif
    bd_var_flcor_.sflag[nb.bufid] = BoundaryStatus::completed;
  }

  // Send polar EMF values
  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    const SimpleNeighborBlock &snb = pbval_->polar_neighbor_north_[n];
    int count = LoadFluxBoundaryBufferToPolar(flux_north_send_[n], snb, true);
    if (snb.rank == Globals::my_rank) { // on the same MPI rank
      CopyPolarBufferSameProcess(snb, count, n, true);
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_flux_north_send_[n]);
#endif
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    const SimpleNeighborBlock &snb = pbval_->polar_neighbor_south_[n];
    int count = LoadFluxBoundaryBufferToPolar(flux_south_send_[n], snb, false);
    if (snb.rank == Globals::my_rank) { // on the same node
      CopyPolarBufferSameProcess(snb, count, n, false);
    }
#ifdef MPI_PARALLEL
    else
      MPI_Start(&req_flux_south_send_[n]);
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetFluxBoundarySameLevel(Real *buf,
//!                                                               const NeighborBlock& nb)
//! \brief Add up the EMF received from a block on the same level
//!        Later they will be divided in the AverageFluxBoundary function

void FaceCenteredBoundaryVariable::SetFluxBoundarySameLevel(Real *buf,
                                                           const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &e1 = pmb->pfield->e.x1e;
  AthenaArray<Real> &e2 = pmb->pfield->e.x2e;
  AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
  int p = 0;
  if (nb.ni.type == NeighborConnect::face) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        // KGF: shearing box
        if (nb.shear) {
          if (nb.fid == BoundaryFace::inner_x1) {
            // store e2 for shearing periodic bcs
            for (int k=pmb->ks; k<=pmb->ke+1; k++) {
              for (int j=pmb->js; j<=pmb->je; j++)
                shear_var_emf_[0].x2e(k,j) += buf[p++];
            }
            // store e3 for shearing periodic bcs
            for (int k=pmb->ks; k<=pmb->ke; k++) {
              for (int j=pmb->js; j<=pmb->je+1; j++)
                shear_var_emf_[0].x3e(k,j) += buf[p++];
            }
          } else if (nb.fid == BoundaryFace::outer_x1) {
            // store e2 for shearing periodic bcs
            for (int k=pmb->ks; k<=pmb->ke+1; k++) {
              for (int j=pmb->js; j<=pmb->je; j++)
                shear_var_emf_[1].x2e(k,j) += buf[p++];
            }
            // store e3 for shearing periodic bcs
            for (int k=pmb->ks; k<=pmb->ke; k++) {
              for (int j=pmb->js; j<=pmb->je+1; j++)
                shear_var_emf_[1].x3e(k,j) += buf[p++];
            }
          }
        } else {
          // unpack e2
          for (int k=pmb->ks; k<=pmb->ke+1; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i) += buf[p++];
          }
          // unpack e3
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je+1; j++)
              e3(k,j,i) += buf[p++];
          }
        } // KGF: shearing box
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        // unpack e1
        Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke+1; k++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i) += sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            e3(k,j,i) += sign*buf[p++];
        }
        // x3 direction
      } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
        int k;
        if (nb.fid == BoundaryFace::inner_x3) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        // unpack e1
        for (int j=pmb->js; j<=pmb->je+1; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i) += buf[p++];
        }
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie+1; i++)
            e2(k,j,i) += buf[p++];
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k = pmb->ks;
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        // KGF: shearing box
        if (pbval_->shearing_box == 1 && nb.shear) {
          if (nb.fid == BoundaryFace::inner_x1) {
            // store e2 for shearing periodic bcs
            for (int j=pmb->js; j<=pmb->je; j++) {
              shear_var_emf_[0].x2e(k+1,j) += buf[p];
              shear_var_emf_[0].x2e(k,j)   += buf[p++];
            }
            // store e3 for shearing periodic bcs
            for (int j=pmb->js; j<=pmb->je+1; j++)
              shear_var_emf_[0].x3e(k,j) += buf[p++];
          } else if (nb.fid == BoundaryFace::outer_x1) {
            // store e2 for shearing periodic bcs
            for (int j=pmb->js; j<=pmb->je; j++) {
              shear_var_emf_[1].x2e(k+1,j) += buf[p];
              shear_var_emf_[1].x2e(k,j)   += buf[p++];
            }
            // store e3 for shearing periodic bcs
            for (int j=pmb->js; j<=pmb->je+1; j++)
              shear_var_emf_[1].x3e(k,j) += buf[p++];
          }
        } else {
          // unpack e2
          for (int j=pmb->js; j<=pmb->je; j++) {
            e2(k+1,j,i) += buf[p];
            e2(k,j,i)   += buf[p++];
          }
          // unpack e3
          for (int j=pmb->js; j<=pmb->je+1; j++)
            e3(k,j,i) += buf[p++];
        } // KGF: shearing box
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        // unpack e1
        for (int i=pmb->is; i<=pmb->ie; i++) {
          e1(k+1,j,i) += buf[p];
          e1(k  ,j,i) += buf[p++];
        }
        // unpack e3
        for (int i=pmb->is; i<=pmb->ie+1; i++)
          e3(k,j,i) += buf[p++];
      }
    } else { // 1D
      int i, j=pmb->js, k = pmb->ks;
      if (nb.fid == BoundaryFace::inner_x1) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      // unpack e2
      e2(k+1,j,i) += buf[p];
      e2(k  ,j,i) += buf[p++];
      // unpack e3
      e3(k,j+1,i) += buf[p];
      e3(k  ,j,i) += buf[p++];
    }
  } else if (nb.ni.type == NeighborConnect::edge) {
    // x1x2 edge (2D and 3D)
    if (nb.eid >= 0 && nb.eid < 4) {
      int i, j;
      if ((nb.eid & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((nb.eid & 2) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      // KGF: shearing box
      if (pbval_->shearing_box == 1 && nb.shear) {
        if (nb.ni.ox1 == -1) {
          // store e3 for shearing periodic bcs
          for (int k = pmb->ks; k<=pmb->ke; k++)
            shear_var_emf_[0].x3e(k,j) += buf[p++];
        } else if (nb.ni.ox1 == 1) {
          // store e3 for shearing periodic bcs
          for (int k = pmb->ks; k<=pmb->ke; k++)
            shear_var_emf_[1].x3e(k,j) += buf[p++];
        }
      } else {
        // unpack e3
        Real sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          e3(k,j,i) += sign*buf[p++];
        }
      } // KGF: shearing box
      // x1x3 edge
    } else if (nb.eid>=4 && nb.eid<8) {
      int i, k;
      if ((nb.eid & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((nb.eid & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      // KGF: shearing box
      if (pbval_->shearing_box == 1 && nb.shear) {
        if (nb.ni.ox1 == -1) {
          // store e2 for shearing periodic bcs
          for (int j=pmb->js; j<=pmb->je; j++)
            shear_var_emf_[0].x2e(k,j) += buf[p++];
        } else if (nb.ni.ox1 == 1) {
          // store e2 for shearing periodic bcs
          for (int j=pmb->js; j<=pmb->je; j++)
            shear_var_emf_[1].x2e(k,j) += buf[p++];
        }
      } else {
        // unpack e2
        for (int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i) += buf[p++];
      } // KGF: shearing box
      // x2x3 edge
    } else if (nb.eid>=8 && nb.eid<12) {
      int j, k;
      if ((nb.eid & 1) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      if ((nb.eid & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      // unpack e1
      Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i) += sign*buf[p++];
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetFluxBoundaryFromFiner(Real *buf,
//!                                                               const NeighborBlock& nb)
//! \brief Add up the EMF received from a block on the finer level
//!        Later they will be divided in the AverageFluxBoundary function

void FaceCenteredBoundaryVariable::SetFluxBoundaryFromFiner(Real *buf,
                                                           const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &e1 = pmb->pfield->e.x1e;
  AthenaArray<Real> &e2 = pmb->pfield->e.x2e;
  AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
  int p = 0;
  if (nb.ni.type == NeighborConnect::face) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i, jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        if (nb.ni.fi1 == 0) {
          ju = pmb->js + pmb->block_size.nx2/2-1;
        } else {
          jl = pmb->js + pmb->block_size.nx2/2;
        }
        if (nb.ni.fi2 == 0) {
          ku = pmb->ks + pmb->block_size.nx3/2-1;
        } else {
          kl = pmb->ks + pmb->block_size.nx3/2;
        }
        // unpack e2
        for (int k=kl; k<=ku+1; k++) {
          for (int j=jl; j<=ju; j++)
            e2(k,j,i) += buf[p++];
        }
        // unpack e3
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju+1; j++)
            e3(k,j,i) += buf[p++];
        }
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j, il = pmb->is, iu = pmb->ie, kl = pmb->ks, ku = pmb->ke;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        if (nb.ni.fi1 == 0) {
          iu = pmb->is + pmb->block_size.nx1/2-1;
        } else {
          il = pmb->is + pmb->block_size.nx1/2;
        }
        if (nb.ni.fi2 == 0) {
          ku = pmb->ks + pmb->block_size.nx3/2-1;
        } else {
          kl = pmb->ks + pmb->block_size.nx3/2;
        }
        // unpack e1
        Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku+1; k++) {
          for (int i=il; i<=iu; i++)
            e1(k,j,i) += sign*buf[p++];
        }
        // unpack e3
        sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k=kl; k<=ku; k++) {
          for (int i=il; i<=iu+1; i++)
            e3(k,j,i) += sign*buf[p++];
        }
        // x3 direction
      } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
        int k, il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je;
        if (nb.fid == BoundaryFace::inner_x3) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        if (nb.ni.fi1 == 0) {
          iu = pmb->is + pmb->block_size.nx1/2-1;
        } else {
          il = pmb->is + pmb->block_size.nx1/2;
        }
        if (nb.ni.fi2 == 0) {
          ju = pmb->js + pmb->block_size.nx2/2-1;
        } else {
          jl = pmb->js + pmb->block_size.nx2/2;
        }
        // unpack e1
        for (int j=jl; j<=ju+1; j++) {
          for (int i=il; i<=iu; i++)
            e1(k,j,i) += buf[p++];
        }
        // unpack e2
        for (int j=jl; j<=ju; j++) {
          for (int i=il; i<=iu+1; i++)
            e2(k,j,i) += buf[p++];
        }
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int k = pmb->ks;
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i, jl = pmb->js, ju = pmb->je;
        if (nb.fid == BoundaryFace::inner_x1) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        if (nb.ni.fi1 == 0) {
          ju = pmb->js + pmb->block_size.nx2/2-1;
        } else {
          jl = pmb->js + pmb->block_size.nx2/2;
        }
        // unpack e2
        for (int j=jl; j<=ju; j++) {
          e2(k+1,j,i) += buf[p];
          e2(k,  j,i) += buf[p++];
        }
        // unpack e3
        for (int j=jl; j<=ju+1; j++)
          e3(k,j,i) += buf[p++];
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j, il = pmb->is, iu = pmb->ie;
        if (nb.fid == BoundaryFace::inner_x2) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        if (nb.ni.fi1 == 0) {
          iu = pmb->is + pmb->block_size.nx1/2-1;
        } else {
          il = pmb->is + pmb->block_size.nx1/2;
        }
        // unpack e1
        for (int i=il; i<=iu; i++) {
          e1(k+1,j,i) += buf[p];
          e1(k  ,j,i) += buf[p++];
        }
        // unpack e3
        for (int i=il; i<=iu+1; i++)
          e3(k,j,i) += buf[p++];
      }
    } else { // 1D
      int i, j = pmb->js, k = pmb->ks;
      if (nb.fid == BoundaryFace::inner_x1) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      // unpack e2
      e2(k+1,j,i) += buf[p];
      e2(k  ,j,i) += buf[p++];
      // unpack e3
      e3(k,j+1,i) += buf[p];
      e3(k  ,j,i) += buf[p++];
    }
  } else if (nb.ni.type == NeighborConnect::edge) {
    if (pmb->block_size.nx3 > 1) { // 3D
      // x1x2 edge
      if (nb.eid>=0 && nb.eid<4) {
        int i, j, kl = pmb->ks, ku = pmb->ke;
        if ((nb.eid & 1) == 0) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        if ((nb.eid & 2) == 0) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        if (nb.ni.fi1 == 0) {
          ku = pmb->ks + pmb->block_size.nx3/2-1;
        } else {
          kl = pmb->ks + pmb->block_size.nx3/2;
        }
        // unpack e3
        Real sign = (nb.polar && flip_across_pole_field[IB3]) ? -1.0 : 1.0;
        for (int k = kl; k<=ku; k++)
          e3(k,j,i) += sign*buf[p++];
        // x1x3 edge
      } else if (nb.eid>=4 && nb.eid<8) {
        int i, k, jl = pmb->js, ju = pmb->je;
        if ((nb.eid & 1) == 0) {
          i = pmb->is;
        } else {
          i = pmb->ie + 1;
        }
        if ((nb.eid & 2) == 0) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        if (nb.ni.fi1 == 0) {
          ju = pmb->js + pmb->block_size.nx2/2-1;
        } else {
          jl = pmb->js + pmb->block_size.nx2/2;
        }
        // unpack e2
        for (int j=jl; j<=ju; j++)
          e2(k,j,i) += buf[p++];
        // x2x3 edge
      } else if (nb.eid>=8 && nb.eid<12) {
        int j, k, il = pmb->is, iu = pmb->ie;
        if ((nb.eid & 1) == 0) {
          j = pmb->js;
        } else {
          j = pmb->je + 1;
        }
        if ((nb.eid & 2) == 0) {
          k = pmb->ks;
        } else {
          k = pmb->ke + 1;
        }
        if (nb.ni.fi1 == 0) {
          iu = pmb->is + pmb->block_size.nx1/2-1;
        } else {
          il = pmb->is + pmb->block_size.nx1/2;
        }
        // unpack e1
        Real sign = (nb.polar && flip_across_pole_field[IB1]) ? -1.0 : 1.0;
        for (int i=il; i<=iu; i++)
          e1(k,j,i) += sign*buf[p++];
      }
    } else if (pmb->block_size.nx2 > 1) { // 2D
      int i, j, k = pmb->ks;
      if ((nb.eid & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((nb.eid & 2) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      // unpack e3
      e3(k,j,i) += buf[p++];
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::SetFluxBoundaryFromPolar(Real **buf_list,
//!                                                           int num_bufs, bool is_north)
//! \brief Overwrite EMF values along polar axis with azimuthal averages

void FaceCenteredBoundaryVariable::SetFluxBoundaryFromPolar(Real **buf_list, int num_bufs,
                                                            bool is_north) {
  MeshBlock *pmb = pmy_block_;
  if (pmb->block_size.nx3 > 1) {
    int j = is_north ? pmb->js : pmb->je + 1;
    int count = 0;
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      Real val = 0.0;
      for (int n = 0; n < num_bufs; ++n)
        val += buf_list[n][count];
      for (int k = pmb->ks-NGHOST; k <= pmb->ke+NGHOST+1; ++k)
        pmb->pfield->e.x1e(k, j, i) = val / num_bufs;
      ++count;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ClearCoarseFluxBoundary()
//! \brief Clear the EMFs on the surface/edge contacting with a finer block

void FaceCenteredBoundaryVariable::ClearCoarseFluxBoundary() {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &e1 = pmb->pfield->e.x1e;
  AthenaArray<Real> &e2 = pmb->pfield->e.x2e;
  AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
  int nl;
  // face
  for (int n=0; n<pbval_->nface_; n++) {
    if (n == BoundaryFace::inner_x1 || n == BoundaryFace::outer_x1) {
      int i;
      if (n == BoundaryFace::inner_x1) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      nl = pbval_->nblevel[1][1][2*n];
      if (nl>pmb->loc.level) { // finer
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i) = 0.0;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i) = 0.0;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i) = e2(pmb->ks+1,j,i) = 0.0;
          for (int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i) = 0.0;
        } else { // 1D
          e2(pmb->ks,pmb->js,i) = e2(pmb->ks+1,pmb->js,i) = 0.0;
          e3(pmb->ks,pmb->js,i) = e3(pmb->ks,pmb->js+1,i) = 0.0;
        }
      }
    }
    if (n == BoundaryFace::inner_x2 || n == BoundaryFace::outer_x2) {
      int j;
      if (n == BoundaryFace::inner_x2) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      nl = pbval_->nblevel[1][2*n-4][1];
      if (nl>pmb->loc.level) { // finer
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i) = 0.0;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i) = 0.0;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i) = e1(pmb->ks+1,j,i) = 0.0;
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i) = 0.0;
        }
      }
    }
    if (n == BoundaryFace::inner_x3 || n == BoundaryFace::outer_x3) {
      int k;
      if (n == BoundaryFace::inner_x3) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      nl = pbval_->nblevel[2*n-8][1][1];
      if (nl>pmb->loc.level) { // finer
        // this is always 3D
        for (int j=pmb->js+1; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i) = 0.0;
        }
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i) = 0.0;
        }
      }
    }
  }
  // edge
  for (int n=0; n<pbval_->nedge_; n++) {
    if (edge_flag_[n]) continue;
    if (n >= 0 && n < 4) {
      int i, j;
      if ((n & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((n & 2) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      for (int k=pmb->ks; k<=pmb->ke; k++)
        e3(k,j,i)=0.0;
      // x1x3 edge
    } else if (n>=4 && n<8) {
      int i, k;
      if ((n & 1) == 0) {
        i = pmb->is;
      } else {
        i = pmb->ie + 1;
      }
      if ((n & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      for (int j=pmb->js; j<=pmb->je; j++)
        e2(k,j,i) = 0.0;
      // x2x3 edge
    } else if (n >= 8 && n < 12) {
      int k, j;
      if ((n & 1) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      if ((n & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i)=0.0;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::AverageFluxBoundary()
//! \brief Set EMF boundary received from a block on the finer level

void FaceCenteredBoundaryVariable::AverageFluxBoundary() {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &e1 = pmb->pfield->e.x1e;
  AthenaArray<Real> &e2 = pmb->pfield->e.x2e;
  AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
  int nl;
  // face
  for (int n=0; n<pbval_->nface_; n++) {
    if ((pbval_->block_bcs[n] != BoundaryFlag::block)
        && (pbval_->block_bcs[n] != BoundaryFlag::periodic)
        && (pbval_->block_bcs[n] != BoundaryFlag::shear_periodic)
        && (pbval_->block_bcs[n] != BoundaryFlag::polar)) continue;
    if (n == BoundaryFace::inner_x1 || n == BoundaryFace::outer_x1) {
      Real div = 0.5;
      int i;
      if (n == BoundaryFace::inner_x1) {
        i = pmb->is;
        if(pbval_->shearing_box == 1 && pbval_->is_shear[0])
          div *= 2.0;
      } else {
        i = pmb->ie + 1;
        if(pbval_->shearing_box == 1 && pbval_->is_shear[1])
          div *= 2.0;
      }
      nl = pbval_->nblevel[1][1][2*n];
      if (nl == pmb->loc.level) { // same ; divide all the face EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++)
              e2(k,j,i) *= div;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js+1; j<=pmb->je; j++)
              e3(k,j,i) *= div;
          }
        } else if (pmb->block_size.nx2 > 1) { // 2D
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(pmb->ks,j,i) *= div, e2(pmb->ks+1,j,i) *= div;
          for (int j=pmb->js+1; j<=pmb->je; j++)
            e3(pmb->ks,j,i) *= div;
        } else { // 1D
          e2(pmb->ks,pmb->js,i) *= 0.5, e2(pmb->ks+1,pmb->js,i) *= 0.5;
          e3(pmb->ks,pmb->js,i) *= 0.5, e3(pmb->ks,pmb->js+1,i) *= 0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          int k = pmb->ks+pmb->block_size.nx3/2;
          for (int j=pmb->js; j<=pmb->je; j++)
            e2(k,j,i) *= 0.5;
        }
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          int j = pmb->js + pmb->block_size.nx2/2;
          for (int k=pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i) *= 0.5;
        }
      }
    }
    if (n == BoundaryFace::inner_x2 || n == BoundaryFace::outer_x2) {
      int j;
      if (n == BoundaryFace::inner_x2) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      nl = pbval_->nblevel[1][2*n-4][1];
      if (nl == pmb->loc.level) { // same ; divide all the face EMFs by 2
        if (pmb->block_size.nx3 > 1) {
          for (int k=pmb->ks+1; k<=pmb->ke; k++) {
            for (int i=pmb->is; i<=pmb->ie; i++)
              e1(k,j,i) *= 0.5;
          }
          for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int i=pmb->is+1; i<=pmb->ie; i++)
              e3(k,j,i) *= 0.5;
          }
        } else if (pmb->block_size.nx2 > 1) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(pmb->ks,j,i) *= 0.5, e1(pmb->ks+1,j,i) *= 0.5;
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e3(pmb->ks,j,i) *= 0.5;
        }
      } else if (nl>pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        if (pmb->block_size.nx3 > 1) { // 3D
          int k = pmb->ks+pmb->block_size.nx3/2;
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i) *= 0.5;
        }
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          int i = pmb->is+pmb->block_size.nx1/2;
          for (int k = pmb->ks; k<=pmb->ke; k++)
            e3(k,j,i) *= 0.5;
        }
      }
    }
    if (n == BoundaryFace::inner_x3 || n == BoundaryFace::outer_x3) {
      int k;
      if (n == BoundaryFace::inner_x3) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      nl = pbval_->nblevel[2*n-8][1][1];
      if (nl == pmb->loc.level) { // same ; divide all the face EMFs by 2
        for (int j=pmb->js+1; j<=pmb->je; j++) {
          for (int i=pmb->is; i<=pmb->ie; i++)
            e1(k,j,i) *= 0.5;
        }
        for (int j=pmb->js; j<=pmb->je; j++) {
          for (int i=pmb->is+1; i<=pmb->ie; i++)
            e2(k,j,i) *= 0.5;
        }
      } else if (nl > pmb->loc.level) { // finer; divide the overlapping EMFs by 2
        // this is always 3D
        int j_fine = pmb->js + pmb->block_size.nx2/2;
        for (int i=pmb->is; i<=pmb->ie; i++)
          e1(k,j_fine,i) *= 0.5;
        int i_fine = pmb->is +pmb->block_size.nx1/2;
        for (int j=pmb->js; j<=pmb->je; j++)
          e2(k,j,i_fine) *= 0.5;
      }
    }
  }
  // edge
  for (int n=0; n<pbval_->nedge_; n++) {
    if (nedge_fine_[n] == 1) continue;
    Real div = 1.0/static_cast<Real>(nedge_fine_[n]);
    Real half_div[2] = {div, div};
    if (pbval_->shearing_box==1) {
      for (int upper=0; upper<2; upper++) {
        if(pbval_->is_shear[upper]) half_div[upper] *= 2.0;
      }
    }
    // x1x2 edge (both 2D and 3D)
    if (n >= 0 && n < 4) {
      int i, j, upper;
      if ((n & 1) == 0) {
        i = pmb->is;
        upper = 0;
      } else {
        i = pmb->ie + 1;
        upper = 1;
      }
      if ((n & 2) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      for (int k=pmb->ks; k<=pmb->ke; k++) {
        e3(k,j,i) *= half_div[upper];
      }
      // x1x3 edge
    } else if (n >= 4 && n < 8) {
      int i, k, upper;
      if ((n & 1) == 0) {
        i = pmb->is;
        upper = 0;
      } else {
        i = pmb->ie + 1;
        upper = 1;
      }
      if ((n & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      for (int j=pmb->js; j<=pmb->je; j++) {
        e2(k,j,i) *= half_div[upper];
      }
      // x2x3 edge
    } else if (n >= 8 && n < 12) {
      int j, k;
      if ((n & 1) == 0) {
        j = pmb->js;
      } else {
        j = pmb->je + 1;
      }
      if ((n & 2) == 0) {
        k = pmb->ks;
      } else {
        k = pmb->ke + 1;
      }
      for (int i=pmb->is; i<=pmb->ie; i++)
        e1(k,j,i) *= div;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::PolarFluxBoundarySingleAzimuthalBlock()
//! \brief polar boundary edge-case:
//!
//! single MeshBlock spans the entire azimuthal (x3) range

void FaceCenteredBoundaryVariable::PolarFluxBoundarySingleAzimuthalBlock() {
  MeshBlock *pmb = pmy_block_;
  if (pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    AthenaArray<Real> &e1 = pmb->pfield->e.x1e;
    AthenaArray<Real> &e3 = pmb->pfield->e.x3e;
    if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
        || pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge) {
      int j = pmb->js;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real tote1 = 0.0;
        for (int k=pmb->ks; k<=pmb->ke; k++)
          tote1 += e1(k,j,i);
        Real e1a = tote1/static_cast<double>(pmb->ke - pmb->ks + 1);
        for (int k=pmb->ks; k<=pmb->ke+1; k++)
          e1(k,j,i) = e1a;
      }
      for (int i=pmb->is; i<=pmb->ie+1; i++) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          pbval_->azimuthal_shift_(k) = e3(k,j,i);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          int k_shift = k;
          k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
          e3(k,j,i) = pbval_->azimuthal_shift_(k_shift);
        }
      }
    }

    if (pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
        || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge) {
      int j = pmb->je + 1;
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real tote1 = 0.0;
        for (int k=pmb->ks; k<=pmb->ke; ++k)
          tote1 += e1(k,j,i);
        Real e1a=tote1/static_cast<double>(pmb->ke - pmb->ks + 1);
        for (int k=pmb->ks; k<=pmb->ke+1; ++k)
          e1(k,j,i) = e1a;
      }
      for (int i=pmb->is; i<=pmb->ie+1; i++) {
        for (int k=pmb->ks; k<=pmb->ke; k++)
          pbval_->azimuthal_shift_(k) = e3(k,j,i);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          int k_shift = k;
          k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
          e3(k,j,i) = pbval_->azimuthal_shift_(k_shift);
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReceiveFluxCorrection()
//! \brief Receive and Apply the surface EMF to the coarse neighbor(s) if needed

bool FaceCenteredBoundaryVariable::ReceiveFluxCorrection() {
  MeshBlock *pmb = pmy_block_;
  bool flag = true;

  // Receive same-level non-polar EMF values
  if (recv_flx_same_lvl_) {
    for (int n=0; n<pbval_->nneighbor; n++) { // first correct the same level
      NeighborBlock& nb = pbval_->neighbor[n];
      if (nb.ni.type != NeighborConnect::face && nb.ni.type != NeighborConnect::edge)
        break;
      if (nb.snb.level!=pmb->loc.level) continue;
      if ((nb.ni.type == NeighborConnect::face)
          || ((nb.ni.type == NeighborConnect::edge) && (edge_flag_[nb.eid]))) {
        if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::completed) continue;
        if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::waiting) {
          if (nb.snb.rank == Globals::my_rank) { // on the same process
            flag = false;
            continue;
          }
#ifdef MPI_PARALLEL
          else { // NOLINT
            int test;
            // probe MPI communications.  This is a bit of black magic that seems to
            // promote communications to top of stack and gets them to complete more
            // quickly
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                       MPI_STATUS_IGNORE);
            MPI_Test(&(bd_var_flcor_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
            if (!static_cast<bool>(test) ) {
              flag = false;
              continue;
            }
            bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::arrived;
          }
#endif
        }
        // boundary arrived; apply EMF correction
        SetFluxBoundarySameLevel(bd_var_flcor_.recv[nb.bufid], nb);
        bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::completed;
      }
    }
    if (!flag) return flag;  // is this flag always false?
    if (pmb->pmy_mesh->multilevel)
      ClearCoarseFluxBoundary();
    recv_flx_same_lvl_ = false;
  }

  // Receive finer non-polar EMF values
  if (pmb->pmy_mesh->multilevel) {
    for (int n=0; n<pbval_->nneighbor; n++) { // then from finer
      NeighborBlock& nb = pbval_->neighbor[n];
      if (nb.ni.type != NeighborConnect::face && nb.ni.type != NeighborConnect::edge)
        break;
      if (nb.snb.level!=pmb->loc.level + 1) continue;
      if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::completed) continue;
      if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::waiting) {
        if (nb.snb.rank == Globals::my_rank) {// on the same process
          flag = false;
          continue;
        }
#ifdef MPI_PARALLEL
        else { // NOLINT
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&(bd_var_flcor_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
          if (!static_cast<bool>(test)) {
            flag = false;
            continue;
          }
          bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::arrived;
        }
#endif
      }
      // boundary arrived; apply EMF correction
      SetFluxBoundaryFromFiner(bd_var_flcor_.recv[nb.bufid], nb);
      bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::completed;
    }
  }

  // Receive polar EMF values
  for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n) {
    const SimpleNeighborBlock &snb = pbval_->polar_neighbor_north_[n];
    if (flux_north_flag_[n]  ==  BoundaryStatus::waiting) {
      if (snb.rank  ==  Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int recv_flag;
        MPI_Test(&req_flux_north_recv_[n], &recv_flag, MPI_STATUS_IGNORE);
        if (!recv_flag) {
          flag = false;
          continue;
        }
        flux_north_flag_[n] = BoundaryStatus::arrived;
      }
#endif
    }
  }
  for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n) {
    const SimpleNeighborBlock &snb = pbval_->polar_neighbor_south_[n];
    if (flux_south_flag_[n] == BoundaryStatus::waiting) {
      if (snb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // NOLINT
        int recv_flag;
        MPI_Test(&req_flux_south_recv_[n], &recv_flag, MPI_STATUS_IGNORE);
        if (!recv_flag) {
          flag = false;
          continue;
        }
        flux_south_flag_[n] = BoundaryStatus::arrived;
      }
#endif
    }
  }

  if (flag) {
    AverageFluxBoundary();
    if (pbval_->num_north_polar_blocks_ > 0)
      SetFluxBoundaryFromPolar(flux_north_recv_, pbval_->num_north_polar_blocks_, true);
    for (int n = 0; n < pbval_->num_north_polar_blocks_; ++n)
      flux_north_flag_[n] = BoundaryStatus::completed;
    if (pbval_->num_south_polar_blocks_ > 0)
      SetFluxBoundaryFromPolar(flux_south_recv_, pbval_->num_south_polar_blocks_, false);
    for (int n = 0; n < pbval_->num_south_polar_blocks_; ++n)
      flux_south_flag_[n] = BoundaryStatus::completed;
    if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
        || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
        || pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge
        || pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge)
      PolarFluxBoundarySingleAzimuthalBlock();
  }
  return flag;
}
