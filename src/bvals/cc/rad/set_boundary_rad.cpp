//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polar_rad.cpp
//  \brief implementation of boundary setting for radiation

// Athena++ headers
#include "bvals_rad.hpp"
#include "../../../athena.hpp"              // Real, indices
#include "../../../athena_arrays.hpp"       // AthenaArray
#include "../../../mesh/mesh.hpp"           // MeshBlock
#include "../../../utils/buffer_utils.hpp"  // BufferUtility

//----------------------------------------------------------------------------------------
// Radiation boundary for blocks at the same refinement level
// Inputs:
//   buf: buffer containing values to be placed in ghost zone
//   nb: neighbor metadata
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in buffer.
//   In case of polar boundary, properly reflects radiation in x2- and x3-directions,
//       needed when going from active cells in a block on one side of the pole to ghost
//       cells in a block on the other side of the pole.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.
//   The special case of a single block wrapping around the pole requires an additional
//       shift of cells in the azimuthal direction, which is done elsewhere (like with all
//       other variables).

void RadBoundaryVariable::SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) {

  // Calculate indices
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  if (nb.ni.ox1 == 0) {
    si = pmb->is;
    ei = pmb->ie;
  } else if (nb.ni.ox1 > 0) {
    si = pmb->ie + 1;
    ei = pmb->ie + NGHOST;
  } else {
    si = pmb->is - NGHOST;
    ei = pmb->is - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->js;
    ej = pmb->je;
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->je + 1;
    ej = pmb->je + NGHOST;
  } else {
    sj = pmb->js - NGHOST;
    ej = pmb->js - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->ks;
    ek = pmb->ke;
  } else if (nb.ni.ox3 > 0) {
    sk = pmb->ke + 1;
    ek = pmb->ke + NGHOST;
  } else {
    sk = pmb->ks - NGHOST;
    ek = pmb->ks - 1;
  }

  // Prepare buffer index
  int p = 0;

  // Set non-polar boundary
  if (not nb.polar) {
    BufferUtility::UnpackData(buf, *var_cc, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // Set polar boundary
  } else {

    // Set polar boundary without accounting for reflections
    for (int n = nl_; n <= nu_; ++n) {
      for (int k = sk; k <= ek; ++k) {
        for (int j = ej; j >= sj; --j) {
          for (int i = si; i <= ei; ++i) {
            (*var_cc)(n,k,j,i) = buf[p++];
          }
        }
      }
    }

    // Account for reflections at north polar boundary
    if (nb.ni.ox2 < 0) {
      for (int k = sk; k <= ek; ++k) {
        for (int dj = 0; dj < NGHOST; ++dj) {
          int j = sj + dj;
          for (int i = si; i <= ei; ++i) {
            for (int n = nl_; n <= nu_; ++n) {
              polar_vals_(n) = (*var_cc)(n,k,j,i);
            }
            for (int n = nl_; n <= nu_; ++n) {
              int ind_a = polar_ind_north_(0,n,k,dj,i);
              int ind_b = polar_ind_north_(1,n,k,dj,i);
              int ind_c = polar_ind_north_(2,n,k,dj,i);
              int ind_d = polar_ind_north_(3,n,k,dj,i);
              Real frac_a = polar_frac_north_(0,n,k,dj,i);
              Real frac_b = polar_frac_north_(1,n,k,dj,i);
              Real frac_c = polar_frac_north_(2,n,k,dj,i);
              Real frac_d = polar_frac_north_(3,n,k,dj,i);
              Real val_a = polar_vals_(ind_a);
              Real val_b = polar_vals_(ind_b);
              Real val_c = polar_vals_(ind_c);
              Real val_d = polar_vals_(ind_d);
              (*var_cc)(n,k,j,i) =
                  frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
            }
          }
        }
      }

    // Account for reflections at south polar boundary
    } else {
      for (int k = sk; k <= ek; ++k) {
        for (int dj = 0; dj < NGHOST; ++dj) {
          int j = sj + dj;
          for (int i = si; i <= ei; ++i) {
            for (int n = nl_; n <= nu_; ++n) {
              polar_vals_(n) = (*var_cc)(n,k,j,i);
            }
            for (int n = nl_; n <= nu_; ++n) {
              int ind_a = polar_ind_south_(0,n,k,dj,i);
              int ind_b = polar_ind_south_(1,n,k,dj,i);
              int ind_c = polar_ind_south_(2,n,k,dj,i);
              int ind_d = polar_ind_south_(3,n,k,dj,i);
              Real frac_a = polar_frac_south_(0,n,k,dj,i);
              Real frac_b = polar_frac_south_(1,n,k,dj,i);
              Real frac_c = polar_frac_south_(2,n,k,dj,i);
              Real frac_d = polar_frac_south_(3,n,k,dj,i);
              Real val_a = polar_vals_(ind_a);
              Real val_b = polar_vals_(ind_b);
              Real val_c = polar_vals_(ind_c);
              Real val_d = polar_vals_(ind_d);
              (*var_cc)(n,k,j,i) =
                  frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation boundary for blocks at a coarser refinement level
// Inputs:
//   buf: buffer containing values to be placed in ghost zone
//   nb: neighbor metadata
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in buffer.
//   In case of polar boundary, properly reflects radiation in x2- and x3-directions,
//       needed when going from active cells in a block on one side of the pole to ghost
//       cells in a block on the other side of the pole.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.
//   The special case of a single block wrapping around the pole requires an additional
//       shift of cells in the azimuthal direction, which is done elsewhere (like with all
//       other variables).

void RadBoundaryVariable::SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) {

  // Calculate indices
  MeshBlock *pmb = pmy_block_;
  int cng = pmb->cnghost;
  int si, sj, sk, ei, ej, ek;
  if (nb.ni.ox1 == 0) {
    si = pmb->cis;
    ei = pmb->cie;
    if ((pmb->loc.lx1 & 1LL) == 0LL) {
      ei += cng;
    } else {
      si -= cng;
    }
  } else if (nb.ni.ox1 > 0) {
    si = pmb->cie + 1;
    ei = pmb->cie + cng;
  } else {
    si = pmb->cis - cng;
    ei = pmb->cis - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->cjs;
    ej = pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL) == 0LL) {
        ej += cng;
      } else {
        sj -= cng;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->cje + 1;
    ej = pmb->cje + cng;
  } else {
    sj = pmb->cjs - cng;
    ej = pmb->cjs - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->cks;
    ek = pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL) == 0LL) {
        ek += cng;
      } else {
        sk -= cng;
      }
    }
  } else if (nb.ni.ox3 > 0) {
    sk = pmb->cke + 1;
    ek = pmb->cke + cng;
  } else {
    sk = pmb->cks - cng;
    ek = pmb->cks - 1;
  }

  // Prepare buffer index
  int p = 0;

  // Set non-polar boundary
  if (not nb.polar) {
    BufferUtility::UnpackData(buf, *coarse_buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // Set polar boundary
  } else {

    // Set polar boundary without accounting for reflections
    for (int n = nl_; n <= nu_; ++n) {
      for (int k = sk; k <= ek; ++k) {
        for (int j = ej; j >= sj; --j) {
          for (int i = si; i <= ei; ++i) {
            (*coarse_buf)(n,k,j,i) = buf[p++];
          }
        }
      }
    }

    // Account for reflections at north polar boundary
    if (nb.ni.ox2 < 0) {
      for (int k = sk; k <= ek; ++k) {
        for (int dj = 0; dj < NGHOST; ++dj) {
          int j = sj + dj;
          for (int i = si; i <= ei; ++i) {
            for (int n = nl_; n <= nu_; ++n) {
              polar_vals_(n) = (*coarse_buf)(n,k,j,i);
            }
            for (int n = nl_; n <= nu_; ++n) {
              int ind_a = polar_ind_north_(0,n,k,dj,i);
              int ind_b = polar_ind_north_(1,n,k,dj,i);
              int ind_c = polar_ind_north_(2,n,k,dj,i);
              int ind_d = polar_ind_north_(3,n,k,dj,i);
              Real frac_a = polar_frac_north_(0,n,k,dj,i);
              Real frac_b = polar_frac_north_(1,n,k,dj,i);
              Real frac_c = polar_frac_north_(2,n,k,dj,i);
              Real frac_d = polar_frac_north_(3,n,k,dj,i);
              Real val_a = polar_vals_(ind_a);
              Real val_b = polar_vals_(ind_b);
              Real val_c = polar_vals_(ind_c);
              Real val_d = polar_vals_(ind_d);
              (*coarse_buf)(n,k,j,i) =
                  frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
            }
          }
        }
      }

    // Account for reflections at south polar boundary
    } else {
      for (int k = sk; k <= ek; ++k) {
        for (int dj = 0; dj < NGHOST; ++dj) {
          int j = sj + dj;
          for (int i = si; i <= ei; ++i) {
            for (int n = nl_; n <= nu_; ++n) {
              polar_vals_(n) = (*coarse_buf)(n,k,j,i);
            }
            for (int n = nl_; n <= nu_; ++n) {
              int ind_a = polar_ind_south_(0,n,k,dj,i);
              int ind_b = polar_ind_south_(1,n,k,dj,i);
              int ind_c = polar_ind_south_(2,n,k,dj,i);
              int ind_d = polar_ind_south_(3,n,k,dj,i);
              Real frac_a = polar_frac_south_(0,n,k,dj,i);
              Real frac_b = polar_frac_south_(1,n,k,dj,i);
              Real frac_c = polar_frac_south_(2,n,k,dj,i);
              Real frac_d = polar_frac_south_(3,n,k,dj,i);
              Real val_a = polar_vals_(ind_a);
              Real val_b = polar_vals_(ind_b);
              Real val_c = polar_vals_(ind_c);
              Real val_d = polar_vals_(ind_d);
              (*coarse_buf)(n,k,j,i) =
                  frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation boundary for blocks at a finer refinement level
// Inputs:
//   buf: buffer containing values to be placed in ghost zone
//   nb: neighbor metadata
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in buffer.
//   In case of polar boundary, properly reflects radiation in x2- and x3-directions,
//       needed when going from active cells in a block on one side of the pole to ghost
//       cells in a block on the other side of the pole.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.
//   The special case of a single block wrapping around the pole requires an additional
//       shift of cells in the azimuthal direction, which is done elsewhere (like with all
//       other variables).

void RadBoundaryVariable::SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) {

  // Calculate indices
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  if (nb.ni.ox1 == 0) {
    si = pmb->is;
    ei = pmb->ie;
    if (nb.ni.fi1 == 1) {
      si += pmb->block_size.nx1/2;
    } else {
      ei -= pmb->block_size.nx1/2;
    }
  } else if (nb.ni.ox1 > 0) {
    si = pmb->ie + 1;
    ei = pmb->ie + NGHOST;
  } else {
    si = pmb->is - NGHOST;
    ei = pmb->is - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->js;
    ej = pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ni.ox1 != 0) {
        if (nb.ni.fi1 == 1) {
          sj += pmb->block_size.nx2/2;
        } else {
          ej -= pmb->block_size.nx2/2;
        }
      } else {
        if (nb.ni.fi2 == 1) {
          sj += pmb->block_size.nx2/2;
        } else {
          ej -= pmb->block_size.nx2/2;
        }
      }
    }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->je + 1;
    ej = pmb->je + NGHOST;
  } else {
    sj = pmb->js - NGHOST;
    ej = pmb->js - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->ks;
    ek = pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ni.ox1 != 0 and nb.ni.ox2 != 0) {
        if (nb.ni.fi1 == 1) {
          sk += pmb->block_size.nx3/2;
        } else {
          ek -= pmb->block_size.nx3/2;
        }
      } else {
        if (nb.ni.fi2 == 1) {
          sk += pmb->block_size.nx3/2;
        } else {
          ek -= pmb->block_size.nx3/2;
        }
      }
    }
  } else if (nb.ni.ox3 > 0) {
    sk = pmb->ke + 1;
    ek = pmb->ke + NGHOST;
  } else {
    sk = pmb->ks - NGHOST;
    ek = pmb->ks - 1;
  }

  // Prepare buffer index
  int p = 0;

  // Set non-polar boundary
  if (not nb.polar) {
    BufferUtility::UnpackData(buf, *var_cc, nl_, nu_, si, ei, sj, ej, sk, ek, p);

  // Set polar boundary
  } else {

    // Set polar boundary without accounting for reflections
    for (int n = nl_; n <= nu_; ++n) {
      for (int k = sk; k <= ek; ++k) {
        for (int j = ej; j >= sj; --j) {
          for (int i = si; i <= ei; ++i) {
            (*var_cc)(n,k,j,i) = buf[p++];
          }
        }
      }
    }

    // Account for reflections at north polar boundary
    if (nb.ni.ox2 < 0) {
      for (int k = sk; k <= ek; ++k) {
        for (int dj = 0; dj < NGHOST; ++dj) {
          int j = sj + dj;
          for (int i = si; i <= ei; ++i) {
            for (int n = nl_; n <= nu_; ++n) {
              polar_vals_(n) = (*var_cc)(n,k,j,i);
            }
            for (int n = nl_; n <= nu_; ++n) {
              int ind_a = polar_ind_north_(0,n,k,dj,i);
              int ind_b = polar_ind_north_(1,n,k,dj,i);
              int ind_c = polar_ind_north_(2,n,k,dj,i);
              int ind_d = polar_ind_north_(3,n,k,dj,i);
              Real frac_a = polar_frac_north_(0,n,k,dj,i);
              Real frac_b = polar_frac_north_(1,n,k,dj,i);
              Real frac_c = polar_frac_north_(2,n,k,dj,i);
              Real frac_d = polar_frac_north_(3,n,k,dj,i);
              Real val_a = polar_vals_(ind_a);
              Real val_b = polar_vals_(ind_b);
              Real val_c = polar_vals_(ind_c);
              Real val_d = polar_vals_(ind_d);
              (*var_cc)(n,k,j,i) =
                  frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
            }
          }
        }
      }

    // Account for reflections at south polar boundary
    } else {
      for (int k = sk; k <= ek; ++k) {
        for (int dj = 0; dj < NGHOST; ++dj) {
          int j = sj + dj;
          for (int i = si; i <= ei; ++i) {
            for (int n = nl_; n <= nu_; ++n) {
              polar_vals_(n) = (*var_cc)(n,k,j,i);
            }
            for (int n = nl_; n <= nu_; ++n) {
              int ind_a = polar_ind_south_(0,n,k,dj,i);
              int ind_b = polar_ind_south_(1,n,k,dj,i);
              int ind_c = polar_ind_south_(2,n,k,dj,i);
              int ind_d = polar_ind_south_(3,n,k,dj,i);
              Real frac_a = polar_frac_south_(0,n,k,dj,i);
              Real frac_b = polar_frac_south_(1,n,k,dj,i);
              Real frac_c = polar_frac_south_(2,n,k,dj,i);
              Real frac_d = polar_frac_south_(3,n,k,dj,i);
              Real val_a = polar_vals_(ind_a);
              Real val_b = polar_vals_(ind_b);
              Real val_c = polar_vals_(ind_c);
              Real val_d = polar_vals_(ind_d);
              (*var_cc)(n,k,j,i) =
                  frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
            }
          }
        }
      }
    }
  }
  return;
}

