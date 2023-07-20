//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weighted_ave.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>  // std::binary_search
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <vector>     // std::vector

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "mesh.hpp"

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::WeightedAve(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
//!                                 AthenaArray<Real> &u_in2, AthenaArray<Real> &u_in3,
//!                                 AthenaArray<Real> &u_in4, const Real wght[5])
//! \brief Compute weighted average of AthenaArrays (including cell-averaged U in time
//!        integrator step)
//!
//! * consider every possible simplified form of weighted sum operator:
//!   U = a*U + b*U1 + c*U2 + d*U3 + e*U4
//! * assuming all 3x arrays are of the same size (or at least u_out is equal or larger
//!   than each input array) in each array dimension, and full range is desired:
//!   nx4*(3D real MeshBlock cells)

void MeshBlock::WeightedAve(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
                            AthenaArray<Real> &u_in2, AthenaArray<Real> &u_in3,
                            AthenaArray<Real> &u_in4, const Real wght[5]) {
  const int nu = u_out.GetDim4() - 1;

  // u_in2, u_in3, and/or u_in4 may be unallocated AthenaArrays if using a
  // 2S time integrator without STS
  if (wght[0] == 1.0) {
    if (wght[4] != 0.0) {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) += wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i)
                                + wght[3]*u_in3(n,k,j,i) + wght[4]*u_in4(n,k,j,i);
            }
          }
        }
      }
    } else { // do not dereference u_in4
      if (wght[3] != 0.0) {
        for (int n=0; n<=nu; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) += wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i)
                                  + wght[3]*u_in3(n,k,j,i);
              }
            }
          }
        }
      } else { // do not dereference u_in3
        if (wght[2] != 0.0) {
          for (int n=0; n<=nu; ++n) {
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  u_out(n,k,j,i) += wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i);
                }
              }
            }
          }
        } else { // do not dereference u_in2
          if (wght[1] != 0.0) {
            for (int n=0; n<=nu; ++n) {
              for (int k=ks; k<=ke; ++k) {
                for (int j=js; j<=je; ++j) {
#pragma omp simd
                  for (int i=is; i<=ie; ++i) {
                    u_out(n,k,j,i) += wght[1]*u_in1(n,k,j,i);
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (wght[0] == 0.0) {
    if (wght[4] != 0.0) {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i)
                               + wght[3]*u_in3(n,k,j,i) + wght[4]*u_in4(n,k,j,i);
            }
          }
        }
      }
    } else {
      if (wght[3] != 0.0) {
        for (int n=0; n<=nu; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) = wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i)
                                 + wght[3]*u_in3(n,k,j,i);
              }
            }
          }
        }
      } else {
        if (wght[2] != 0.0) {
          for (int n=0; n<=nu; ++n) {
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  u_out(n,k,j,i) = wght[1]*u_in1(n,k,j,i) + wght[2]*u_in2(n,k,j,i);
                }
              }
            }
          }
        } else {
          if (wght[1] == 1.0) {
            // just deep copy
            for (int n=0; n<=nu; ++n) {
              for (int k=ks; k<=ke; ++k) {
                for (int j=js; j<=je; ++j) {
#pragma omp simd
                  for (int i=is; i<=ie; ++i) {
                    u_out(n,k,j,i) = u_in1(n,k,j,i);
                  }
                }
              }
            }
          } else {
            for (int n=0; n<=nu; ++n) {
              for (int k=ks; k<=ke; ++k) {
                for (int j=js; j<=je; ++j) {
#pragma omp simd
                  for (int i=is; i<=ie; ++i) {
                    u_out(n,k,j,i) = wght[1]*u_in1(n,k,j,i);
                  }
                }
              }
            }
          }
        }
      }
    }
  } else {
    if (wght[4] != 0.0) {
      for (int n=0; n<=nu; ++n) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i)
                               + wght[2]*u_in2(n,k,j,i) + wght[3]*u_in3(n,k,j,i)
                               + wght[4]*u_in4(n,k,j,i);
            }
          }
        }
      }
    } else { // do not dereference u_in4
      if (wght[3] != 0.0) {
        for (int n=0; n<=nu; ++n) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i)
                                 + wght[2]*u_in2(n,k,j,i) + wght[3]*u_in3(n,k,j,i);
              }
            }
          }
        }
      } else { // do not dereference u_in3
        if (wght[2] != 0.0) {
          for (int n=0; n<=nu; ++n) {
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i)
                               + wght[2]*u_in2(n,k,j,i);
                }
              }
            }
          }
        } else { // do not dereference u_in2
          if (wght[1] != 0.0) {
            for (int n=0; n<=nu; ++n) {
              for (int k=ks; k<=ke; ++k) {
                for (int j=js; j<=je; ++j) {
#pragma omp simd
                  for (int i=is; i<=ie; ++i) {
                    u_out(n,k,j,i) = wght[0]*u_out(n,k,j,i) + wght[1]*u_in1(n,k,j,i);
                  }
                }
              }
            }
          } else { // do not dereference u_in1
            for (int n=0; n<=nu; ++n) {
              for (int k=ks; k<=ke; ++k) {
                for (int j=js; j<=je; ++j) {
#pragma omp simd
                  for (int i=is; i<=ie; ++i) {
                    u_out(n,k,j,i) *= wght[0];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}



//override function for arrays with different order
void MeshBlock::WeightedAve(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
                            AthenaArray<Real> &u_in2, const Real wght[3], int flag) {
  // consider every possible simplified form of weighted sum operator:
  // U = a*U + b*U1 + c*U2

  // assuming all 3x arrays are of the same size (or at least u_out is equal or larger
  // than each input array) in each array dimension, and full range is desired:
  // nx4*(3D real MeshBlock cells)

  if (flag == 1) {
    const int nu = u_out.GetDim1() - 1;
    // u_in2 may be an unallocated AthenaArray if using a 2S time integrator
    if (wght[0] == 1.0) {
      if (wght[2] != 0.0) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
#pragma omp simd
              for (int n=0; n<=nu; ++n) {
                u_out(k,j,i,n) += wght[1]*u_in1(k,j,i,n) + wght[2]*u_in2(k,j,i,n);
              }
            }
          }
        }
      } else { // do not dereference u_in2
        if (wght[1] != 0.0) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie; ++i) {
#pragma omp simd
                for (int n=0; n<=nu; ++n) {
                  u_out(k,j,i,n) += wght[1]*u_in1(k,j,i,n);
                }
              }
            }
          }
        }
      }
    } else if (wght[0] == 0.0) {
      if (wght[2] != 0.0) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
#pragma omp simd
              for (int n=0; n<=nu; ++n) {
                u_out(k,j,i,n) = wght[1]*u_in1(k,j,i,n) + wght[2]*u_in2(k,j,i,n);
              }
            }
          }
        }
      } else if (wght[1] == 1.0) {
        // just deep copy
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
#pragma omp simd
              for (int n=0; n<=nu; ++n) {
                u_out(k,j,i,n) = u_in1(k,j,i,n);
              }
            }
          }
        }
      } else {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
#pragma omp simd
              for (int n=0; n<=nu; ++n) {
                u_out(k,j,i,n) = wght[1]*u_in1(k,j,i,n);
              }
            }
          }
        }
      }
    } else {
      if (wght[2] != 0.0) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
#pragma omp simd
              for (int n=0; n<=nu; ++n) {
                u_out(k,j,i,n) = wght[0]*u_out(k,j,i,n) + wght[1]*u_in1(k,j,i,n)
                                 + wght[2]*u_in2(k,j,i,n);
              }
            }
          }
        }
      } else { // do not dereference u_in2
        if (wght[1] != 0.0) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie; ++i) {
#pragma omp simd
                for (int n=0; n<=nu; ++n) {
                  u_out(k,j,i,n) = wght[0]*u_out(k,j,i,n) + wght[1]*u_in1(k,j,i,n);
                }
              }
            }
          }
        } else { // do not dereference u_in1
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie; ++i) {
#pragma omp simd
                for (int n=0; n<=nu; ++n) {
                  u_out(k,j,i,n) *= wght[0];
                }
              }
            }
          }
        }
      }
    }
  } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in MeshBlock::WeightedAve" << std::endl
          << "flag=" << flag << " not supported!."
          << std::endl;
      ATHENA_ERROR(msg);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::WeightedAve(FaceField &b_out, FaceField &b_in1,
//!                                 FaceField &b_in2, FaceField &b_in3,
//!                                 FaceField &b_in4, const Real wght[5])
//! \brief Compute weighted average of face-averaged B in time integrator step

void MeshBlock::WeightedAve(FaceField &b_out, FaceField &b_in1,
                            FaceField &b_in2, FaceField &b_in3,
                            FaceField &b_in4, const Real wght[5]) {
  int jl=js; int ju=je+1;
  // move these limit modifications outside the loop
  if (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
      || pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge)
    jl=js+1;
  if (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
      || pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge)
    ju=je;

  // Note: these loops can be combined now that they avoid curl terms
  // Only need to separately account for the final longitudinal face in each loop limit
  if (wght[0] == 1.0) {
    if (wght[4] != 0.0) {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) += wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i)
                                + wght[3]*b_in3.x1f(k,j,i) + wght[4]*b_in4.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) += wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i)
                                + wght[3]*b_in3.x3f(k,j,i) + wght[4]*b_in4.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) += wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i)
                                + wght[3]*b_in3.x3f(k,j,i) + wght[4]*b_in4.x3f(k,j,i);
          }
        }
      }
    } else { // do not dereference u_in4
      if (wght[3] != 0.0) {
      //---- B1
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie+1; ++i) {
              b_out.x1f(k,j,i) += wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i)
                                  + wght[3]*b_in3.x1f(k,j,i);
            }
          }
        }
        //---- B2
        for (int k=ks; k<=ke; ++k) {
          for (int j=jl; j<=ju; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x2f(k,j,i) += wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i)
                                   + wght[3]*b_in3.x2f(k,j,i);
            }
          }
        }
        //---- B3
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x3f(k,j,i) += wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i)
                                  + wght[3]*b_in3.x3f(k,j,i);
            }
          }
        }
      } else { // do not dereference u_in3
        if (wght[2] != 0.0) {
          //---- B1
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie+1; ++i) {
                b_out.x1f(k,j,i) += wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i);
              }
            }
          }
          //---- B2
          for (int k=ks; k<=ke; ++k) {
            for (int j=jl; j<=ju; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                b_out.x2f(k,j,i) += wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i);
              }
            }
          }
          //---- B3
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                b_out.x3f(k,j,i) += wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i);
              }
            }
          }
        } else { // do not dereference u_in2
          if (wght[1] != 0.0) {
            //---- B1
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie+1; ++i) {
                  b_out.x1f(k,j,i) += wght[1]*b_in1.x1f(k,j,i);
                }
              }
            }
            //---- B2
            for (int k=ks; k<=ke; ++k) {
              for (int j=jl; j<=ju; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x2f(k,j,i) += wght[1]*b_in1.x2f(k,j,i);
                }
              }
            }
            //---- B3
            for (int k=ks; k<=ke+1; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x3f(k,j,i) += wght[1]*b_in1.x3f(k,j,i);
                }
              }
            }
          }
        }
      }
    }
  } else if (wght[0] == 0.0) {
    if (wght[4] != 0.0) {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) = wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i)
                               + wght[3]*b_in3.x1f(k,j,i) + wght[4]*b_in4.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) = wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i)
                               + wght[3]*b_in3.x2f(k,j,i) + wght[4]*b_in4.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) = wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i)
                               + wght[3]*b_in3.x3f(k,j,i) + wght[4]*b_in4.x3f(k,j,i);
          }
        }
      }
    } else {
      if (wght[3] != 0.0) {
        //---- B1
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie+1; ++i) {
              b_out.x1f(k,j,i) = wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i)
                                 + wght[3]*b_in3.x1f(k,j,i);
            }
          }
        }
        //---- B2
        for (int k=ks; k<=ke; ++k) {
          for (int j=jl; j<=ju; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x2f(k,j,i) = wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i)
                                 + wght[3]*b_in3.x2f(k,j,i);
            }
          }
        }
        //---- B3
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x3f(k,j,i) = wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i)
                                 + wght[3]*b_in3.x3f(k,j,i);
            }
          }
        }
      } else {
        if (wght[2] != 0.0) {
          //---- B1
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie+1; ++i) {
                b_out.x1f(k,j,i) = wght[1]*b_in1.x1f(k,j,i) + wght[2]*b_in2.x1f(k,j,i);
              }
            }
          }
          //---- B2
          for (int k=ks; k<=ke; ++k) {
            for (int j=jl; j<=ju; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                b_out.x2f(k,j,i) = wght[1]*b_in1.x2f(k,j,i) + wght[2]*b_in2.x2f(k,j,i);
              }
            }
          }
          //---- B3
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                b_out.x3f(k,j,i) = wght[1]*b_in1.x3f(k,j,i) + wght[2]*b_in2.x3f(k,j,i);
              }
            }
          }
        } else {
          if (wght[1] == 1.0) {
            // just deep copy
            //---- B1
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie+1; ++i) {
                  b_out.x1f(k,j,i) = b_in1.x1f(k,j,i);
                }
              }
            }
            //---- B2
            for (int k=ks; k<=ke; ++k) {
              for (int j=jl; j<=ju; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x2f(k,j,i) = b_in1.x2f(k,j,i);
                }
              }
            }
            //---- B3
            for (int k=ks; k<=ke+1; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x3f(k,j,i) = b_in1.x3f(k,j,i);
                }
              }
            }
          } else {
            //---- B1
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie+1; ++i) {
                  b_out.x1f(k,j,i) = wght[1]*b_in1.x1f(k,j,i);
                }
              }
            }
            //---- B2
            for (int k=ks; k<=ke; ++k) {
              for (int j=jl; j<=ju; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x2f(k,j,i) = wght[1]*b_in1.x2f(k,j,i);
                }
              }
            }
            //---- B3
            for (int k=ks; k<=ke+1; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x3f(k,j,i) = wght[1]*b_in1.x3f(k,j,i);
                }
              }
            }
          }
        }
      }
    }
  } else {
    if (wght[4] != 0.0) {
      //---- B1
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i)
                               + wght[2]*b_in2.x1f(k,j,i) + wght[3]*b_in3.x1f(k,j,i)
                               + wght[4]*b_in4.x1f(k,j,i);
          }
        }
      }
      //---- B2
      for (int k=ks; k<=ke; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i)
                               + wght[2]*b_in2.x2f(k,j,i) + wght[3]*b_in3.x2f(k,j,i)
                               + wght[4]*b_in4.x2f(k,j,i);
          }
        }
      }
      //---- B3
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i)
                               + wght[2]*b_in2.x3f(k,j,i) + wght[3]*b_in3.x3f(k,j,i)
                               + wght[4]*b_in4.x3f(k,j,i);
          }
        }
      }
    } else { // do not dereference u_in4
      if (wght[3] != 0.0) {
        //---- B1
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie+1; ++i) {
              b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i)
                                 + wght[2]*b_in2.x1f(k,j,i) + wght[3]*b_in3.x1f(k,j,i);
            }
          }
        }
        //---- B2
        for (int k=ks; k<=ke; ++k) {
          for (int j=jl; j<=ju; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i)
                                 + wght[2]*b_in2.x2f(k,j,i) + wght[3]*b_in3.x2f(k,j,i);
            }
          }
        }
        //---- B3
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i)
                                 + wght[2]*b_in2.x3f(k,j,i) + wght[3]*b_in3.x3f(k,j,i);
            }
          }
        }
      } else { // do not dereference u_in3
        if (wght[2] != 0.0) {
          //---- B1
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie+1; ++i) {
                b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i)
                                   + wght[2]*b_in2.x1f(k,j,i);
              }
            }
          }
          //---- B2
          for (int k=ks; k<=ke; ++k) {
            for (int j=jl; j<=ju; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i)
                                   + wght[2]*b_in2.x2f(k,j,i);
              }
            }
          }
          //---- B3
          for (int k=ks; k<=ke+1; ++k) {
            for (int j=js; j<=je; ++j) {
#pragma omp simd
              for (int i=is; i<=ie; ++i) {
                b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i)
                                   + wght[2]*b_in2.x3f(k,j,i);
              }
            }
          }
        } else { // do not dereference u_in2
          if (wght[1] != 0.0) {
            //---- B1
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie+1; ++i) {
                  b_out.x1f(k,j,i) = wght[0]*b_out.x1f(k,j,i) + wght[1]*b_in1.x1f(k,j,i);
                }
              }
            }
            //---- B2
            for (int k=ks; k<=ke; ++k) {
              for (int j=jl; j<=ju; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x2f(k,j,i) = wght[0]*b_out.x2f(k,j,i) + wght[1]*b_in1.x2f(k,j,i);
                }
              }
            }
            //---- B3
            for (int k=ks; k<=ke+1; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x3f(k,j,i) = wght[0]*b_out.x3f(k,j,i) + wght[1]*b_in1.x3f(k,j,i);
                }
              }
            }
          } else { // do not dereference u_in1
            //---- B1
            for (int k=ks; k<=ke; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie+1; ++i) {
                  b_out.x1f(k,j,i) *= wght[0];
                }
              }
            }
            //---- B2
            for (int k=ks; k<=ke; ++k) {
              for (int j=jl; j<=ju; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x2f(k,j,i) *= wght[0];
                }
              }
            }
            //---- B3
            for (int k=ks; k<=ke+1; ++k) {
              for (int j=js; j<=je; ++j) {
#pragma omp simd
                for (int i=is; i<=ie; ++i) {
                  b_out.x3f(k,j,i) *= wght[0];
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}
