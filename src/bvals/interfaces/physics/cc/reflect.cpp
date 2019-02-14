//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX1(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void CellCenteredBoundaryVariable::ReflectInnerX1(
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=nl; i<=nu; ++i) {
          prim(IVX,k,j,il-i) = -prim(IVX,k,j,(il+i-1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=nl; i<=nu; ++i) {
          prim(n,k,j,il-i) = prim(n,k,j,(il+i-1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        b.x1f(k,j,(il-i)) = -b.x1f(k,j,(il+i  ));  // reflect 1-field
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        b.x2f(k,j,(il-i)) =  b.x2f(k,j,(il+i-1));
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        b.x3f(k,j,(il-i)) =  b.x3f(k,j,(il+i-1));
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX1(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void CellCenteredBoundaryVariable::ReflectOuterX1(
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=nl; i<=nu; ++i) {
          prim(IVX,k,j,iu+i) = -prim(IVX,k,j,(iu-i+1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=nl; i<=nu; ++i) {
          prim(n,k,j,iu+i) = prim(n,k,j,(iu-i+1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        b.x1f(k,j,(iu+i+1)) = -b.x1f(k,j,(iu-i+1));  // reflect 1-field
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        b.x2f(k,j,(iu+i  )) =  b.x2f(k,j,(iu-i+1));
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=nl; i<=nu; ++i) {
        b.x3f(k,j,(iu+i  )) =  b.x3f(k,j,(iu-i+1));
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflecInnerX2(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void CellCenteredBoundaryVariable::ReflectInnerX2(
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=kl; k<=ku; ++k) {
      for (int j=nl; j<=nu; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(IVY,k,jl-j,i) = -prim(IVY,k,jl+j-1,i);  // reflect 2-velocity
        }
      }}
    } else {
      for (int k=kl; k<=ku; ++k) {
      for (int j=nl; j<=nu; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,jl-j,i) = prim(n,k,jl+j-1,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f(k,(jl-j),i) =  b.x1f(k,(jl+j-1),i);
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f(k,(jl-j),i) = -b.x2f(k,(jl+j  ),i);  // reflect 2-field
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f(k,(jl-j),i) =  b.x3f(k,(jl+j-1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX2(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void CellCenteredBoundaryVariable::ReflectOuterX2(
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=kl; k<=ku; ++k) {
      for (int j=nl; j<=nu; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(IVY,k,ju+j,i) = -prim(IVY,k,ju-j+1,i);  // reflect 2-velocity
        }
      }}
    } else {
      for (int k=kl; k<=ku; ++k) {
      for (int j=nl; j<=nu; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,ju+j,i) = prim(n,k,ju-j+1,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f(k,(ju+j  ),i) =  b.x1f(k,(ju-j+1),i);
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f(k,(ju+j+1),i) = -b.x2f(k,(ju-j+1),i);  // reflect 2-field
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=nl; j<=nu; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f(k,(ju+j  ),i) =  b.x3f(k,(ju-j+1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX3(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void CellCenteredBoundaryVariable::ReflectInnerX3(
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=nl; k<=nu; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(IVZ,kl-k,j,i) = -prim(IVZ,kl+k-1,j,i);  // reflect 3-velocity
        }
      }}
    } else {
      for (int k=nl; k<=nu; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,kl-k,j,i) = prim(n,kl+k-1,j,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f((kl-k),j,i) =  b.x1f((kl+k-1),j,i);
      }
    }}

    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f((kl-k),j,i) =  b.x2f((kl+k-1),j,i);
      }
    }}

    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f((kl-k),j,i) = -b.x3f((kl+k  ),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX3(
//                         FaceField &b, const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void CellCenteredBoundaryVariable::ReflectOuterX3(
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int nl, int nu) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=nl; k<=nu; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(IVZ,ku+k,j,i) = -prim(IVZ,ku-k+1,j,i);  // reflect 3-velocity
        }
      }}
    } else {
      for (int k=nl; k<=nu; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,ku+k,j,i) = prim(n,ku-k+1,j,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f((ku+k  ),j,i) =  b.x1f((ku-k+1),j,i);
      }
    }}

    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f((ku+k  ),j,i) =  b.x2f((ku-k+1),j,i);
      }
    }}

    for (int k=nl; k<=nu; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f((ku+k+1),j,i) = -b.x3f((ku-k+1),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}
