//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "bvals.hpp"

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(IVX,k,j,is-i) = -a(IVX,k,j,(is+i-1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(n,k,j,is-i) = a(n,k,j,(is+i-1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) { 
        b.x1f(k,j,(is-i)) = -b.x1f(k,j,(is+i  ));  // reflect 1-field
      } 
    }}
  
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,(is+i-1));
      }
    }}  
        
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) =  b.x3f(k,j,(is+i-1));
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(IVX,k,j,ie+i) = -a(IVX,k,j,(ie-i+1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          a(n,k,j,ie+i) = a(n,k,j,(ie-i+1));
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(ie+i+1)) = -b.x1f(k,j,(ie-i+1));  // reflect 1-field
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i  )) =  b.x2f(k,j,(ie-i+1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i  )) =  b.x3f(k,j,(ie-i+1));
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflecInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVY,k,js-j,i) = -a(IVY,k,js+j-1,i);  // reflect 2-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,k,js-j,i) = a(n,k,js+j-1,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) =  b.x1f(k,(js+j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = -b.x2f(k,(js+j  ),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) =  b.x3f(k,(js+j-1),i);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVY)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVY,k,je+j,i) = -a(IVY,k,je-j+1,i);  // reflect 2-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,k,je+j,i) = a(n,k,je-j+1,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) =  b.x1f(k,(je-j+1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = -b.x2f(k,(je-j+1),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) =  b.x3f(k,(je-j+1),i);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVZ,ks-k,j,i) = -a(IVZ,ks+k-1,j,i);  // reflect 3-velocity
        }
      }}
    } else {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,ks-k,j,i) = a(n,ks+k-1,j,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ks-k),j,i) =  b.x1f((ks+k-1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f((ks-k),j,i) =  b.x2f((ks+k-1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ks-k),j,i) = -b.x3f((ks+k  ),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(IVZ,ke+k,j,i) = -a(IVZ,ke-k+1,j,i);  // reflect 3-velocity
        }
      }}
    } else {
      for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          a(n,ke+k,j,i) = a(n,ke-k+1,j,i);
        }
      }}
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ke+k  ),j,i) =  b.x1f((ke-k+1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
        b.x2f((ke+k  ),j,i) =  b.x2f((ke-k+1),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ke+k+1),j,i) = -b.x3f((ke-k+1),j,i);  // reflect 3-field
      }
    }}
  }

  return;
}
