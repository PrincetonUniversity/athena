//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow.cpp
//  \brief implementation of outflow BCs in each dimension

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x1 boundary

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        prim(n,k,j,is-i) = prim(n,k,j,is);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                         FaceField &b, Real time, Real dt,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x1 boundary

void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        prim(n,k,j,ie+i) = prim(n,k,j,ie);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x2 boundary

void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(n,k,js-j,i) = prim(n,k,js,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) = b.x1f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = b.x2f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) = b.x3f(k,js,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x2 boundary

void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(n,k,je+j,i) = prim(n,k,je,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) = b.x1f(k,(je  ),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) = b.x3f(k,(je  ),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x3 boundary

void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(n,ks-k,j,i) = prim(n,ks,j,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ks-k),j,i) = b.x1f(ks,j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f((ks-k),j,i) = b.x2f(ks,j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ks-k),j,i) = b.x3f(ks,j,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x3 boundary

void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(n,ke+k,j,i) = prim(n,ke,j,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f((ke+k  ),j,i) = b.x1f((ke  ),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f((ke+k  ),j,i) = b.x2f((ke  ),j,i);
      }
    }}

    for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f((ke+k+1),j,i) = b.x3f((ke+1),j,i);
      }
    }}
  }

  return;
}




// begin fourth-order MHD outflow BC functions-------------------------


void OutflowCellFieldInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &bc,
							 int is, int ie, int js, int je, int ks, int ke) {
  for (int n=0; n<(NFIELD); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=(NGHOST); ++i) {
          bc(n,k,j,is-i) = bc(n,k,j,is);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowFieldInnerX1(MeshBlock *pmb, Coordinates *pco, FaceField &b,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x1 boundary. Applied to arbitrary FaceField
//         variable

void OutflowFaceFieldInnerX1(MeshBlock *pmb, Coordinates *pco, FaceField &b, int is,
                             int ie, int js, int je, int ks, int ke) {
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }
  return;
}


void OutflowCellFieldOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &bc,
							 int is, int ie, int js, int je, int ks, int ke) {
  for (int n=0; n<(NFIELD); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=(NGHOST); ++i) {
          bc(n,k,j,ie+i) = bc(n,k,j,ie);
        }
      }
    }
  }

  return;
}
//----------------------------------------------------------------------------------------
//! \fn void OutflowFieldOuterX1(MeshBlock *pmb, Coordinates *pco, FaceField &b,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x1 boundary. Applied to arbitrary
//         FaceField variable

void OutflowFaceFieldOuterX1(MeshBlock *pmb, Coordinates *pco,FaceField &b, int is,
                             int ie, int js, int je, int ks, int ke) {
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=(NGHOST); ++i) {
          b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
        }
      }
    }
  }

  return;
}

void OutflowCellFieldInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &bc,
							 int is, int ie, int js, int je, int ks, int ke) {
  for (int n=0; n<(NFIELD); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
        for (int i=is; i<=ie; ++i) {
          bc(n,k,js-j,i) = bc(n,k,js,i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowFieldInnerX2(MeshBlock *pmb, Coordinates *pco, FaceField &b,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x2 boundary. Applied to arbitrary FaceField
//         variable

void OutflowFaceFieldInnerX2(MeshBlock *pmb, Coordinates *pco, FaceField &b, int is,
                             int ie, int js, int je, int ks, int ke) {
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(js-j),i) = b.x1f(k,js,i);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(js-j),i) = b.x2f(k,js,i);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(js-j),i) = b.x3f(k,js,i);
        }
      }
    }
  }
  return;
}


void OutflowCellFieldOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &bc,
							 int is, int ie, int js, int je, int ks, int ke) {
  for (int n=0; n<(NFIELD); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
        for (int i=is; i<=ie; ++i) {
          bc(n,k,je+j,i) = bc(n,k,je,i);
        }
      }
    }
  }

  return;
}
//----------------------------------------------------------------------------------------
//! \fn void OutflowFieldOuterX2(MeshBlock *pmb, Coordinates *pco, FaceField &b,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x2 boundary. Applied to arbitrary
//         FaceField variable

void OutflowFaceFieldOuterX2(MeshBlock *pmb, Coordinates *pco,FaceField &b, int is,
                             int ie, int js, int je, int ks, int ke) {
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x1f(k,(je+1+j),i) = b.x1f(k,(je+1),i);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(je+1+j),i) = b.x2f(k,(je+1),i);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(je+1+j),i) = b.x3f(k,(je+1),i);
        }
      }
    }
  }

  return;
}


// -------- end fourth-order MHD outflow BC functions
