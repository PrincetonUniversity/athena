//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file coordinates.cpp
//  \brief implements functions for Coordinates abstract base class

// C/C++ headers
#include <algorithm>

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../bvals/bvals.hpp"

//----------------------------------------------------------------------------------------
// Coordinates constructor: sets coordinates and coordinate spacing of cell FACES

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin, bool flag) {
  pmy_block = pmb;
  coarse_flag=flag;
  int is, ie, js, je, ks, ke, ng;
  if (coarse_flag==true) {
    is = pmb->cis; js = pmb->cjs; ks = pmb->cks;
    ie = pmb->cie; je = pmb->cje; ke = pmb->cke;
    ng=pmb->cnghost;
  } else {
    is = pmb->is; js = pmb->js; ks = pmb->ks;
    ie = pmb->ie; je = pmb->je; ke = pmb->ke;
    ng=NGHOST;
  }
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // allocate arrays for face-centered coordinates and coordinate spacing
  int ncells1 = (ie-is+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (je-js+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ke-ks+1) + 2*ng;

  // note extra cell for face-positions
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));

  int64_t nrootmesh, noffset;
  int64_t &lx1=pmy_block->loc.lx1;
  int64_t &lx2=pmy_block->loc.lx2;
  int64_t &lx3=pmy_block->loc.lx3;
  int &ll=pmy_block->loc.level;

//--- X1-DIRECTION: initialize coordinates and spacing of cell FACES (x1f,dx1f)

  nrootmesh=mesh_size.nx1*(1L<<(ll-pm->root_level));

  // use nonuniform or user-defined meshgen fn
  if (pm->use_uniform_meshgen_fn_[X1DIR]==false) {
    for (int i=is-ng; i<=ie+ng+1; ++i) {
      // if there are too many levels, this won't work or be precise enough
      if (coarse_flag == false) {
        noffset = static_cast<int64_t>(i-is + lx1*block_size.nx1);
      } else {
        noffset = static_cast<int64_t>((i-is)*2 + lx1*block_size.nx1);
      }
      Real rx = ComputeMeshGeneratorX(noffset, nrootmesh, false);
      x1f(i) = pm->MeshGenerator_[X1DIR](rx, mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
    for (int i=is-ng; i<=ie+ng; ++i) {
      dx1f(i) = x1f(i+1) - x1f(i);
    }

    // check that coordinate spacing is reasonable
    Real rmax=1.0, rmin=1.0;
    for (int i=is; i<=ie; i++) {
      rmax=std::max(dx1f(i+1)/dx1f(i),rmax);
      rmin=std::min(dx1f(i+1)/dx1f(i),rmin);
    }
    if (rmax > 1.1 || rmin  < 1.0/1.1) {
       std::cout << "### Warning in Coordinates constructor" << std::endl
         << "Neighboring cell sizes differ by more than 10% in the x1 direction."
         << std::endl;
    }

  } else {
    // uniform grid: use UniformMeshGeneratorX1()
    Real dx=(block_size.x1max-block_size.x1min)/(ie-is+1);
    for (int i=is-ng; i<=ie+ng+1; ++i) {
      // if there are too many levels, this won't work or be precise enough
      if (coarse_flag == false) {
        noffset = static_cast<int64_t>(i-is + lx1*block_size.nx1);
      } else {
        noffset = static_cast<int64_t>((i-is)*2 + lx1*block_size.nx1);
      }
      Real rx = ComputeMeshGeneratorX(noffset, nrootmesh, true);
      x1f(i) = pm->MeshGenerator_[X1DIR](rx, mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;

    for (int i=is-ng; i<=ie+ng; ++i) {
      dx1f(i)=dx;
    }
  }

  // correct cell face coordinates in ghost zones for reflecting boundary condition
  if (pmy_block->pbval->block_bcs[INNER_X1] == REFLECTING_BNDRY) {
    for (int i=1; i<=ng; ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }
  if (pmy_block->pbval->block_bcs[OUTER_X1] == REFLECTING_BNDRY) {
    for (int i=1; i<=ng; ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

//--- X2-DIRECTION: initialize coordinates and spacing of cell FACES (x2f,dx2f)

  if (ncells2 > 1) {
    nrootmesh=mesh_size.nx2*(1L<<(ll-pm->root_level));

    // use nonuniform or user-defined meshgen fn
    if (pm->use_uniform_meshgen_fn_[X2DIR]==false) {
      for (int j=js-ng; j<=je+ng+1; ++j) {
        // if there are too many levels, this won't work or be precise enough
        if (coarse_flag == false) {
          noffset = static_cast<int64_t>(j-js + lx2*block_size.nx2);
        } else {
          noffset = static_cast<int64_t>((j-js)*2 + lx2*block_size.nx2);
        }
        Real rx = ComputeMeshGeneratorX(noffset, nrootmesh, false);
        x2f(j) = pm->MeshGenerator_[X2DIR](rx, mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
      for (int j=js-ng; j<=je+ng; ++j) {
        dx2f(j)=x2f(j+1)-x2f(j);
      }

      // check that coordinate spacing is reasonable
      Real rmax=1.0, rmin=1.0;
      for (int j=pmy_block->js; j<=pmy_block->je; j++) {
        rmax=std::max(dx2f(j+1)/dx2f(j),rmax);
        rmin=std::min(dx2f(j+1)/dx2f(j),rmin);
      }
      if (rmax > 1.1 || rmin  < 1.0/1.1) {
         std::cout << "### Warning in Coordinates constructor" << std::endl
           << "Neighboring cell sizes differ by more than 10% in the x2 direction."
           << std::endl;
      }

    } else {
      // uniform grid: use UniformMeshGeneratorX2()
      Real dx=(block_size.x2max-block_size.x2min)/(je-js+1);
      for (int j=js-ng; j<=je+ng+1; ++j) {
        if (coarse_flag == false) {
          noffset = static_cast<int64_t>(j-js + lx2*block_size.nx2);
        } else {
          noffset = static_cast<int64_t>((j-js)*2 + lx2*block_size.nx2);
        }
        Real rx = ComputeMeshGeneratorX(noffset, nrootmesh, true);
        x2f(j) = pm->MeshGenerator_[X2DIR](rx, mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;

      for (int j=js-ng; j<=je+ng; ++j) {
        dx2f(j) = dx;
      }
    }

    // correct cell face coordinates in ghost zones for reflecting boundary condition
    if (pmy_block->pbval->block_bcs[INNER_X2] == REFLECTING_BNDRY
     || pmy_block->pbval->block_bcs[INNER_X2] == POLAR_BNDRY) { // also polar boundary
      for (int j=1; j<=ng; ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }
    if (pmy_block->pbval->block_bcs[OUTER_X2] == REFLECTING_BNDRY
     || pmy_block->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY) { // also polar boundary
      for (int j=1; j<=ng; ++j) {
        dx2f(je+j  ) = dx2f(je-j+1);
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    }

  // 1D problem
  } else {
    dx2f(js  ) = block_size.x2max-block_size.x2min;
    x2f (js  ) = block_size.x2min;
    x2f (je+1) = block_size.x2max;
  }

//--- X3-DIRECTION: initialize coordinates and spacing of cell FACES (x3f,dx3f)

  if (ncells3 > 1) {
    nrootmesh=mesh_size.nx3*(1L<<(ll-pm->root_level));

    // use nonuniform or user-defined meshgen fn
    if (pm->use_uniform_meshgen_fn_[X3DIR]==false) {
      for (int k=ks-ng; k<=ke+ng+1; ++k) {
        // if there are too many levels, this won't work or be precise enough
        if (coarse_flag == false) {
          noffset = static_cast<int64_t>(k-ks + lx3*block_size.nx3);
        } else {
          noffset = static_cast<int64_t>((k-ks)*2 + lx3*block_size.nx3);
        }
        Real rx = ComputeMeshGeneratorX(noffset, nrootmesh, false);
        x3f(k) = pm->MeshGenerator_[X3DIR](rx, mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
      for (int k=ks-ng; k<=ke+ng; ++k) {
        dx3f(k)=x3f(k+1)-x3f(k);
      }

      // check that coordinate spacing is reasonable
      Real rmax=1.0, rmin=1.0;
      for (int k=pmy_block->ks; k<=pmy_block->ke; k++) {
        rmax=std::max(dx3f(k+1)/dx3f(k),rmax);
        rmin=std::min(dx3f(k+1)/dx3f(k),rmin);
      }
      if (rmax > 1.1 || rmin  < 1.0/1.1) {
         std::cout << "### Warning in Coordinates constructor" << std::endl
           << "Neighboring cell sizes differ by more than 10% in the x3 direction."
           << std::endl;
      }

    } else {
      // uniform grid: use UniformMeshGeneratorX3()
      Real dx=(block_size.x3max-block_size.x3min)/(ke-ks+1);
      for (int k=ks-ng; k<=ke+ng+1; ++k) {
        if (coarse_flag == false) {
          noffset = static_cast<int64_t>(k-ks + lx3*block_size.nx3);
        } else {
          noffset = static_cast<int64_t>((k-ks)*2 + lx3*block_size.nx3);
        }
        Real rx = ComputeMeshGeneratorX(noffset, nrootmesh, true);
        x3f(k) = pm->MeshGenerator_[X3DIR](rx, mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;

      for (int k=ks-ng; k<=ke+ng; ++k) {
        dx3f(k) = dx;
      }
    }

    // correct cell face coordinates in ghost zones for reflecting boundary condition
    if (pmy_block->pbval->block_bcs[INNER_X3] == REFLECTING_BNDRY) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }
    if (pmy_block->pbval->block_bcs[OUTER_X3] == REFLECTING_BNDRY) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ke+k  ) = dx3f(ke-k+1);
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    }

  // 1D or 2D problem
  } else {
    dx3f(ks) = block_size.x3max-block_size.x3min;
    x3f(ks  ) = block_size.x3min;
    x3f(ke+1) = block_size.x3max;
  }

}

// destructor

Coordinates::~Coordinates() {
  dx1f.DeleteAthenaArray();
  dx2f.DeleteAthenaArray();
  dx3f.DeleteAthenaArray();
  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))

void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    len(i) = dx1f(i);
  }
  return;
}

// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    len(i) = dx2f(j);
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    len(i) = dx3f(k);
  }
  return;
}


//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X at (i,j,k)

Real Coordinates::GetEdge1Length(const int k, const int j, const int i) {
  return dx1f(i);
}

Real Coordinates::GetEdge2Length(const int k, const int j, const int i) {
  return dx2f(j);
}

Real Coordinates::GetEdge3Length(const int k, const int j, const int i) {
  return dx3f(k);
}

//----------------------------------------------------------------------------------------
// VolCenterXLength functions: compute physical length connecting cell centers as vector
// VolCenter1(i,j,k) located at (i+1/2,j,k), i.e. (x1f(i+1), x2v(j), x3v(k))
void Coordinates::VolCenter1Length(const int k, const int j, const int il, const int iu,
                                   AthenaArray<Real> &len) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
        len(i) = dx1v(i);
    }
    return;
}
void Coordinates::VolCenter2Length(const int k, const int j, const int il, const int iu,
                                   AthenaArray<Real> &len) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
        len(i) = dx2v(j);
    }
    return;
}
void Coordinates::VolCenter3Length(const int k, const int j, const int il, const int iu,
                                   AthenaArray<Real> &len) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
        len(i) = dx3v(k);
    }
    return;
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void Coordinates::CenterWidth1(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx1) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx1(i) = dx1f(i);
  }
  return;
}

void Coordinates::CenterWidth2(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx2) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx2(i) = dx2f(j);
  }
  return;
}

void Coordinates::CenterWidth3(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx3) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx3(i) = dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area) {
#pragma nounroll
  for (int i=il; i<=iu; ++i) {
    // area1 = dy dz
    Real& area_i = area(i);
    area_i = dx2f(j)*dx3f(k);
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area) {
#pragma nounroll
  for (int i=il; i<=iu; ++i) {
    // area2 = dx dz
    Real& area_i = area(i);
    area_i = dx1f(i)*dx3f(k);
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area) {
#pragma nounroll
  for (int i=il; i<=iu; ++i) {
    // area3 = dx dy
    Real& area_i = area(i);
    area_i = dx1f(i)*dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real Coordinates::GetFace1Area(const int k, const int j, const int i) {
  return dx2f(j)*dx3f(k);
}

Real Coordinates::GetFace2Area(const int k, const int j, const int i) {
  return dx1f(i)*dx3f(k);
}

Real Coordinates::GetFace3Area(const int k, const int j, const int i) {
  return dx1f(i)*dx2f(j);
}

//----------------------------------------------------------------------------------------
// VolCenterFaceXArea functions: compute area of face with normal in X-dir as vector
// where the faces are joined by cell centers (for non-ideal MHD)

void Coordinates::VolCenterFace1Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx2v(j)*dx3v(k);
  }
  return;
}

void Coordinates::VolCenterFace2Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx1v(i)*dx3v(k);
  }
  return;
}

void Coordinates::VolCenterFace3Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    Real& area_i = area(i);
    area_i = dx1v(i)*dx2v(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // volume = dx dy dz
    Real& vol_i = vol(i);
    vol_i = dx1f(i)*dx2f(j)*dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real Coordinates::GetCellVolume(const int k, const int j, const int i) {
  return dx1f(i)*dx2f(j)*dx3f(k);
}

//-------------------------------------------------------------------------------------
// Laplacian: calculate total Laplacian of 4D scalar array s() to second order accuracy
// may need to replace dx*f with dx*v for nonuniform coordinates for some applications

void Coordinates::Laplacian(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                            const int il, const int iu, const int jl, const int ju,
                            const int kl, const int ku, const int nl, const int nu) {
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          delta_s(n,k,j,i) = (s(n,k,j,i-1) - 2.0*s(n,k,j,i) + s(n,k,j,i+1))
              /(dx1f(i)*dx1f(i));
        }
        if (pmy_block->block_size.nx2 > 1) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) += (s(n,k,j-1,i) - 2.0*s(n,k,j,i) + s(n,k,j+1,i))
                /(dx2f(j)*dx2f(j));
          }
        }
        if (pmy_block->block_size.nx3 > 1) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) += (s(n,k-1,j,i) - 2.0*s(n,k,j,i) + s(n,k+1,j,i))
                /(dx3f(k)*dx3f(k));
          }
        }
      }
    }
  }
  return;
}

//-------------------------------------------------------------------------------------
// LaplacianX* functions: calculate Laplacian in subspaces orthogonal to X-dir

void Coordinates::LaplacianX1(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                              const int il, const int iu, const int jl, const int ju,
                              const int kl, const int ku, const int nl, const int nu) {
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        if (pmy_block->block_size.nx2 > 1) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) = (s(n,k,j-1,i) - 2.0*s(n,k,j,i) + s(n,k,j+1,i))
                /(dx2f(j)*dx2f(j));
          }
        } else { // 1D domain
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) = 0.0;
          }
        }
        if (pmy_block->block_size.nx3 > 1) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) += (s(n,k-1,j,i) - 2.0*s(n,k,j,i) + s(n,k+1,j,i))
                /(dx3f(k)*dx3f(k));
          }
        }
      }
    }
  }
  return;
}

void Coordinates::LaplacianX2(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                              const int il, const int iu, const int jl, const int ju,
                              const int kl, const int ku, const int nl, const int nu) {
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          delta_s(n,k,j,i) = (s(n,k,j,i-1) - 2.0*s(n,k,j,i) + s(n,k,j,i+1))
              /(dx1f(i)*dx1f(i));
        }
        if (pmy_block->block_size.nx3 > 1) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) += (s(n,k-1,j,i) - 2.0*s(n,k,j,i) + s(n,k+1,j,i))
                /(dx3f(k)*dx3f(k));
          }
        }
      }
    }
  }
  return;
}

void Coordinates::LaplacianX3(const AthenaArray<Real> &s, AthenaArray<Real> &delta_s,
                              const int il, const int iu, const int jl, const int ju,
                              const int kl, const int ku, const int nl, const int nu) {
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          delta_s(n,k,j,i) = (s(n,k,j,i-1) - 2.0*s(n,k,j,i) + s(n,k,j,i+1))
              /(dx1f(i)*dx1f(i));
        }
        if (pmy_block->block_size.nx2 > 1) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            delta_s(n,k,j,i) += (s(n,k,j-1,i) - 2.0*s(n,k,j,i) + s(n,k,j+1,i))
                /(dx2f(j)*dx2f(j));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function
void Coordinates::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u) {
  return;
}

//----------------------------------------------------------------------------------------
// Function for determining if index corresponds to a polar boundary
// Inputs:
//   j: x2-index
// Outputs:
//   returned value: true if face indexed with j is on a pole; false otherwise

bool Coordinates::IsPole(int j) {
  if ((pmy_block->pbval->block_bcs[INNER_X2] == POLAR_BNDRY ||
       pmy_block->pbval->block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) &&
      j == pmy_block->js) {
    return true;
  }
  if ((pmy_block->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY ||
       pmy_block->pbval->block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) &&
      j == pmy_block->je+1) {
    return true;
  }
  return false;
}

//----------------------------------------------------------------------------------------
// Function for implementing user-defined metric
// Inputs:
//   x1,x2,x3: spatial location of point
//   pin: pointer to runtime inputs
// Outputs:
//   g,g_inv: arrays of metric covariant and contravariant components
//   dg_dx1,dg_dx2,dg_dx3: arrays of spatial derivatives of covariant components

void Coordinates::Metric(Real x1, Real x2, Real x3, ParameterInput *pin,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv, AthenaArray<Real> &dg_dx1,
    AthenaArray<Real> &dg_dx2, AthenaArray<Real> &dg_dx3) {
  pmy_block->pmy_mesh->UserMetric_(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);
  return;
}
