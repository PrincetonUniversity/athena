#ifndef COORDINATES_HPP
#define COORDINATES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file coordinates.hpp
//  \brief defines Coordinates class used to compute/store geometrical factors (areas,
//  volumes, source terms) related to a Mesh
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"

// forward declarations
class MeshBlock;
class ParameterInput;

//! \class Coordinates
//  \brief coordinate data and functions

class Coordinates {
public:
  friend class HydroSourceTerms;
  Coordinates(MeshBlock *pmb, ParameterInput *pin);
  ~Coordinates();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Coordinates

  // coordinate arrays
  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions

  // arrays for mesh refinement
  AthenaArray<Real> coarse_x1f,  coarse_x2f,  coarse_x3f;
  AthenaArray<Real> coarse_x1v,  coarse_x2v,  coarse_x3v;

  // arrays for MHD mesh refinement
  //note: only x1s2 and x1s3 in spherical and schwarzschild coordinates are non trivial
  AthenaArray<Real> x1s2, x1s3, x2s1, x2s3, x3s1, x3s2; // area averaged positions
  AthenaArray<Real> coarse_x1s2, coarse_x1s3, coarse_x2s1,
                    coarse_x2s3, coarse_x3s1, coarse_x3s2;

  void AllocateAndSetBasicCoordinates(void);
  void DeleteBasicCoordinates(void);

  // Functions for returning private variables
  Real GetMass() const {return bh_mass_;}
  Real GetSpin() const {return bh_spin_;}

// functions to compute length of edges
  void Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);

  Real GetEdge2Length(const int k, const int j, const int i);
  Real GetEdge3Length(const int k, const int j, const int i);

// functions to compute physical width at cell center
  Real CenterWidth1(const int k, const int j, const int i);
  Real CenterWidth2(const int k, const int j, const int i);
  Real CenterWidth3(const int k, const int j, const int i);

// functions to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);

  Real GetFace1Area(const int k, const int j, const int i);

// function to compute volume of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

// function to compute geometrical source terms
  void CoordSrcTerms(const int k, const int j, const Real dt,
    const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

// Functions to calculate covariant derivatives at faces, for viscosity calculations  
  void FaceXdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceXdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceXdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);
  void FaceZdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len);

// function to compute Divv
  void Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv);

// function to compute viscous source terms
  void VisSrcTermsX1(const int k, const int j, const Real dt,
    const AthenaArray<Real> &flx,
    const AthenaArray<Real> &prim, AthenaArray<Real> &u);

  void VisSrcTermsX2(const int k, const int j, const Real dt,
    const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
    const AthenaArray<Real> &prim, AthenaArray<Real> &u);

  void VisSrcTermsX3(const int k, const int j, const Real dt,
    const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
    const AthenaArray<Real> &prim, AthenaArray<Real> &u);

  // Functions for use in general relativity
  #if GENERAL_RELATIVITY  // declare, but do not define, in GR case
    void CellMetric(const int k, const int j, const int il, const int iu,
        AthenaArray<Real> &g, AthenaArray<Real> &gi);
    void Face1Metric(const int k, const int j, const int il, const int iu,
        AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
    void Face2Metric(const int k, const int j, const int il, const int iu,
        AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
    void Face3Metric(const int k, const int j, const int il, const int iu,
        AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
    void PrimToLocal1(const int k, const int j, const int il, const int iu,
        const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
        AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
    void PrimToLocal2(const int k, const int j, const int il, const int iu,
        const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
        AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
    void PrimToLocal3(const int k, const int j, const int il, const int iu,
        const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
        AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
    void FluxToGlobal1(const int k, const int j, const int il, const int iu,
        const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
        AthenaArray<Real> &flux);
    void FluxToGlobal2(const int k, const int j, const int il, const int iu,
        const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
        AthenaArray<Real> &flux);
    void FluxToGlobal3(const int k, const int j, const int il, const int iu,
        const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
        AthenaArray<Real> &flux);
    Real DistanceBetweenPoints(Real a1, Real a2, Real a3, Real bx, Real by, Real bz);
    void MinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3,
        Real *pt, Real *px, Real *py, Real *pz);
    void TransformVectorCell(Real at, Real ax, Real ay, Real az, int k, int j, int i,
        Real *a0, Real *a1, Real *a2, Real *a3);
    void TransformVectorFace1(Real at, Real ax, Real ay, Real az, int k, int j, int i,
        Real *a0, Real *a1, Real *a2, Real *a3);
    void TransformVectorFace2(Real at, Real ax, Real ay, Real az, int k, int j, int i,
        Real *a0, Real *a1, Real *a2, Real *a3);
    void TransformVectorFace3(Real at, Real ax, Real ay, Real az, int k, int j, int i,
        Real *a0, Real *a1, Real *a2, Real *a3);
    void LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
        Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3);
    void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3,
        Real *pr, Real *ptheta, Real *pphi);
  #else  // define no-op functions otherwise not defined in non-GR case
    void CellMetric(const int, const int, const int, const int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void Face1Metric(const int, const int, const int, const int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void Face2Metric(const int, const int, const int, const int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void Face3Metric(const int, const int, const int, const int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void PrimToLocal1(const int, const int, const int, const int,
        const AthenaArray<Real> &, AthenaArray<Real> &, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void PrimToLocal2(const int, const int, const int, const int,
        const AthenaArray<Real> &, AthenaArray<Real> &, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void PrimToLocal3(const int, const int, const int, const int,
        const AthenaArray<Real> &, AthenaArray<Real> &, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FluxToGlobal1(const int, const int, const int, const int,
        const AthenaArray<Real> &, const AthenaArray<Real> &, AthenaArray<Real> &)
        {return;}
    void FluxToGlobal2(const int, const int, const int, const int,
        const AthenaArray<Real> &, const AthenaArray<Real> &, AthenaArray<Real> &)
        {return;}
    void FluxToGlobal3(const int, const int, const int, const int,
        const AthenaArray<Real> &, const AthenaArray<Real> &, AthenaArray<Real> &)
        {return;}
    Real DistanceBetweenPoints(Real, Real, Real, Real, Real, Real) {return 0.0;}
    void MinkowskiCoordinates(Real, Real, Real, Real, Real *, Real *, Real *, Real *)
        {return;}
    void TransformVectorCell(Real, Real, Real, Real, int, int, int, Real *, Real *,
        Real *, Real *) {return;}
    void TransformVectorFace1(Real, Real, Real, Real, int, int, int, Real *, Real *,
        Real *, Real *) {return;}
    void TransformVectorFace2(Real, Real, Real, Real, int, int, int, Real *, Real *,
        Real *, Real *) {return;}
    void TransformVectorFace3(Real, Real, Real, Real, int, int, int, Real *, Real *,
        Real *, Real *) {return;}
    void LowerVectorCell(Real, Real, Real, Real, int, int, int, Real *, Real *, Real *,
        Real *) {return;}
  #endif  // GENERAL_RELATIVITY

private:

  // Constant parameters for various coordinate systems
  Real bh_mass_;          // M: Schwarzschild
  Real bh_spin_;          // a (dimensionless): Schwarzschild
  Real sinu_amplitude_;   // a: sinusoidal
  Real sinu_wavenumber_;  // k: sinusoidal
  Real tilted_a_;         // a: tilted

  // Scratch arrays for coordinate factors
  // Format: coord_<type>[<direction>]_<index>[<count>]_
  //   type: vol[ume], area, etc.
  //   direction: 1/2/3 depending on which face, edge, etc. is in play
  //   index: i/j/k indicating which coordinates index array
  //   count: 1/2/... in case multiple arrays are needed for different terms
  AthenaArray<Real> coord_vol_i_, coord_vol_i1_, coord_vol_i2_;
  AthenaArray<Real> coord_vol_j_, coord_vol_j1_, coord_vol_j2_;
  AthenaArray<Real> coord_vol_k1_;
  AthenaArray<Real> coord_area1_i_, coord_area1_i1_;
  AthenaArray<Real> coord_area1_j_, coord_area1_j1_, coord_area1_j2_;
  AthenaArray<Real> coord_area1_k1_;
  AthenaArray<Real> coord_area2_i_, coord_area2_i1_, coord_area2_i2_;
  AthenaArray<Real> coord_area2_j_, coord_area2_j1_, coord_area2_j2_;
  AthenaArray<Real> coord_area2_k1_;
  AthenaArray<Real> coord_area3_i_, coord_area3_i1_, coord_area3_i2_;
  AthenaArray<Real> coord_area3_j1_, coord_area3_j2_;
  AthenaArray<Real> coord_len1_i1_, coord_len1_i2_;
  AthenaArray<Real> coord_len1_j1_, coord_len1_j2_;
  AthenaArray<Real> coord_len2_i1_;
  AthenaArray<Real> coord_len2_j1_, coord_len2_j2_;
  AthenaArray<Real> coord_len3_i1_;
  AthenaArray<Real> coord_len3_j1_, coord_len3_j2_;
  AthenaArray<Real> coord_len3_k1_;
  AthenaArray<Real> coord_width1_i1_;
  AthenaArray<Real> coord_width2_i1_;
  AthenaArray<Real> coord_width2_j1_;
  AthenaArray<Real> coord_width3_j1_, coord_width3_j2_, coord_width3_j3_;
  AthenaArray<Real> coord_width3_k1_;
  AthenaArray<Real> coord_width3_ji1_;
  AthenaArray<Real> coord_src_i1_, coord_src_i2_, coord_src_i3_, coord_src_i4_;
  AthenaArray<Real> coord_src_j1_, coord_src_j2_, coord_src_j3_;
  AthenaArray<Real> coord_src1_i_;
  AthenaArray<Real> coord_src1_j_;
  AthenaArray<Real> coord_src2_i_;
  AthenaArray<Real> coord_src2_j_;

  // Scratch arrays for physical source terms
  AthenaArray<Real> phy_src1_i_, phy_src2_i_;

  // GR-specific scratch arrays
  AthenaArray<Real> metric_cell_i1_, metric_cell_i2_;
  AthenaArray<Real> metric_cell_j1_, metric_cell_j2_;
  AthenaArray<Real> metric_face1_i1_, metric_face1_i2_;
  AthenaArray<Real> metric_face1_j1_, metric_face1_j2_;
  AthenaArray<Real> metric_face2_i1_, metric_face2_i2_;
  AthenaArray<Real> metric_face2_j1_, metric_face2_j2_;
  AthenaArray<Real> metric_face3_i1_, metric_face3_i2_;
  AthenaArray<Real> metric_face3_j1_, metric_face3_j2_;
  AthenaArray<Real> trans_face1_i1_, trans_face1_i2_;
  AthenaArray<Real> trans_face1_j1_;
  AthenaArray<Real> trans_face1_ji1_, trans_face1_ji2_, trans_face1_ji3_,
      trans_face1_ji4_, trans_face1_ji5_, trans_face1_ji6_, trans_face1_ji7_;
  AthenaArray<Real> trans_face2_i1_, trans_face2_i2_;
  AthenaArray<Real> trans_face2_j1_;
  AthenaArray<Real> trans_face2_ji1_, trans_face2_ji2_, trans_face2_ji3_,
      trans_face2_ji4_, trans_face2_ji5_, trans_face2_ji6_;
  AthenaArray<Real> trans_face3_i1_, trans_face3_i2_;
  AthenaArray<Real> trans_face3_j1_;
  AthenaArray<Real> trans_face3_ji1_, trans_face3_ji2_, trans_face3_ji3_,
      trans_face3_ji4_, trans_face3_ji5_, trans_face3_ji6_;
  AthenaArray<Real> g_, gi_;
};


inline void Coordinates::AllocateAndSetBasicCoordinates(void)
{
  RegionSize& block_size = pmy_block->block_size;
  Mesh *pm=pmy_block->pmy_mesh;

  // allocate arrays for sizes and positions of cells
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);

  // cell sizes
  dx1f.NewAthenaArray(ncells1);
  dx2f.NewAthenaArray(ncells2);
  dx3f.NewAthenaArray(ncells3);
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);

  // cell positions. Note the extra element for cell face positions
  x1f.NewAthenaArray((ncells1+1));
  x2f.NewAthenaArray((ncells2+1));
  x3f.NewAthenaArray((ncells3+1));
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  if(pm->multilevel==true) {
    int cnghost=pmy_block->cnghost;
    int ncc1=block_size.nx1/2+2*cnghost;
    int ncc2=1, ncc3=1;
    if(block_size.nx2>1) // 2D or 3D
      ncc2=block_size.nx2/2+2*cnghost;
    if(block_size.nx3>1) // 3D
      ncc3=block_size.nx3/2+2*cnghost;

    // cell positions. Note the extra element for cell face positions
    coarse_x1f.NewAthenaArray(ncc1+1);
    coarse_x2f.NewAthenaArray(ncc2+1);
    coarse_x3f.NewAthenaArray(ncc3+1);
    coarse_x1v.NewAthenaArray(ncc1);
    coarse_x2v.NewAthenaArray(ncc2);
    coarse_x3v.NewAthenaArray(ncc3);

    if (MAGNETIC_FIELDS_ENABLED) {
      // area weighted positions for AMR/SMR MHD
      x1s2.NewAthenaArray(ncells1);
      x1s3.NewAthenaArray(ncells1);
      x2s1.NewAthenaArray(ncells2);
      x2s3.NewAthenaArray(ncells2);
      x3s1.NewAthenaArray(ncells3);
      x3s2.NewAthenaArray(ncells3);
      coarse_x1s2.NewAthenaArray(ncc1);
      coarse_x1s3.NewAthenaArray(ncc1);
      coarse_x2s1.NewAthenaArray(ncc2);
      coarse_x2s3.NewAthenaArray(ncc2);
      coarse_x3s1.NewAthenaArray(ncc3);
      coarse_x3s2.NewAthenaArray(ncc3);
    }
  }

  int is=pmy_block->is, js=pmy_block->js, ks=pmy_block->ks;
  int ie=pmy_block->ie, je=pmy_block->je, ke=pmy_block->ke;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  long long nrootmesh, noffset;
  int root_level;
  long int &lx1=pmy_block->loc.lx1;
  long int &lx2=pmy_block->loc.lx2;
  long int &lx3=pmy_block->loc.lx3;
  int &ll=pmy_block->loc.level;

// X1-DIRECTION: initialize sizes and positions of cell FACES (dx1f,x1f)
  nrootmesh=mesh_size.nx1*(1L<<(ll-pm->root_level));
  if(block_size.x1rat == 1.0) { // uniform
    Real dx=(block_size.x1max-block_size.x1min)/block_size.nx1;
    for(int i=is-NGHOST; i<=ie+NGHOST; ++i)
      dx1f(i)=dx;
    x1f(is-NGHOST)=block_size.x1min-NGHOST*dx;
    for(int i=is-NGHOST+1;i<=ie+NGHOST+1;i++)
      x1f(i)=x1f(i-1)+dx;
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
  }
  else {
    for (int i=is-NGHOST; i<=ie+NGHOST+1; ++i) {
       // if there are too many levels, this won't work or be precise enough
      noffset=(i-is)+(long long)lx1*block_size.nx1;
      Real rx=(Real)noffset/(Real)nrootmesh;
      x1f(i)=pm->MeshGeneratorX1(rx,mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
    for(int i=is-NGHOST; i<=ie+NGHOST; ++i)
      dx1f(i)=x1f(i+1)-x1f(i);
  }

// correct cell face positions in ghost zones for reflecting boundary condition
  if (pmy_block->block_bcs[inner_x1] == 1) {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }
  if (pmy_block->block_bcs[outer_x1] == 1) {
    for (int i=1; i<=(NGHOST); ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

// X2-DIRECTION: initialize spacing and positions of cell FACES (dx2f,x2f)
  if(block_size.nx2 > 1) {
    nrootmesh=mesh_size.nx2*(1L<<(ll-root_level));
    if(block_size.x2rat == 1.0) { // uniform
      Real dx=(block_size.x2max-block_size.x2min)/block_size.nx2;
      for(int j=js-NGHOST; j<=je+NGHOST; ++j)
        dx2f(j)=dx;
      x2f(js-NGHOST)=block_size.x2min-NGHOST*dx;
      for(int j=js-NGHOST+1;j<=je+NGHOST+1;j++)
        x2f(j)=x2f(j-1)+dx;
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
    }
    else {
      for (int j=js-NGHOST; j<=je+NGHOST+1; ++j) {
         // if there are too many levels, this won't work or be precise enough
        noffset=(j-js)+(long long)lx2*block_size.nx2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x2f(j)=pm->MeshGeneratorX2(rx,mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
      for(int j=js-NGHOST; j<=je+NGHOST; ++j)
        dx2f(j)=x2f(j+1)-x2f(j);
    }

  // correct cell face positions in ghost zones for reflecting boundary condition
    if (pmy_block->block_bcs[inner_x2] == 1) {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }
    if (pmy_block->block_bcs[outer_x2] == 1) {
      for (int j=1; j<=(NGHOST); ++j) {
        dx2f(je+j  ) = dx2f(je-j+1);
         x2f(je+j+1) =  x2f(je+j) + dx2f(je+j);
      }
    }
  }
  else {
    dx2f(js) = block_size.x2max-block_size.x2min;
    x2f(js  ) = block_size.x2min;
    x2f(je+1) = block_size.x2max;
  }


// X3-DIRECTION: initialize spacing and positions of cell FACES (dx3f,x3f)
  if(block_size.nx3 > 1) {
    nrootmesh=mesh_size.nx3*(1L<<(ll-root_level));
    if(block_size.x3rat == 1.0) { // uniform
      Real dx=(block_size.x3max-block_size.x3min)/block_size.nx3;
      for(int k=ks-NGHOST; k<=ke+NGHOST; ++k)
        dx3f(k)=dx;
      x3f(ks-NGHOST)=block_size.x3min-NGHOST*dx;
      for(int k=ks-NGHOST+1;k<=ke+NGHOST+1;k++)
        x3f(k)=x3f(k-1)+dx;
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
    }
    else {
      for (int k=ks-NGHOST; k<=ke+NGHOST+1; ++k) {
         // if there are too many levels, this won't work or be precise enough
        noffset=(k-ks)+(long long)lx3*block_size.nx3;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x3f(k)=pm->MeshGeneratorX3(rx,mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
      for(int k=ks-NGHOST; k<=ke+NGHOST; ++k)
        dx3f(k)=x3f(k+1)-x3f(k);
    }

  // correct cell face positions in ghost zones for reflecting boundary condition
    if (pmy_block->block_bcs[inner_x3] == 1) {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }
    if (pmy_block->block_bcs[outer_x3] == 1) {
      for (int k=1; k<=(NGHOST); ++k) {
        dx3f(ke+k  ) = dx3f(ke-k+1);
         x3f(ke+k+1) =  x3f(ke+k) + dx3f(ke+k);
      }
    }
  }
  else {
    dx3f(ks) = block_size.x3max-block_size.x3min;
    x3f(ks  ) = block_size.x3min;
    x3f(ke+1) = block_size.x3max;
  }

  if(pm->multilevel==true) {// set coarse coordinates for SMR/AMR
    int cnghost=pmy_block->cnghost;
    int cis=pmy_block->cis, cjs=pmy_block->cjs, cks=pmy_block->cks;
    int cie=pmy_block->cie, cje=pmy_block->cje, cke=pmy_block->cke;

    // x1
    for(int i=cis; i<=cie; i++) { // active region
      int ifl=(i-cis)*2+is;
      coarse_x1f(i)=x1f(ifl);
    }
    coarse_x1f(cie+1)=x1f(ie+1);
    // left ghost zone
    if(pmy_block->block_bcs[inner_x1]==1) { // reflecting
      for(int i=1; i<=cnghost; i++) {
        Real dx1=coarse_x1f(cis+i)-coarse_x1f(cis+i-1);
        coarse_x1f(cis-i) =  coarse_x1f(cis-i+1) - dx1;
      }
    }
    else {
      if(block_size.x1rat==1.0) { // uniform
        for(int i=1; i<=cnghost; i++)
          coarse_x1f(cis-i) =  2.0* coarse_x1f(cis-i+1) - coarse_x1f(cis-i+2);
      }
      else {
        for(int i=1; i<=cnghost; i++) { 
          noffset=(long long)lx1*block_size.nx1-i*2;
          Real rx=(Real)noffset/(Real)nrootmesh;
          coarse_x1f(cis-i)=pm->MeshGeneratorX1(rx,pm->mesh_size);
        }
      }
    }
    // right ghost zone
    if(pmy_block->block_bcs[outer_x1]==1) { // reflecting
      for(int i=1; i<=cnghost; i++) {
        Real dx1=coarse_x1f(cie-i+2)-coarse_x1f(cie-i+1);
        coarse_x1f(cie+i+1) =  coarse_x1f(cie+i) + dx1;
      }
    }
    else {
      if(block_size.x1rat==1.0) { // uniform
        for(int i=1; i<=cnghost; i++)
          coarse_x1f(cie+i+1) = 2.0*coarse_x1f(cie+i) - coarse_x1f(cie+i-1);
      }
      else {
        for(int i=1; i<=cnghost; i++) { 
          noffset=(long long)(lx1+1L)*block_size.nx1+i*2;
          Real rx=(Real)noffset/(Real)nrootmesh;
          coarse_x1f(cie+i+1)=pm->MeshGeneratorX1(rx,pm->mesh_size);
        }
      }
    }

    // x2
    for(int j=cjs; j<=cje; j++) { // active region
      int jfl=(j-cjs)*2+js;
      coarse_x2f(j)=x2f(jfl);
    }
    coarse_x2f(cje+1)=x2f(je+1);
    // left ghost zone
    if(block_size.nx2 > 1) {
      if(pmy_block->block_bcs[inner_x2]==1) { // reflecting
        for(int j=1; j<=cnghost; j++) {
          Real dx2=coarse_x2f(cjs+j)-coarse_x2f(cjs+j-1);
          coarse_x2f(cjs-j) =  coarse_x2f(cjs-j+1) - dx2;
        }
      }
      else {
        if(block_size.x2rat==1.0) { // uniform
          for(int j=1; j<=cnghost; j++)
            coarse_x2f(cjs-j) =  2.0* coarse_x2f(cjs-j+1) - coarse_x2f(cjs-j+2);
        }
        else {
          for(int j=1; j<=cnghost; j++) { 
            noffset=(long long)lx2*block_size.nx2-j*2;
            Real rx=(Real)noffset/(Real)nrootmesh;
            coarse_x2f(cjs-j)=pm->MeshGeneratorX2(rx,pm->mesh_size);
          }
        }
      }
      // right ghost zone
      if(pmy_block->block_bcs[outer_x2]==1) { // reflecting
        for(int j=1; j<=cnghost; j++) {
          Real dx2=coarse_x2f(cje-j+2)-coarse_x2f(cje-j+1);
          coarse_x2f(cje+j+1) =  coarse_x2f(cje+j) + dx2;
        }
      }
      else {
        if(block_size.x2rat==1.0) { // uniform
          for(int j=1; j<=cnghost; j++)
            coarse_x2f(cje+j+1) = 2.0*coarse_x2f(cje+j) - coarse_x2f(cje+j-1);
        }
        else {
          for(int j=1; j<=cnghost; j++) { 
            noffset=(long long)(lx2+1L)*block_size.nx2+j*2;
            Real rx=(Real)noffset/(Real)nrootmesh;
            coarse_x2f(cje+j+1)=pm->MeshGeneratorX2(rx,pm->mesh_size);
          }
        }
      }
    }

    // x3
    for(int k=cks; k<=cke; k++) { // active region
      int kfl=(k-cks)*2+ks;
      coarse_x3f(k)=x3f(kfl);
    }
    coarse_x3f(cke+1)=x3f(ke+1);
    if(block_size.nx3 > 1) {
      // left ghost zone
      if(pmy_block->block_bcs[inner_x3]==1) { // reflecting
        for(int k=1; k<=cnghost; k++) {
          Real dx3=coarse_x3f(cks+k)-coarse_x3f(cks+k-1);
          coarse_x3f(cks-k) =  coarse_x3f(cks-k+1) - dx3;
        }
      }
      else if(block_size.nx3 > 1) {
        if(block_size.x3rat==1.0) { // uniform
          for(int k=1; k<=cnghost; k++)
            coarse_x3f(cks-k) =  2.0* coarse_x3f(cks-k+1) - coarse_x3f(cks-k+2);
        }
        else {
          for(int k=1; k<=cnghost; k++) { 
            noffset=(long long)lx3*block_size.nx3-k*2;
            Real rx=(Real)noffset/(Real)nrootmesh;
            coarse_x3f(cks-k)=pm->MeshGeneratorX3(rx,pm->mesh_size);
          }
        }
      }
      // right ghost zone
      if(pmy_block->block_bcs[outer_x3]==1) { // reflecting
        for(int k=1; k<=cnghost; k++) {
          Real dx3=coarse_x3f(cke-k+2)-coarse_x3f(cke-k+1);
          coarse_x3f(cke+k+1) =  coarse_x3f(cke+k) + dx3;
        }
      }
      else {
        if(block_size.x3rat==1.0) { // uniform
          for(int k=1; k<=cnghost; k++)
            coarse_x3f(cke+k+1) = 2.0*coarse_x3f(cke+k) - coarse_x3f(cke+k-1);
        }
        else {
          for(int k=1; k<=cnghost; k++) { 
            noffset=(long long)(lx3+1L)*block_size.nx3+k*2;
            Real rx=(Real)noffset/(Real)nrootmesh;
            coarse_x3f(cke+k+1)=pm->MeshGeneratorX3(rx,pm->mesh_size);
          }
        }
      }
    }
  }
  return;
}


inline void Coordinates::DeleteBasicCoordinates(void)
{
  //destroy AthenaArrays
  dx1f.DeleteAthenaArray();
  dx2f.DeleteAthenaArray();
  dx3f.DeleteAthenaArray();
  dx1v.DeleteAthenaArray();
  dx2v.DeleteAthenaArray();
  dx3v.DeleteAthenaArray();
  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  if(pmy_block->pmy_mesh->multilevel==true) {
    coarse_x1f.DeleteAthenaArray();
    coarse_x2f.DeleteAthenaArray();
    coarse_x3f.DeleteAthenaArray();
    coarse_x1v.DeleteAthenaArray();
    coarse_x2v.DeleteAthenaArray();
    coarse_x3v.DeleteAthenaArray();

    if (MAGNETIC_FIELDS_ENABLED) {
      x1s2.DeleteAthenaArray();
      x1s3.DeleteAthenaArray();
      x2s1.DeleteAthenaArray();
      x2s3.DeleteAthenaArray();
      x3s1.DeleteAthenaArray();
      x3s2.DeleteAthenaArray();
      coarse_x1s2.DeleteAthenaArray();
      coarse_x1s3.DeleteAthenaArray();
      coarse_x2s1.DeleteAthenaArray();
      coarse_x2s3.DeleteAthenaArray();
      coarse_x3s1.DeleteAthenaArray();
      coarse_x3s2.DeleteAthenaArray();
    }
  }
}

#endif // COORDINATES_HPP
