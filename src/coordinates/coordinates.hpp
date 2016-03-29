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
#include <iostream>

// forward declarations
class MeshBlock;
class ParameterInput;

//! \class Coordinates
//  \brief coordinate data and functions

class Coordinates {
public:
  friend class HydroSourceTerms;
  Coordinates(MeshBlock *pmb, ParameterInput *pin, int flag = 0);
  ~Coordinates();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Coordinates

  // coordinate arrays
  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions

  // arrays for MHD mesh refinement
  //note: only x1s2 and x1s3 in spherical and schwarzschild coordinates are non trivial
  AthenaArray<Real> x1s2, x1s3, x2s1, x2s3, x3s1, x3s2; // area averaged positions

  void AllocateAndSetBasicCoordinates(void);
  void DeleteBasicCoordinates(void);
  void CheckMeshSpacing(void);

  // Function for checking for poles
  bool IsPole(int j);

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
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

// Functions to calculate covariant derivatives at faces, for viscosity calculations  
  #if !GENERAL_RELATIVITY  // declare, but do not define, in non-GR case
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
  #else  // GR: define as no-op in all coordinate systems
    void FaceXdx(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceXdy(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceXdz(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceYdx(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceYdy(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceYdz(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceZdx(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceZdy(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
    void FaceZdz(const int, const int, const int, const int, const AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}
  #endif  // !GENERAL_RELATIVITY

// function to compute Divv
  #if !GENERAL_RELATIVITY  // declare, but do not define, in non-GR case
    void Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv);
  #else  // GR: define as no-op in all coordinate systems
    void Divv(const AthenaArray<Real> &, AthenaArray<Real> &) {return;}
  #endif  // !GENERAL_RELATIVITY

// function to compute viscous source terms
  #if !GENERAL_RELATIVITY  // declare, but do not define, in non-GR case
    void VisSrcTermsX1(const int k, const int j, const Real dt,
      const AthenaArray<Real> &flx,
      const AthenaArray<Real> &prim, AthenaArray<Real> &u);
    void VisSrcTermsX2(const int k, const int j, const Real dt,
      const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
      const AthenaArray<Real> &prim, AthenaArray<Real> &u);
    void VisSrcTermsX3(const int k, const int j, const Real dt,
      const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
      const AthenaArray<Real> &prim, AthenaArray<Real> &u);
  #else  // GR: define as no-op in all coordinate systems
    void VisSrcTermsX1(const int, const int, const Real, const AthenaArray<Real> &,
        const AthenaArray<Real> &, AthenaArray<Real> &) {return;}
    void VisSrcTermsX2(const int, const int, const Real, const AthenaArray<Real> &,
        const AthenaArray<Real> &, const AthenaArray<Real> &, AthenaArray<Real> &)
        {return;}
    void VisSrcTermsX3(const int, const int, const Real, const AthenaArray<Real> &,
        const AthenaArray<Real> &, const AthenaArray<Real> &, AthenaArray<Real> &)
        {return;}
  #endif  // !GENERAL_RELATIVITY

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
    void RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j, int i,
        Real *pa0, Real *pa1, Real *pa2, Real *pa3);
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
    void RaiseVectorCell(Real, Real, Real, Real, int, int, int, Real *, Real *, Real *,
        Real *) {return;}
    void LowerVectorCell(Real, Real, Real, Real, int, int, int, Real *, Real *, Real *,
        Real *) {return;}
  #endif  // GENERAL_RELATIVITY

private:

  int cflag;
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
  AthenaArray<Real> coord_src3_j_;

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
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // allocate arrays for sizes and positions of cells
  int is, ie, js, je, ks, ke, ng;
  if(cflag==0) {
    is = pmy_block->is; js = pmy_block->js; ks = pmy_block->ks;
    ie = pmy_block->ie; je = pmy_block->je; ke = pmy_block->ke;
    ng=NGHOST;
  }
  else {
    is = pmy_block->cis; js = pmy_block->cjs; ks = pmy_block->cks;
    ie = pmy_block->cie; je = pmy_block->cje; ke = pmy_block->cke;
    ng=pmy_block->cnghost;
  }
  int ncells1 = (ie-is+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (je-js+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ke-ks+1) + 2*ng;

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

  if((pm->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    // area weighted positions for AMR/SMR MHD
    x1s2.NewAthenaArray(ncells1);
    x1s3.NewAthenaArray(ncells1);
    x2s1.NewAthenaArray(ncells2);
    x2s3.NewAthenaArray(ncells2);
    x3s1.NewAthenaArray(ncells3);
    x3s2.NewAthenaArray(ncells3);
  }

  long long nrootmesh, noffset;
  long int &lx1=pmy_block->loc.lx1;
  long int &lx2=pmy_block->loc.lx2;
  long int &lx3=pmy_block->loc.lx3;
  int &ll=pmy_block->loc.level;

  // X1-DIRECTION: initialize sizes and positions of cell FACES (dx1f,x1f)
  nrootmesh=mesh_size.nx1*(1L<<(ll-pm->root_level));
  if(pm->user_meshgen_[x1dir]==false) { // uniform
    Real dx=(block_size.x1max-block_size.x1min)/(ie-is+1);
    for(int i=is-ng; i<=ie+ng; ++i)
      dx1f(i)=dx;
    x1f(is-ng)=block_size.x1min-ng*dx;
    for(int i=is-ng+1;i<=ie+ng+1;i++)
      x1f(i)=x1f(i-1)+dx;
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
  }
  else {
    for (int i=is-ng; i<=ie+ng+1; ++i) {
      // if there are too many levels, this won't work or be precise enough
      noffset=((i-is)<<cflag)+(long long)lx1*block_size.nx1;
      Real rx=(Real)noffset/(Real)nrootmesh;
      x1f(i)=pm->MeshGenerator_[x1dir](rx,mesh_size);
    }
    x1f(is) = block_size.x1min;
    x1f(ie+1) = block_size.x1max;
    for(int i=is-ng; i<=ie+ng; ++i)
      dx1f(i)=x1f(i+1)-x1f(i);
  }

  // correct cell face positions in ghost zones for reflecting boundary condition
  if (pmy_block->block_bcs[INNER_X1] == 1) {
    for (int i=1; i<=ng; ++i) {
      dx1f(is-i) = dx1f(is+i-1);
       x1f(is-i) =  x1f(is-i+1) - dx1f(is-i);
    }
  }
  if (pmy_block->block_bcs[OUTER_X1] == 1) {
    for (int i=1; i<=ng; ++i) {
      dx1f(ie+i  ) = dx1f(ie-i+1);
       x1f(ie+i+1) =  x1f(ie+i) + dx1f(ie+i);
    }
  }

  // X2-DIRECTION: initialize spacing and positions of cell FACES (dx2f,x2f)
  if(ncells2 > 1) {
    nrootmesh=mesh_size.nx2*(1L<<(ll-pm->root_level));
    if(pm->user_meshgen_[x2dir]==false) { // uniform
      Real dx=(block_size.x2max-block_size.x2min)/(je-js+1);
      for(int j=js-ng; j<=je+ng; ++j)
        dx2f(j)=dx;
      x2f(js-ng)=block_size.x2min-ng*dx;
      for(int j=js-ng+1;j<=je+ng+1;j++)
        x2f(j)=x2f(j-1)+dx;
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
    }
    else {
      for (int j=js-ng; j<=je+ng+1; ++j) {
        // if there are too many levels, this won't work or be precise enough
        noffset=((j-js)<<cflag)+(long long)lx2*block_size.nx2;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x2f(j)=pm->MeshGenerator_[x2dir](rx,mesh_size);
      }
      x2f(js) = block_size.x2min;
      x2f(je+1) = block_size.x2max;
      for(int j=js-ng; j<=je+ng; ++j)
        dx2f(j)=x2f(j+1)-x2f(j);
    }

    // correct cell face positions in ghost zones for reflecting boundary condition
    if (pmy_block->block_bcs[INNER_X2] == 1) {
      for (int j=1; j<=ng; ++j) {
        dx2f(js-j) = dx2f(js+j-1);
         x2f(js-j) =  x2f(js-j+1) - dx2f(js-j);
      }
    }
    if (pmy_block->block_bcs[OUTER_X2] == 1) {
      for (int j=1; j<=ng; ++j) {
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
  if(ncells3 > 1) {
    nrootmesh=mesh_size.nx3*(1L<<(ll-pm->root_level));
    if(pm->user_meshgen_[x3dir]==false) { // uniform
      Real dx=(block_size.x3max-block_size.x3min)/(ke-ks+1);
      for(int k=ks-ng; k<=ke+ng; ++k)
        dx3f(k)=dx;
      x3f(ks-ng)=block_size.x3min-ng*dx;
      for(int k=ks-ng+1;k<=ke+ng+1;k++)
        x3f(k)=x3f(k-1)+dx;
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
    }
    else {
      for (int k=ks-ng; k<=ke+ng+1; ++k) {
        // if there are too many levels, this won't work or be precise enough
        noffset=((k-ks)<<cflag)+(long long)lx3*block_size.nx3;
        Real rx=(Real)noffset/(Real)nrootmesh;
        x3f(k)=pm->MeshGenerator_[x3dir](rx,mesh_size);
      }
      x3f(ks) = block_size.x3min;
      x3f(ke+1) = block_size.x3max;
      for(int k=ks-ng; k<=ke+ng; ++k)
        dx3f(k)=x3f(k+1)-x3f(k);
    }

    // correct cell face positions in ghost zones for reflecting boundary condition
    if (pmy_block->block_bcs[INNER_X3] == 1) {
      for (int k=1; k<=ng; ++k) {
        dx3f(ks-k) = dx3f(ks+k-1);
         x3f(ks-k) =  x3f(ks-k+1) - dx3f(ks-k);
      }
    }
    if (pmy_block->block_bcs[OUTER_X3] == 1) {
      for (int k=1; k<=ng; ++k) {
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
  if((pmy_block->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    x1s2.DeleteAthenaArray();
    x1s3.DeleteAthenaArray();
    x2s1.DeleteAthenaArray();
    x2s3.DeleteAthenaArray();
    x3s1.DeleteAthenaArray();
    x3s2.DeleteAthenaArray();
  }
}


inline void Coordinates::CheckMeshSpacing(void)
{
  Real rmax=1.0, rmin=1.0;
  for(int i=pmy_block->is; i<=pmy_block->ie; i++) {
    rmax=std::max(dx1v(i+1)/dx1v(i),rmax);
    rmin=std::min(dx1v(i+1)/dx1v(i),rmin);
  }
  if(rmax > 1.1 || rmin  < 1.0/1.1) {
     std::cout << "### Warning in Coordinates::CheckMeshSpacing" << std::endl
       << "Neighboring cell sizes differ by more than 10% in the x1 direction."
       << std::endl;
  }
  if(pmy_block->je!=pmy_block->js) {
    rmax=1.0, rmin=1.0;
    for(int j=pmy_block->js; j<=pmy_block->je; j++) {
      rmax=std::max(dx2v(j+1)/dx2v(j),rmax);
      rmin=std::min(dx2v(j+1)/dx2v(j),rmin);
    }
    if(rmax > 1.1 || rmin  < 1.0/1.1) {
       std::cout << "### Warning in Coordinates::CheckMeshSpacing" << std::endl
         << "Neighboring cell sizes differ by more than 10% in the x2 direction."
         << std::endl;
    }
  }
  if(pmy_block->ke!=pmy_block->ks) {
    rmax=1.0, rmin=1.0;
    for(int k=pmy_block->ks; k<=pmy_block->ke; k++) {
      rmax=std::max(dx3v(k+1)/dx3v(k),rmax);
      rmin=std::min(dx3v(k+1)/dx3v(k),rmin);
    }
    if(rmax > 1.1 || rmin  < 1.0/1.1) {
       std::cout << "### Warning in Coordinates::CheckMeshSpacing" << std::endl
         << "Neighboring cell sizes differ by more than 10% in the x3 direction."
         << std::endl;
    }
  }
  return;
}

// Function for determining if index corresponds to a polar boundary
// Inputs:
//   j: x2-index
// Outputs:
//   returned value: true if face indexed with j is on a pole; false otherwise
inline bool Coordinates::IsPole(int j)
{
  if (pmy_block->block_bcs[INNER_X2] == POLAR_BNDRY and j == pmy_block->js)
    return true;
  if (pmy_block->block_bcs[OUTER_X2] == POLAR_BNDRY and j == pmy_block->je+1)
    return true;
  return false;
}

#endif // COORDINATES_HPP
