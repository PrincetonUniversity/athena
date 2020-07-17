#ifndef SPHERICAL_GRID_HPP_
#define SPHERICAL_GRID_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_grid.hpp
//  \brief defines SphericalGrid and SphericalPatch

#include <vector>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../lagrange_interp.hpp"
#include "mesh.hpp"
#include "geodesic_grid.hpp"

//! \class SphericalGrid
//! \brief A class representing a grid on a topological sphere (wrapping around GeodesicGrid)
class SphericalGrid: public GeodesicGrid {
  public:
    //! Creates a geodetic grid with nlev levels
    SphericalGrid(
      int nlev,               //! [in] number of levels of the grid
      Real rad = 1.0          //! [in] radius
    );

    //! Returns the x,y,z coordinates of a vertex on the primal grid
    void Position(
      int ic,                 //! [in] cell number
      Real * x,               //! [out] x-position
      Real * y,               //! [out] y-position
      Real * z                //! [out] z-position
    ) const;
    //! Get the surface area of a given cell in the dual grid
    Real ComputeWeight(
      int ic                  //! [in] cell index
    ) const;
    //! Get the arc length of the face between cells ic1 and ic2
    /*!
     *  This is computed only to 2nd order accuracy at the moment.
     *
     *  An assert error is thrown if ic1 and ic2 are not neighbors or if either
     *  index is invalid.
     */
    Real ArcLength(
      int ic1,                //! [in] index of first cell
      int ic2                 //! [in] index of second cell
    ) const;
  public:
    AthenaArray<Real> rad;
};

//! \class SphericalPatch
//! \brief This class represents the intersection between a spherical grid and a mesh block
//!  Note: this class assumes that the MeshBlock is Cartesian uniformly spaced
class SphericalPatch {
  public:
    enum collocation_t {cell, vertex};
  public:
    SphericalPatch(SphericalGrid const * psphere, MeshBlock const * pblock, collocation_t coll);
    ~SphericalPatch();
    //! Interpolate a group of arrays defined on the MeshBlock to the SphericalPatch
    //  The destination array should be allocated
    void InterpToSpherical(AthenaArray<Real> const & src,    //! [in] data defined on the MeshBlock
                           AthenaArray<Real> * dst) const;   //! [out] data defined on the SphericalGrid
    //! Merge data into a global array defined on the whole sphere
    void MergeData(AthenaArray<Real> const & src,            //! [in] data defined on the SphericalPatch
                   AthenaArray<Real> * dst) const;           //! [out] data defined on the SphericalGrid
    //! Number of points on the spherical patch
    inline int NumPoints() const {
      return n;
    }
    //! Map patch degrees of freedom to the corresponding index in the full SphericalGrid
    inline int idxMap(int idx) const {
      return map[idx];
    }
  private:
    //! Interpolate an array defined on the MeshBlock to the SphericalPatch
    //  The destination array should be allocated
    void interpToSpherical(Real const * src,    //! [in] 1D data defined on the MeshBlock
                           Real * dst) const;   //! [out] 1D data defined on the SphericalGrid
    //! Merge data arrays into a global arrays defined on the whole sphere
    void mergeData(Real const * src,            //! [in] 1D data defined on the SphericalPatch
                   Real * dst) const;           //! [out] 1D data defined on the SphericalGrid
  public:
    //! Type of collocation
    collocation_t const coll;
    //! Parent spherical grid
    SphericalGrid const * psphere;
    // Parent mesh block
    MeshBlock const * pblock;
  private:
    // Number of points in the spherical patch
    int n;
    // Maps local indices to global indices on the SphericalGrid
    std::vector<int> map;
    // Interpolating polynomials
    LagrangeInterpND<2*NGHOST-1, 3> ** pinterp;
};

#endif
