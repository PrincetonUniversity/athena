#ifndef MESH_GEODESIC_GRID_HPP_
#define MESH_GEODESIC_GRID_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_grid.hpp
//  \brief defines GeodesicGrid, a geodesic grid on the unit sphere

#include "../athena.hpp"                  // Real, indices, function prototypes
#include "../athena_arrays.hpp"           // AthenaArray

class GeodesicGrid {
  public:
    static const int EMPTY = -1;
    static const int NRGHOST = 1;

  public:
    //! Creates a geodesic grid with nlev levels
    //! The number of grid points is \f$ 10 {\rm nlev}^2 + 2\f$.
    GeodesicGrid(
        int nlev                 //! [in] number of levels of the grid
        );

    //! Returns the number of vertices in the primal grid
    int NumVertices() const;

    //! Returns the number of neighbors of a given cell in the dual grid
    /*!
     *  On output the array neighbors should contain the indices of all
     *  neighbors. If there are only 5 neighbors this method should return 5
     *  and set the last value to GeodesicGrid::EMPTY
     */
    int NumNeighbors(
        int ic,                  //! [in] cell number
        int neighbors[6]         //! [out] on output this arrays contains the indices of the neighbors
        ) const;

    //! Returns the theta and phi coordinate of a vertex on the primal grid
    void PositionPolar(
        int ic,                  //! [in] cell number
        Real * theta,            //! [out] latitude angle
        Real * phi               //! [out] azimuth angle
        ) const;

    //! Returns the x,y,z coordinates of a vertex on the primal grid
    void Position(
        int ic,                  //! [in] cell number
        Real * x,                //! [out] x-position
        Real * y,                //! [out] y-position
        Real * z                 //! [out] z-position
        ) const;

    //! Returns the x,y,z coordinates of a midpoint between vertices on a unit sphere
    void PositionMid(
        int ic1,                 //! [in] index of first cell
        int ic2,                 //! [in] index of second cell
        Real * x,                //! [out] x-position of midpoint
        Real * y,                //! [out] y-position of midpoint
        Real * z                 //! [out] z-position of midpoint
        ) const;

    //! Returns the components of a unit vector in the direction of the flux
    //! (it is tangent to the great circle connecting the two grid points
    //!  and originates at the midpoint between them on a unit sphere)
    void UnitFluxDir(
        int ic1,                //! [in] index of first cell
        int ic2,                //! [in] index of second cell
        Real * dtheta,          //! [out] theta component of the vector
        Real * dphi             //! [out] phi component of the vector
        ) const;

    //! Get the surface area of a given (hexagonal or pentagonal) cell in the
    //! dual grid and the lengths of its edges
    Real ComputeWeightAndDualEdges(
        int ic,                  //! [in] cell index
        Real length[6]           //! [out] lengths of the edges of the dual grid
        ) const;

    //! Get the surface area (weight) of the cell only
    Real ComputeWeight(
        int ic                   //! [in] cell index
        ) const;

    //! Returns true if two cells are neighbors, false otherwise
    /*!
     *  This method is used for debugging only and it is not very efficient
     */
    bool AreNeighbors(
        int ic1,                 //! [in] index of first cell
        int ic2                  //! [in] index of second cell
        ) const;

    //! Computes xi and eta coordinates for WENO scheme
    //! Equation (18) of Florinski et al. 2013 (arXiv:1302.2087v1)
    void ComputeXiEta(
        int ic,                  //! [in] cell index
        Real xi[6],
        Real eta[6]
        ) const;

    //! Get the arc length between the center of cells ic1 and ic2
    /*!
     *  An assert error is thrown if ic1 and ic2 are not neighbors or if either
     *  index is invalid.
     */
    Real ArcLength(
        int ic1,                 //! [in] index of first cell
        int ic2                  //! [in] index of second cell
        ) const;



    //! Finds the coordinates of the circumcenter of a triangle
    //! and projects it on a unit sphere
    void CircumcenterNormalized(Real x1, Real x2, Real x3, Real y1,
        Real y2, Real y3, Real z1, Real z2, Real z3,
        Real * x_cc, Real * y_cc, Real * z_cc) const;

    //! Finds the parameters of the great circle passing through two points
    void GreatCircleParam(Real zeta1, Real zeta2, Real psi1, Real psi2,
        Real * apar, Real * psi0) const;


  private:

    int numpoints_;
    int nlev_;

    AthenaArray<Real> geo_grid_xcomp; // x coordinate of the vertex on a unit sphere
    AthenaArray<Real> geo_grid_ycomp; // y coordinate of the vertex on a unit sphere
    AthenaArray<Real> geo_grid_zcomp; // z coordinate of the vertex on a unit sphere

    AthenaArray<Real> geo_grid_pol_xcomp; // x coordinate of the poles
    AthenaArray<Real> geo_grid_pol_ycomp; // y coordinate of the poles
    AthenaArray<Real> geo_grid_pol_zcomp; // z coordinate of the poles

    AthenaArray<int> geo_grid_index;
    AthenaArray<int> geo_grid_pol_index;


    //! Fills in the ghost regions of the blocks in geodesic grid
    template<typename T>
    void FillInGhostRegions(int nlev,AthenaArray<T> &blocks,
        AthenaArray<T> &poles) const {
      for (int bl = 0; bl < 5; ++bl){
        for (int k = 0; k < nlev; ++k){
          blocks(bl,0,k+1)           = blocks((bl+4)%5,k+1,1);
          blocks(bl,0,k+nlev+1)      = blocks((bl+4)%5,nlev,k+1);
          blocks(bl,k+1,2*nlev+1)    = blocks((bl+4)%5,nlev,k+nlev+1);
          blocks(bl,k+2,0)           = blocks((bl+1)%5,1,k+1);
          blocks(bl,nlev+1,k+1)      = blocks((bl+1)%5,1,k+nlev+1);
          blocks(bl,nlev+1,k+nlev+1) = blocks((bl+1)%5,k+2,2*nlev);
        }
        blocks(bl,1,0)           = poles(0);
        blocks(bl,nlev+1,2*nlev) = poles(1);
        blocks(bl,0,2*nlev+1)    = blocks(bl,0,2*nlev);
      }
    }

};

#endif
