#ifndef GRAVITY_MGGRAVITY_HPP_
#define GRAVITY_MGGRAVITY_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mggravity.hpp
//  \brief defines MGGravity class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../multigrid/multigrid.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Multigrid;

//! \class MGGravity
//  \brief Multigrid gravity solver for each block

class MGGravity : public Multigrid {
public:
  MGGravity(MultigridDriver *pmd, LogicalLocation iloc, int igid, int ilid,
            RegionSize isize, MGBoundaryFunc_t *MGBoundary,
            enum BoundaryFlag *input_bcs, bool root)
   : Multigrid(pmd,iloc,igid,ilid,1,1,isize,MGBoundary,input_bcs,root), omega_(1.15)
  { btype=BNDRY_MGGRAV; btypef=BNDRY_MGGRAVF; };
  ~MGGravity() {}
  void Smooth(int color);
  void CalculateDefect(void);

private:
  const Real omega_;
};


//! \class MGGravityDriver
//  \brief Multigrid gravity solver

class MGGravityDriver : public MultigridDriver{
public:
  MGGravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary, ParameterInput *pin);
  ~MGGravityDriver() {}
  void Solve(int stage);

private:
  Real four_pi_G_;
};

#endif // GRAVITY_MGGRAVITY_HPP_
