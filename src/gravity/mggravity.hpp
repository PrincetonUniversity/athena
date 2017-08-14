#ifndef MGGRAVITY_HPP
#define MGGRAVITY_HPP

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
  MGGravity(Mesh *pm, MeshBlock *pmb, int nx, int ny, int nz,
            RegionSize isize, MGBoundaryFunc_t *MGBoundary)
   : Multigrid(pm,pmb,1,nx,ny,nz,1,isize,MGBoundary),omega_(1.15)
  { btype=BNDRY_MGGRAV; btypef=BNDRY_MGGRAVF; };
  ~MGGravity() {};
  void Smooth(int color);
  void CalculateDefect(void);

private:
  const Real omega_;
};


//! \class GravityDriver
//  \brief Multigrid gravity solver

class GravityDriver : public MultigridDriver{
public:
  GravityDriver(Mesh *pm, MeshBlock *pmb, MGBoundaryFunc_t *MGBoundary,
                ParameterInput *pin);
  Multigrid* GetMultigridBlock (MeshBlock *pmb);
  void LoadSourceAndData(void);

private:
  Real four_pi_G_;
};

#endif // MGGRAVITY_HPP
