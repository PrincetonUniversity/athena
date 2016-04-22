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
//! \file default_pgen.cpp
//  \brief Provides default (empty) versions of all functions in problem generator files
//  This means user does not have to implement these functions if they are not needed.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//======================================================================================

void __attribute__((weak)) Mesh::InitUserMeshData(ParameterInput *pin)
{
  // do nothing
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Should be used to set initial conditions.
//======================================================================================

void __attribute__((weak)) MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // In practice, this function should *always* be replaced by a version 
  // that sets the initial conditions for the problem of interest.
  return;
}

//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//======================================================================================

void __attribute__((weak)) MeshBlock::UserWorkInLoop(void)
{
  // do nothing
  return;
}

//======================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Function called after main loop is finished for user-defined work.
//======================================================================================

void __attribute__((weak)) Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  // do nothing
  return;
}
