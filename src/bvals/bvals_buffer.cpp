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
//! \file bvals_buffer.cpp
//  \brief utility functions for BoundaryValues buffers
//======================================================================================

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/buffer_utils.hpp"

// this class header
#include "bvals.hpp"

//--------------------------------------------------------------------------------------
//! \fn unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
//  \brief calculate a buffer identifier
unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
{
  unsigned int ux1=(unsigned)(ox1+1);
  unsigned int ux2=(unsigned)(ox2+1);
  unsigned int ux3=(unsigned)(ox3+1);
  return (ux1<<6) | (ux2<<4) | (ux3<<2) | (fi1<<1) | fi2;
}


//--------------------------------------------------------------------------------------
//! \fn int BufferID(int dim, bool multilevel)
//  \brief calculate neighbor indexes and target buffer IDs
int BufferID(int dim, bool multilevel)
{
  int nf1=1, nf2=1;
  if(multilevel==true) {
    if(dim>=2) nf1=2;
    if(dim>=3) nf2=2;
  }
  int b=0;
  // x1 face
  for(int n=-1; n<=1; n+=2) {
    for(int f2=0;f2<nf2;f2++) {
      for(int f1=0;f1<nf1;f1++) {
        BoundaryValues::ni[b].ox1=n;
        BoundaryValues::ni[b].ox2=0;
        BoundaryValues::ni[b].ox3=0;
        BoundaryValues::ni[b].fi1=f1;
        BoundaryValues::ni[b].fi2=f2;
        BoundaryValues::ni[b].type=NEIGHBOR_FACE;
        b++;
      }
    }
  }
  // x2 face
  if(dim>=2) {
    for(int n=-1; n<=1; n+=2) {
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          BoundaryValues::ni[b].ox1=0;
          BoundaryValues::ni[b].ox2=n;
          BoundaryValues::ni[b].ox3=0;
          BoundaryValues::ni[b].fi1=f1;
          BoundaryValues::ni[b].fi2=f2;
          BoundaryValues::ni[b].type=NEIGHBOR_FACE;
          b++;
        }
      }
    }
  }
  if(dim==3) {
    // x3 face
    for(int n=-1; n<=1; n+=2) {
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          BoundaryValues::ni[b].ox1=0;
          BoundaryValues::ni[b].ox2=0;
          BoundaryValues::ni[b].ox3=n;
          BoundaryValues::ni[b].fi1=f1;
          BoundaryValues::ni[b].fi2=f2;
          BoundaryValues::ni[b].type=NEIGHBOR_FACE;
          b++;
        }
      }
    }
  }
  // edges
  // x1x2
  if(dim>=2) {
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf2;f1++) {
          BoundaryValues::ni[b].ox1=n;
          BoundaryValues::ni[b].ox2=m;
          BoundaryValues::ni[b].ox3=0;
          BoundaryValues::ni[b].fi1=f1;
          BoundaryValues::ni[b].fi2=0;
          BoundaryValues::ni[b].type=NEIGHBOR_EDGE;
          b++;
        }
      }
    }
  }
  if(dim==3) {
    // x1x3
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf1;f1++) {
          BoundaryValues::ni[b].ox1=n;
          BoundaryValues::ni[b].ox2=0;
          BoundaryValues::ni[b].ox3=m;
          BoundaryValues::ni[b].fi1=f1;
          BoundaryValues::ni[b].fi2=0;
          BoundaryValues::ni[b].type=NEIGHBOR_EDGE;
          b++;
        }
      }
    }
    // x2x3
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        for(int f1=0;f1<nf1;f1++) {
          BoundaryValues::ni[b].ox1=0;
          BoundaryValues::ni[b].ox2=n;
          BoundaryValues::ni[b].ox3=m;
          BoundaryValues::ni[b].fi1=f1;
          BoundaryValues::ni[b].fi2=0;
          BoundaryValues::ni[b].type=NEIGHBOR_EDGE;
          b++;
        }
      }
    }
    // corners
    for(int l=-1; l<=1; l+=2) {
      for(int m=-1; m<=1; m+=2) {
        for(int n=-1; n<=1; n+=2) {
          BoundaryValues::ni[b].ox1=n;
          BoundaryValues::ni[b].ox2=m;
          BoundaryValues::ni[b].ox3=l;
          BoundaryValues::ni[b].fi1=0;
          BoundaryValues::ni[b].fi2=0;
          BoundaryValues::ni[b].type=NEIGHBOR_CORNER;
          b++;
        }
      }
    }
  }

  for(int n=0;n<b;n++)
    BoundaryValues::bufid[n]=CreateBufferID(BoundaryValues::ni[n].ox1,
      BoundaryValues::ni[n].ox2, BoundaryValues::ni[n].ox3, BoundaryValues::ni[n].fi1,
      BoundaryValues::ni[n].fi2);

  return b;
}

//--------------------------------------------------------------------------------------
//! \fn int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax)
//  \brief
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax)
{
  int bid=CreateBufferID(ox1, ox2, ox3, fi1, fi2);

  for(int i=0;i<bmax;i++) {
    if(bid==BoundaryValues::bufid[i]) return i;
  }
  return -1;
}

//--------------------------------------------------------------------------------------
//! \fn unsigned int CreateBvalsMPITag(int lid, int phys, int bufid)
//  \brief calculate an MPI tag for Bval communications
// tag = local id of destination (19) + bufid(7) + physics(5)

unsigned int CreateBvalsMPITag(int lid, int phys, int bufid)
{
  return (lid<<12) | (bufid<<5) | phys;
}
