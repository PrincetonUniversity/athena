#ifndef PLIMPTON_PACK_3D
#define PLIMPTON_PACK_3D

/* parallel pack functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

/* loop counters for doing a pack/unpack */

struct pack_plan_2d {
  int nfast;                 /* # of elements in fast index */
  int nslow;                 /* # of elements in slow index */
  int nstride;               /* stride between succesive slow indices */
  int nqty;                  /* # of values/element */
};

/* function prototypes */

void pack_2d(double *, double *, struct pack_plan_2d *);
void unpack_2d(double *, double *, struct pack_plan_2d *);
void unpack_2d_permute_1(double *, double *, struct pack_plan_2d *);
void unpack_2d_permute_2(double *, double *, struct pack_plan_2d *);
void unpack_2d_permute_n(double *, double *, struct pack_plan_2d *);

#endif
