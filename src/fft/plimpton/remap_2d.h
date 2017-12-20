#ifndef PLIMPTON_REMAP_2D
#define PLIMPTON_REMAP_2D

/* parallel remap functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

/* details of how to do a 2d remap */

struct remap_plan_2d {
  double *sendbuf;                  /* buffer for MPI sends */
  double *scratch;                  /* scratch buffer for MPI recvs */
  void (*pack)(double *, double *, struct pack_plan_2d *);                   /* which pack function to use */
  void (*unpack)(double *, double *, struct pack_plan_2d *);                 /* which unpack function to use */
  int *send_offset;                 /* extraction loc for each send */
  int *send_size;                   /* size of each send message */
  int *send_proc;                   /* proc to send each message to */
  struct pack_plan_2d *packplan;    /* pack plan for each send message */
  int *recv_offset;                 /* insertion loc for each recv */
  int *recv_size;                   /* size of each recv message */
  int *recv_proc;                   /* proc to recv each message from */
  int *recv_bufloc;                 /* offset in scratch buf for each recv */
  MPI_Request *request;             /* MPI request for each posted recv */
  struct pack_plan_2d *unpackplan;  /* unpack plan for each recv message */
  int nrecv;                        /* # of recvs from other procs */
  int nsend;                        /* # of sends to other procs */
  int self;                         /* whether I send/recv with myself */
  int memory;                       /* user provides scratch space or not */
  MPI_Comm comm;                    /* group of procs performing remap */
};

/* collision between 2 regions */

struct extent_2d {
  int ilo,ihi,isize;
  int jlo,jhi,jsize;
};

/* function prototypes */

void remap_2d(double *, double *, double *, struct remap_plan_2d *);
struct remap_plan_2d *remap_2d_create_plan(MPI_Comm,
  int, int, int, int, int, int, int, int,
  int, int, int, int);
void remap_2d_destroy_plan(struct remap_plan_2d *);
int remap_2d_collide(struct extent_2d *, 
		     struct extent_2d *, struct extent_2d *);

/* machine specifics */

#ifdef T3E_KLUDGE

#define remap_2d_ REMAP_2D
#define remap_2d_create_plan_ REMAP_2D_CREATE_PLAN
#define remap_2d_destroy_plan_ REMAP_2D_DESTROY_PLAN

#endif

#endif
