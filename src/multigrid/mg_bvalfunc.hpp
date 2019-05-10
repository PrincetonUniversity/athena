#ifndef MG_BVALFUNC_HPP_
#define MG_BCALFUNC_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_bvalfunc.hpp
//  \brief defines the boundary functions for Multigrid

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

void MGPeriodicInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

void MGZeroGradientInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

void MGZeroFixedInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

#endif // MG_BVALFUNC_HPP_
