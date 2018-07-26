//========================================================================================
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>    // memset
#include <iostream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "./multigrid.hpp"

//----------------------------------------------------------------------------------------
//! \fn Multigrid::Multigrid(MultigridDriver *pmd, LogicalLocation iloc, int igid,
//                           int ilid, int invar, int nghost, RegionSize isize,
//                           MGBoundaryFunc_t *MGBoundary, enum BoundaryFlag *input_bcs,
//                           bool root);
//  \brief Multigrid constructor

Multigrid::Multigrid(MultigridDriver *pmd, LogicalLocation iloc, int igid, int ilid,
           int invar, int nghost, RegionSize isize, MGBoundaryFunc_t *MGBoundary,
           enum BoundaryFlag *input_bcs, bool root) {
}


//----------------------------------------------------------------------------------------
//! \fn Multigrid::~Multigrid
//  \brief Multigrid destroctor

Multigrid::~Multigrid() {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh)
//  \brief Fill the inital guess in the active zone of the finest level
void Multigrid::LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh,
//                                 Real fac)
//  \brief Fill the source in the active zone of the finest level
void Multigrid::LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RestrictFMGSource(void)
//  \brief restrict the source through all the multigrid levels

void Multigrid::RestrictFMGSource(void) {
}

//----------------------------------------------------------------------------------------
//! \fn void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh)
//  \brief Set the result, including the ghost zone
void Multigrid::RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ZeroClearData(void)
//  \brief Clear the data array with zero
void Multigrid::ZeroClearData(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::Restrict(void)
//  \brief Restrict the defect to the source
void Multigrid::Restrict(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::ProlongateAndCorrect(void)
//  \brief Prolongate the potential using tri-linear interpolation
void Multigrid::ProlongateAndCorrect(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::FMGProlongate(void)
//  \brief Prolongate the potential for Full Multigrid cycle
void Multigrid::FMGProlongate(void) {
}


//----------------------------------------------------------------------------------------
//! \fn void Multigrid::SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck)
//  \brief Load the data from the root grid
void Multigrid::SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck) {
}




//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the residual norm

Real Multigrid::CalculateDefectNorm(int n, int nrm) {
  return 0.0;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::CalculateTotal(int type, int n)
//  \brief calculate the sum of the array (type: 0=src, 1=u)

Real Multigrid::CalculateTotal(int type, int n) {
  return 0.0;
}


//----------------------------------------------------------------------------------------
//! \fn Real Multigrid::SubtractAverage(int type, int n, Real ave)
//  \brief subtract the average value (type: 0=source, 1=u)

void Multigrid::SubtractAverage(int type, int n, Real ave) {
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X1 direction

void MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X1 direction

void MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X2 direction

void MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X2 direction

void MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the inner-X3 direction

void MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
}


//----------------------------------------------------------------------------------------
//! \fn MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, int nvar,
//                int is, int ie, int js, int je, int ks, int ke, int ngh,
//                Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief Periodic (default) boundary condition in the outer-X3 direction

void MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
}
