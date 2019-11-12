//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file history.cpp
//  \brief writes history output data, volume-averaged quantities that are output
//         frequently in time to trace their history.

// C headers

// C++ headers
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
// BD: new problem
#include "../wave/wave.hpp"
// -BD
#include "../advection/advection.hpp"
#include "../z4c/z4c.hpp"
#include "../mesh/mesh.hpp"
#include "../scalars/scalars.hpp"
#include "outputs.hpp"

// NEW_OUTPUT_TYPES:

// "3" for 1-KE, 2-KE, 3-KE additional columns (come before tot-E)
// BD: new problem
#define NHISTORY_VARS (((NHYDRO) + 3) * (FLUID_ENABLED) + (SELF_GRAVITY_ENABLED) + \
                       (NFIELD) + (NSCALARS) + \
                       2 * (WAVE_ENABLED) + \
                       2 * (ADVECTION_ENABLED))
// -BD

//----------------------------------------------------------------------------------------
//! \fn void OutputType::HistoryFile()
//  \brief Writes a history file

void HistoryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  MeshBlock *pmb = pm->pblock;
  Real real_max = std::numeric_limits<Real>::max();
  Real real_min = std::numeric_limits<Real>::min();
  AthenaArray<Real> vol(pmb->ncells1);
  const int nhistory_output = NHISTORY_VARS + pm->nuser_history_output_;
  std::unique_ptr<Real[]> hst_data(new Real[nhistory_output]);
  // initialize built-in variable sums to 0.0
  for (int n=0; n<NHISTORY_VARS; ++n) hst_data[n] = 0.0;
  // initialize user-defined history outputs depending on the requested operation
  for (int n=0; n<pm->nuser_history_output_; n++) {
    switch (pm->user_history_ops_[n]) {
      case UserHistoryOperation::sum:
        hst_data[NHISTORY_VARS+n] = 0.0;
        break;
      case UserHistoryOperation::max:
        hst_data[NHISTORY_VARS+n] = real_min;
        break;
      case UserHistoryOperation::min:
        hst_data[NHISTORY_VARS+n] = real_max;
        break;
    }
  }

  // Loop over MeshBlocks
  while (pmb != nullptr) {
    Hydro *phyd = pmb->phydro;
    Field *pfld = pmb->pfield;
    PassiveScalars *psclr = pmb->pscalars;
    Gravity *pgrav = pmb->pgrav;
    // BD: new problem
    Wave  *pwave = pmb->pwave;
    // -BD
    Advection  *padv = pmb->padv;
    Z4c *pz4c = pmb->pz4c;

    // Sum history variables over cells.  Note ghost cells are never included in sums
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          // NEW_OUTPUT_TYPES:

          int isum = 0;
          if (FLUID_ENABLED) {
            // Hydro conserved variables:
            Real& u_d  = phyd->u(IDN,k,j,i);
            Real& u_mx = phyd->u(IM1,k,j,i);
            Real& u_my = phyd->u(IM2,k,j,i);
            Real& u_mz = phyd->u(IM3,k,j,i);

            hst_data[isum++] += vol(i)*u_d;
            hst_data[isum++] += vol(i)*u_mx;
            hst_data[isum++] += vol(i)*u_my;
            hst_data[isum++] += vol(i)*u_mz;
            // + partitioned KE by coordinate direction:
            hst_data[isum++] += vol(i)*0.5*SQR(u_mx)/u_d;
            hst_data[isum++] += vol(i)*0.5*SQR(u_my)/u_d;
            hst_data[isum++] += vol(i)*0.5*SQR(u_mz)/u_d;

            if (NON_BAROTROPIC_EOS) {
              Real& u_e = phyd->u(IEN,k,j,i);;
              hst_data[isum++] += vol(i)*u_e;
            }
            // Graviatational potential energy:
            if (SELF_GRAVITY_ENABLED) {
              Real& phi = pgrav->phi(k,j,i);
              hst_data[isum++] += vol(i)*0.5*u_d*phi;
            }
            // Cell-centered magnetic energy, partitioned by coordinate direction:
            if (MAGNETIC_FIELDS_ENABLED) {
              Real& bcc1 = pfld->bcc(IB1,k,j,i);
              Real& bcc2 = pfld->bcc(IB2,k,j,i);
              Real& bcc3 = pfld->bcc(IB3,k,j,i);
              // constexpr int prev_out = NHYDRO + 3 + SELF_GRAVITY_ENABLED;
              hst_data[isum++] += vol(i)*0.5*bcc1*bcc1;
              hst_data[isum++] += vol(i)*0.5*bcc2*bcc2;
              hst_data[isum++] += vol(i)*0.5*bcc3*bcc3;
            }
            // (conserved variable) Passive scalars:
            for (int n=0; n<NSCALARS; n++) {
              Real& s = psclr->s(n,k,j,i);
              // constexpr int prev_out = NHYDRO + 3 + SELF_GRAVITY_ENABLED + NFIELD;
              hst_data[isum++] += vol(i)*s;
            }
          }

          // BD: new problem
          if (WAVE_ENABLED) {
            Real& wave_error = pwave->error(k,j,i);
            hst_data[isum++] += vol(i)*wave_error;
            hst_data[isum++] += vol(i)*SQR(wave_error);
          }
          // -BD

          if (ADVECTION_ENABLED) {
            Real& adv_error = padv->error(k,j,i);
            hst_data[isum++] += vol(i)*adv_error;
            hst_data[isum++] += vol(i)*SQR(adv_error);
          }

          if (Z4C_ENABLED) {
            Real const H_err  = std::abs(pz4c->con.H(k,j,i));
            Real const M2_err = std::abs(pz4c->con.M(k,j,i));
            Real const Mx_err = std::abs(pz4c->con.M_d(0,k,j,i));
            Real const My_err = std::abs(pz4c->con.M_d(1,k,j,i));
            Real const Mz_err = std::abs(pz4c->con.M_d(2,k,j,i));
            Real const Z2_err = std::abs(pz4c->con.Z(k,j,i));
            Real const theta  = std::abs(pz4c->z4c.Theta(k,j,i));
            Real const C2_err = std::abs(pz4c->con.C(k,j,i));

            hst_data[isum++] += vol(i)*SQR(H_err);
            hst_data[isum++] += vol(i)*M2_err; //M is already squared
            hst_data[isum++] += vol(i)*SQR(Mx_err);
            hst_data[isum++] += vol(i)*SQR(My_err);
            hst_data[isum++] += vol(i)*SQR(Mz_err);
            hst_data[isum++] += vol(i)*Z2_err; //Z is already squared
            hst_data[isum++] += vol(i)*SQR(theta);
            hst_data[isum++] += vol(i)*C2_err; //C is already squared
          }

        }
      }
    }
    for (int n=0; n<pm->nuser_history_output_; n++) { // user-defined history outputs
      if (pm->user_history_func_[n] != nullptr) {
        Real usr_val = pm->user_history_func_[n](pmb, n);
        switch (pm->user_history_ops_[n]) {
          case UserHistoryOperation::sum:
            // TODO(felker): this should automatically volume-weight the sum, like the
            // built-in variables. But existing user-defined .hst fns are currently
            // weighting their returned values.
            hst_data[NHISTORY_VARS+n] += usr_val;
            break;
          case UserHistoryOperation::max:
            hst_data[NHISTORY_VARS+n] = std::max(usr_val, hst_data[NHISTORY_VARS+n]);
            break;
          case UserHistoryOperation::min:
            hst_data[NHISTORY_VARS+n] = std::min(usr_val, hst_data[NHISTORY_VARS+n]);
            break;
        }
      }
    }
    pmb = pmb->next;
  }  // end loop over MeshBlocks

#ifdef MPI_PARALLEL
  // sum built-in/predefined hst_data[] over all ranks
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, hst_data.get(), NHISTORY_VARS, MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(hst_data.get(), hst_data.get(), NHISTORY_VARS, MPI_ATHENA_REAL, MPI_SUM,
               0, MPI_COMM_WORLD);
  }
  // apply separate chosen operations to each user-defined history output
  for (int n=0; n<pm->nuser_history_output_; n++) {
    Real *usr_hst_data = hst_data.get() + NHISTORY_VARS + n;
    MPI_Op usr_op;
    switch (pm->user_history_ops_[n]) {
      case UserHistoryOperation::sum:
        usr_op = MPI_SUM;
        break;
      case UserHistoryOperation::max:
        usr_op = MPI_MAX;
        break;
      case UserHistoryOperation::min:
        usr_op = MPI_MIN;
        break;
    }
    if (Globals::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, usr_hst_data, 1, MPI_ATHENA_REAL, usr_op, 0,
                 MPI_COMM_WORLD);
    } else {
      MPI_Reduce(usr_hst_data, usr_hst_data, 1, MPI_ATHENA_REAL, usr_op, 0,
                 MPI_COMM_WORLD);
    }
  }
#endif

  // only the master rank writes the file
  // create filename: "file_basename" + ".hst".  There is no file number.
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign(output_params.file_basename);
    fname.append(".hst");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
      msg << "### FATAL ERROR in function [OutputType::HistoryFile]" << std::endl
          << "Output file '" << fname << "' could not be opened";
      ATHENA_ERROR(msg);
    }

    // If this is the first output, write header
    if (output_params.file_number == 0) {
      // NEW_OUTPUT_TYPES:

      int iout = 1;
      std::fprintf(pfile,"# Athena++ history data\n"); // descriptor is first line
      std::fprintf(pfile,"# [%d]=time     ", iout++);
      std::fprintf(pfile,"[%d]=dt       ", iout++);
      if (FLUID_ENABLED) {
        std::fprintf(pfile,"[%d]=mass     ", iout++);
        std::fprintf(pfile,"[%d]=1-mom    ", iout++);
        std::fprintf(pfile,"[%d]=2-mom    ", iout++);
        std::fprintf(pfile,"[%d]=3-mom    ", iout++);
        std::fprintf(pfile,"[%d]=1-KE     ", iout++);
        std::fprintf(pfile,"[%d]=2-KE     ", iout++);
        std::fprintf(pfile,"[%d]=3-KE     ", iout++);
        if (NON_BAROTROPIC_EOS) std::fprintf(pfile,"[%d]=tot-E   ", iout++);
        if (SELF_GRAVITY_ENABLED) std::fprintf(pfile,"[%d]=grav-E   ", iout++);
        if (MAGNETIC_FIELDS_ENABLED) {
          std::fprintf(pfile,"[%d]=1-ME    ", iout++);
          std::fprintf(pfile,"[%d]=2-ME    ", iout++);
          std::fprintf(pfile,"[%d]=3-ME    ", iout++);
        }
        for (int n=0; n<NSCALARS; n++) {
          std::fprintf(pfile,"[%d]=%d-scalar    ", iout++, n);
        }
      }

      // BD: new problem
      if (WAVE_ENABLED) {
        std::fprintf(pfile,"[%d]=err-norm1 ", iout++);
        std::fprintf(pfile,"[%d]=err-norm2 ", iout++);
      }
      // -BD

      if (ADVECTION_ENABLED) {
        std::fprintf(pfile,"[%d]=err-norm1 ", iout++);
        std::fprintf(pfile,"[%d]=err-norm2 ", iout++);
      }

      if (Z4C_ENABLED) {
        std::fprintf(pfile,"[%d]=H-norm2 ",     iout++);
        std::fprintf(pfile,"[%d]=M-norm2 ",     iout++);
        std::fprintf(pfile,"[%d]=Mx-norm2 ",    iout++);
        std::fprintf(pfile,"[%d]=My-norm2 ",    iout++);
        std::fprintf(pfile,"[%d]=Mz-norm2 ",    iout++);
        std::fprintf(pfile,"[%d]=Z-norm2 ",     iout++);
        std::fprintf(pfile,"[%d]=Theta-norm2 ", iout++);
        std::fprintf(pfile,"[%d]=C-norm2 ",     iout++);
      }

      for (int n=0; n<pm->nuser_history_output_; n++)
        std::fprintf(pfile,"[%d]=%-8s", iout++,
                     pm->user_history_output_names_[n].c_str());
      std::fprintf(pfile,"\n");                              // terminate line
    }

    // write history variables
    std::fprintf(pfile, output_params.data_format.c_str(), pm->time);
    std::fprintf(pfile, output_params.data_format.c_str(), pm->dt);
    for (int n=0; n<nhistory_output; ++n)
      std::fprintf(pfile, output_params.data_format.c_str(), hst_data[n]);
    std::fprintf(pfile,"\n"); // terminate line
    std::fclose(pfile);
  }

  // increment counters, clean up
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  return;
}
