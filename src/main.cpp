//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//================================= Athena++ Main Program ================================
//! \file main.cpp
//! \brief Athena++ main program
//!
//! Based on the Athena MHD code (Cambridge version), originally written in 2002-2005 by
//! Jim Stone, Tom Gardiner, and Peter Teuben, with many important contributions by many
//! other developers after that, i.e. 2005-2014.
//!
//! Athena++ was started in Jan 2014.  The core design was finished during 4-7/2014 at the
//! KITP by Jim Stone.  GR was implemented by Chris White and AMR by Kengo Tomida during
//! 2014-2016.  Contributions from many others have continued to the present.
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <csignal>    // ISO C/C++ signal() and sigset_t, sigemptyset() POSIX C extensions
#include <cstdint>    // int64_t
#include <cstdio>     // sscanf()
#include <cstdlib>    // strtol
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <exception>  // exception
#include <iomanip>    // setprecision()
#include <iostream>   // cout, endl
#include <limits>     // max_digits10
#include <new>        // bad_alloc
#include <string>     // string

// Athena++ headers
#include "athena.hpp"
#include "fft/turbulence.hpp"
#include "globals.hpp"
#include "gravity/fft_gravity.hpp"
#include "gravity/mg_gravity.hpp"
#include "mesh/mesh.hpp"
#include "outputs/io_wrapper.hpp"
#include "outputs/outputs.hpp"
#include "parameter_input.hpp"
#include "utils/utils.hpp"

// MPI/OpenMP headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn int main(int argc, char *argv[])
//! \brief Athena++ main program

int main(int argc, char *argv[]) {
  std::string athena_version = "version 21.0 - January 2021";
  char *input_filename = nullptr, *restart_filename = nullptr;
  char *prundir = nullptr;
  int res_flag = 0;   // set to 1 if -r        argument is on cmdline
  int narg_flag = 0;  // set to 1 if -n        argument is on cmdline
  int iarg_flag = 0;  // set to 1 if -i <file> argument is on cmdline
  int mesh_flag = 0;  // set to <nproc> if -m <nproc> argument is on cmdline
  int wtlim = 0;
  std::uint64_t mbcnt = 0;

  //--- Step 1. --------------------------------------------------------------------------
  // Initialize MPI environment, if necessary

#ifdef MPI_PARALLEL
#ifdef OPENMP_PARALLEL
  int mpiprv;
  if (MPI_SUCCESS != MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpiprv)) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "MPI Initialization failed." << std::endl;
    return(0);
  }
  if (mpiprv != MPI_THREAD_MULTIPLE) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "MPI_THREAD_MULTIPLE must be supported for the hybrid parallelzation. "
              << MPI_THREAD_MULTIPLE << " : " << mpiprv
              << std::endl;
    MPI_Finalize();
    return(0);
  }
#else  // no OpenMP
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "MPI Initialization failed." << std::endl;
    return(0);
  }
#endif  // OPENMP_PARALLEL
  // Get process id (rank) in MPI_COMM_WORLD
  if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &(Globals::my_rank))) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "MPI_Comm_rank failed." << std::endl;
    MPI_Finalize();
    return(0);
  }

  // Get total number of MPI processes (ranks)
  if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &Globals::nranks)) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "MPI_Comm_size failed." << std::endl;
    MPI_Finalize();
    return(0);
  }
#else  // no MPI
  Globals::my_rank = 0;
  Globals::nranks  = 1;
#endif  // MPI_PARALLEL

  //--- Step 2. --------------------------------------------------------------------------
  // Check for command line options and respond.

  for (int i=1; i<argc; i++) {
    // If argv[i] is a 2 character string of the form "-?" then:
    if (*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0') {
      // check validity of command line options + arguments:
      char opt_letter = *(argv[i]+1);
      switch(opt_letter) {
        // options that do not take arguments:
        case 'n':
        case 'c':
        case 'h':
          break;
          // options that require arguments:
        default:
          if ((i+1 >= argc) // flag is at the end of the command line options
              || (*argv[i+1] == '-') ) { // flag is followed by another flag
            if (Globals::my_rank == 0) {
              std::cout << "### FATAL ERROR in main" << std::endl
                        << "-" << opt_letter << " must be followed by a valid argument\n";
#ifdef MPI_PARALLEL
              MPI_Finalize();
#endif
              return(0);
            }
          }
      }
      switch(*(argv[i]+1)) {
        case 'i':                      // -i <input_filename>
          input_filename = argv[++i];
          iarg_flag = 1;
          break;
        case 'r':                      // -r <restart_file>
          res_flag = 1;
          restart_filename = argv[++i];
          break;
        case 'd':                      // -d <run_directory>
          prundir = argv[++i];
          break;
        case 'n':
          narg_flag = 1;
          break;
        case 'm':                      // -m <nproc>
          mesh_flag = static_cast<int>(std::strtol(argv[++i], nullptr, 10));
          break;
        case 't':                      // -t <hh:mm:ss>
          int wth, wtm, wts;
          std::sscanf(argv[++i], "%d:%d:%d", &wth, &wtm, &wts);
          wtlim = wth*3600 + wtm*60 + wts;
          break;
        case 'c':
          if (Globals::my_rank == 0) ShowConfig();
#ifdef MPI_PARALLEL
          MPI_Finalize();
#endif
          return(0);
          break;
        case 'h':
        default:
          if (Globals::my_rank == 0) {
            std::cout << "Athena++ " << athena_version << std::endl;
            std::cout << "Usage: " << argv[0] << " [options] [block/par=value ...]\n";
            std::cout << "Options:" << std::endl;
            std::cout << "  -i <file>       specify input file [athinput]\n";
            std::cout << "  -r <file>       restart with this file\n";
            std::cout << "  -d <directory>  specify run dir [current dir]\n";
            std::cout << "  -n              parse input file and quit\n";
            std::cout << "  -c              show configuration and quit\n";
            std::cout << "  -m <nproc>      output mesh structure and quit\n";
            std::cout << "  -t hh:mm:ss     wall time limit for final output\n";
            std::cout << "  -h              this help\n";
            ShowConfig();
          }
#ifdef MPI_PARALLEL
          MPI_Finalize();
#endif
          return(0);
          break;
      }
    } // else if argv[i] not of form "-?" ignore it here (tested in ModifyFromCmdline)
  }

  if (restart_filename == nullptr && input_filename == nullptr) {
    // no input file is given
    std::cout << "### FATAL ERROR in main" << std::endl
              << "No input file or restart file is specified." << std::endl;
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }

  // Set up the signal handler
  SignalHandler::SignalHandlerInit();
  if (Globals::my_rank == 0 && wtlim > 0)
    SignalHandler::SetWallTimeAlarm(wtlim);

  // Note steps 3-6 are protected by a simple error handler
  //--- Step 3. --------------------------------------------------------------------------
  // Construct object to store input parameters, then parse input file and command line.
  // With MPI, the input is read by every process in parallel using MPI-IO.

  ParameterInput *pinput;
  IOWrapper infile, restartfile;
#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    pinput = new ParameterInput;
    if (res_flag == 1) {
      restartfile.Open(restart_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(restartfile);
      // make sure next_time gets corrected in case -i input file or cmdline args change
      // the output next_time, dt, etc.
      // This needs to be corrected on the restart file because we need the old dt.
      pinput->RollbackNextTime();
      // leave the restart file open for later use
    }
    if (iarg_flag == 1) {
      // if both -r and -i are specified, override the parameters using the input file
      infile.Open(input_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(infile);
      infile.Close();
    }
    pinput->ModifyFromCmdline(argc ,argv);
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class ParameterInput: "
              << ba.what() << std::endl;
    if (res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    if (res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  //--- Step 4. --------------------------------------------------------------------------
  // Construct and initialize Mesh

  Mesh *pmesh;
#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    if (res_flag == 0) {
      pmesh = new Mesh(pinput, mesh_flag);
    } else {
      pmesh = new Mesh(pinput, restartfile, mesh_flag);
    }
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class Mesh: "
              << ba.what() << std::endl;
    if (res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    if (res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  // With current mesh time possibly read from restart file, correct next_time for outputs
  if (res_flag == 1) {
    // ensure that next_time  >= mesh_time - dt, in case input file or command line
    // overrides it
    pinput->ForwardNextTime(pmesh->time);
  }

  // Dump input parameters and quit if code was run with -n option.
  if (narg_flag) {
    if (Globals::my_rank == 0) pinput->ParameterDump(std::cout);
    if (res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }

  if (res_flag == 1) restartfile.Close(); // close the restart file here

  // Quit if -m was on cmdline.  This option builds and outputs mesh structure.
  if (mesh_flag > 0) {
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }

  //--- Step 5. --------------------------------------------------------------------------
  // Construct and initialize TaskList

  TimeIntegratorTaskList *ptlist;
#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    ptlist = new TimeIntegratorTaskList(pinput, pmesh);
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl << "memory allocation failed "
              << "in creating task list " << ba.what() << std::endl;
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  SuperTimeStepTaskList *pststlist = nullptr;
  if (STS_ENABLED) {
#ifdef ENABLE_EXCEPTIONS
    try {
#endif
      pststlist = new SuperTimeStepTaskList(pinput, pmesh, ptlist);
#ifdef ENABLE_EXCEPTIONS
    }
    catch(std::bad_alloc& ba) {
      std::cout << "### FATAL ERROR in main" << std::endl << "memory allocation failed "
                << "in creating task list " << ba.what() << std::endl;
#ifdef MPI_PARALLEL
      MPI_Finalize();
#endif
      return(0);
    }
#endif // ENABLE_EXCEPTIONS
  }

  //--- Step 6. --------------------------------------------------------------------------
  // Set initial conditions by calling problem generator, or reading restart file

#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    pmesh->Initialize(res_flag, pinput);
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl << "memory allocation failed "
              << "in problem generator " << ba.what() << std::endl;
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  //--- Step 7. --------------------------------------------------------------------------
  // Change to run directory, initialize outputs object, and make output of ICs

  Outputs *pouts;
#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    ChangeRunDir(prundir);
    pouts = new Outputs(pmesh, pinput);
    if (res_flag == 0) pouts->MakeOutputs(pmesh, pinput);
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed setting initial conditions: "
              << ba.what() << std::endl;
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  //=== Step 8. === START OF MAIN INTEGRATION LOOP =======================================
  // For performance, there is no error handler protecting this step (except outputs)

  if (Globals::my_rank == 0) {
    std::cout << "\nSetup complete, entering main loop...\n" << std::endl;
  }

  clock_t tstart = clock();
#ifdef OPENMP_PARALLEL
  double omp_start_time = omp_get_wtime();
#endif

  while ((pmesh->time < pmesh->tlim) &&
         (pmesh->nlim < 0 || pmesh->ncycle < pmesh->nlim)) {
    if (Globals::my_rank == 0)
      pmesh->OutputCycleDiagnostics();

    if (STS_ENABLED) {
      pmesh->sts_loc = TaskType::op_split_before;
      // compute nstages for this STS
      if (pmesh->sts_integrator == "rkl2") { // default
        pststlist->nstages =
            static_cast<int>
              (0.5*(-1. + std::sqrt(9. + 16.*(0.5*pmesh->dt)/pmesh->dt_parabolic))) + 1;
      } else { // rkl1
        pststlist->nstages =
            static_cast<int>
              (0.5*(-1. + std::sqrt(1. + 8.*pmesh->dt/pmesh->dt_parabolic))) + 1;
      }
      if (pststlist->nstages % 2 == 0) { // guarantee odd nstages for STS
        pststlist->nstages += 1;
      }
      // take super-timestep
      for (int stage=1; stage<=pststlist->nstages; ++stage)
        pststlist->DoTaskListOneStage(pmesh, stage);

      pmesh->sts_loc = TaskType::main_int;
    }

    if (pmesh->turb_flag > 1) pmesh->ptrbd->Driving(); // driven turbulence

    for (int stage=1; stage<=ptlist->nstages; ++stage) {
      ptlist->DoTaskListOneStage(pmesh, stage);
      if (ptlist->CheckNextMainStage(stage)) {
        if (SELF_GRAVITY_ENABLED == 1) // fft (0: discrete kernel, 1: continuous kernel)
          pmesh->pfgrd->Solve(stage, 0);
        else if (SELF_GRAVITY_ENABLED == 2) // multigrid
          pmesh->pmgrd->Solve(stage);
      }
    }

    if (STS_ENABLED && pmesh->sts_integrator == "rkl2") {
      pmesh->sts_loc = TaskType::op_split_after;
      // take super-timestep
      for (int stage=1; stage<=pststlist->nstages; ++stage)
        pststlist->DoTaskListOneStage(pmesh, stage);
    }

    pmesh->UserWorkInLoop();

    pmesh->ncycle++;
    pmesh->time += pmesh->dt;
    mbcnt += pmesh->nbtotal;
    pmesh->step_since_lb++;

    pmesh->LoadBalancingAndAdaptiveMeshRefinement(pinput);

    pmesh->NewTimeStep();

#ifdef ENABLE_EXCEPTIONS
    try {
#endif
      if (pmesh->time < pmesh->tlim) // skip the final output as it happens later
        pouts->MakeOutputs(pmesh,pinput);
#ifdef ENABLE_EXCEPTIONS
    }
    catch(std::bad_alloc& ba) {
      std::cout << "### FATAL ERROR in main" << std::endl
                << "memory allocation failed during output: " << ba.what() <<std::endl;
#ifdef MPI_PARALLEL
      MPI_Finalize();
#endif
      return(0);
    }
    catch(std::exception const& ex) {
      std::cout << ex.what() << std::endl;  // prints diagnostic message
#ifdef MPI_PARALLEL
      MPI_Finalize();
#endif
      return(0);
    }
#endif // ENABLE_EXCEPTIONS

    // check for signals
    if (SignalHandler::CheckSignalFlags() != 0) {
      break;
    }
  } // END OF MAIN INTEGRATION LOOP ======================================================
  // Make final outputs, print diagnostics, clean up and terminate

  if (Globals::my_rank == 0 && wtlim > 0)
    SignalHandler::CancelWallTimeAlarm();


  //--- Step 9. --------------------------------------------------------------------------
  // Output the final cycle diagnostics and make the final outputs

  if (Globals::my_rank == 0)
    pmesh->OutputCycleDiagnostics();

  pmesh->UserWorkAfterLoop(pinput);

#ifdef ENABLE_EXCEPTIONS
  try {
#endif
    pouts->MakeOutputs(pmesh,pinput,true);
#ifdef ENABLE_EXCEPTIONS
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed during output: " << ba.what() <<std::endl;
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    return(0);
  }
#endif // ENABLE_EXCEPTIONS

  //--- Step 10. -------------------------------------------------------------------------
  // Print diagnostic messages related to the end of the simulation

  if (Globals::my_rank == 0) {
    if (SignalHandler::GetSignalFlag(SIGTERM) != 0) {
      std::cout << std::endl << "Terminating on Terminate signal" << std::endl;
    } else if (SignalHandler::GetSignalFlag(SIGINT) != 0) {
      std::cout << std::endl << "Terminating on Interrupt signal" << std::endl;
    } else if (SignalHandler::GetSignalFlag(SIGALRM) != 0) {
      std::cout << std::endl << "Terminating on wall-time limit" << std::endl;
    } else if (pmesh->ncycle == pmesh->nlim) {
      std::cout << std::endl << "Terminating on cycle limit" << std::endl;
    } else {
      std::cout << std::endl << "Terminating on time limit" << std::endl;
    }

    std::cout << "time=" << pmesh->time << " cycle=" << pmesh->ncycle << std::endl;
    std::cout << "tlim=" << pmesh->tlim << " nlim=" << pmesh->nlim << std::endl;

    if (pmesh->adaptive) {
      std::cout << std::endl << "Number of MeshBlocks = " << pmesh->nbtotal
                << "; " << pmesh->nbnew << "  created, " << pmesh->nbdel
                << " destroyed during this simulation." << std::endl;
    }

    // Calculate and print the zone-cycles/cpu-second and wall-second
#ifdef OPENMP_PARALLEL
    double omp_time = omp_get_wtime() - omp_start_time;
#endif
    clock_t tstop = clock();
    double cpu_time = (tstop>tstart ? static_cast<double> (tstop-tstart) :
                       1.0)/static_cast<double> (CLOCKS_PER_SEC);
    std::uint64_t zonecycles = mbcnt
      *static_cast<std::uint64_t> (pmesh->my_blocks(0)->GetNumberOfMeshBlockCells());
    double zc_cpus = static_cast<double> (zonecycles) / cpu_time;

    std::cout << std::endl << "zone-cycles = " << zonecycles << std::endl;
    std::cout << "cpu time used  = " << cpu_time << std::endl;
    std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
    double zc_omps = static_cast<double> (zonecycles) / omp_time;
    std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
    std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
  }

  delete pinput;
  delete pmesh;
  delete ptlist;
  delete pouts;

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif

  return(0);
}
