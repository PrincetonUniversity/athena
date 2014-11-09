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

// Primary header
#include "athena.hpp"

// C headers
#include <stdint.h>  // int64_t

// C++ headers
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <exception>  // exception
#include <iomanip>    // setprecision()
#include <iostream>   // cout, endl
#include <new>        // bad_alloc
#include <string>     // string

// Athena headers
#include "mesh.hpp"             // Mesh
#include "parameter_input.hpp"  // ParameterInput
#include "outputs/outputs.hpp"  // Outputs

// OpenMP headers
#include <omp.h>

// function prototypes
void ShowConfig();
void ChangeToRunDir(const char *pdir);

//======================================================================================
/////////////////////////////////// Athena++ Main Program \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//! \file main.cpp
//  \brief Athena++ main program file.
//
//  Based on the Athena MHD code in C.  History:
//    - 2003-2013: development of v1.0-v4.2 of Athena code in C
//    - Jan 2014: start of Athena++ code project
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn int main(int argc, char *argv[]) 
//  \brief Athena++ main program   */

int main(int argc, char *argv[])
{
  std::string athena_version = "version 0.1 - February 2014";
  std::string input_file = "athinput";
  char *prestart_file = NULL;
  char *prundir = NULL;
  int res_flag=0;     // gets set to 1 if -r        argument is on cmdline
  int narg_flag=0;    // gets set to 1 if -n        argument is on cmdline
  int iarg_flag=0;    // gets set to 1 if -i <file> argument is on cmdline

//--- Step 1. --------------------------------------------------------------------------
// Check for command line options and respond. 

  for (int i=1; i<argc; i++) {

// If argv[i] is a 2 character string of the form "-?" then:

    if(*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0'){
      switch(*(argv[i]+1)) {
      case 'i':                      // -i <input_file>
        input_file = argv[++i];
        iarg_flag = 1;
      break;
      case 'r':                      // -r <restart_file>
        res_flag = 1;
        prestart_file = argv[++i];
        if(iarg_flag) input_file = prestart_file; // use restart if input file not set 
        break;
      case 'd':                      // -d <run_directory>
        prundir = argv[++i];
        break;
      case 'n':
        narg_flag = 1;
        break;
      case 'c':
        ShowConfig();
        return(0);
      break;
      case 'h':
      default:
        std::cout<<"Athena++ "<< athena_version << std::endl;
        std::cout<<"Usage: "<< argv[0] <<" [options] [block/par=value ...]"<< std::endl;
        std::cout<<"Options:" << std::endl;
        std::cout<<"  -i <file>       specify input file [athinput]"<< std::endl;
        std::cout<<"  -r <file>       restart with this file"<< std::endl;
        std::cout<<"  -d <directory>  specify run dir [current dir]"<< std::endl;
        std::cout<<"  -n              parse input file and quit"<< std::endl;
        std::cout<<"  -c              show configuration and quit"<< std::endl;
        std::cout<<"  -h              this help"<< std::endl;
        ShowConfig();
        return(0);
        break;
      }
    } // else if argv[i] not of form "-?" ignore it
  }

//--- Step 2. --------------------------------------------------------------------------
// Construct object to store input parameters, then parse input file and command line
// Note memory allocations and parameter input are protected by a simple error handler

  ParameterInput *pinput;
  try {
    pinput = new ParameterInput;
    pinput->LoadFromFile(input_file);
    pinput->ModifyFromCmdline(argc,argv);
  } 
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class ParameterInput: " 
              << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message  
    return(0);
  }

// Dump input parameters and quit if code was run with -n option.

  if (narg_flag){
    pinput->ParameterDump(std::cout);
    return(0);
  }

//--- Step 3. --------------------------------------------------------------------------
// Initialize MPI environment, distribute input parameters to all ranks

//  g_comm_world_id = 0;

// Note steps 4-6 are protected by a simple error handler
//--- Step 4. --------------------------------------------------------------------------
// Construct and initialize Mesh

  Mesh *pmesh;
  try {
    pmesh = new Mesh(pinput);
  }
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed initializing class Mesh: " 
              << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message
    return(0);
  }

//--- Step 5. --------------------------------------------------------------------------
// Set initial conditions by calling problem generator on each MeshDomain/MeshBlock

  try {
    pmesh->ForAllDomains(pgen,pinput);
  } 
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl << "memory allocation failed "
              << "in problem generator " << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message 
    return(0);
  }

// apply BCs, compute primitive from conserved variables, compute first timestep

  pmesh->ForAllDomains( fluid_bcs_n,pinput);
  if (MAGNETIC_FIELDS_ENABLED) pmesh->ForAllDomains(bfield_bcs_n,pinput);
  pmesh->ForAllDomains(primitives_n,pinput);
  pmesh->ForAllDomains(new_timestep,pinput);

//--- Step 6. --------------------------------------------------------------------------
// Change to run directory, initialize outputs object, and make output of ICs

  Outputs *pouts;
  try {
    ChangeToRunDir(prundir);
    pouts = new Outputs(pmesh, pinput);
    pouts->MakeOutputs(pmesh);
  } 
  catch(std::bad_alloc& ba) {
    std::cout << "### FATAL ERROR in main" << std::endl
              << "memory allocation failed setting initial conditions: " 
              << ba.what() << std::endl;
    return(0);
  }
  catch(std::exception const& ex) {
    std::cout << ex.what() << std::endl;  // prints diagnostic message  
    return(0);
  }

//--- Step 9. === START OF MAIN INTEGRATION LOOP =======================================
// For performance, there is no error handler protecting this step (except outputs)

  std::cout<<std::endl<< "Setup complete, entering main loop..." <<std::endl<<std::endl;
  clock_t tstart = clock();
#ifdef OPENMP_PARALLEL
  double omp_time = omp_get_wtime();
#endif

  while ((pmesh->time < pmesh->tlim) && 
         (pmesh->nlim < 0 || pmesh->ncycle < pmesh->nlim)){
    std::cout << "cycle=" << pmesh->ncycle << std::scientific << std::setprecision(5)
              << " time=" << pmesh->time << " dt=" << pmesh->dt << std::endl;

// predict step

    pmesh->ForAllDomains( fluid_predict  ,pinput);
    pmesh->ForAllDomains( fluid_bcs_nhalf,pinput);

    if (MAGNETIC_FIELDS_ENABLED) {
      pmesh->ForAllDomains(bfield_predict  ,pinput);
      pmesh->ForAllDomains(bfield_bcs_nhalf,pinput);
    }

    pmesh->ForAllDomains(primitives_nhalf,pinput);

// correct step

    pmesh->ForAllDomains( fluid_correct,pinput);
    pmesh->ForAllDomains( fluid_bcs_n,  pinput);

    if (MAGNETIC_FIELDS_ENABLED) {
      pmesh->ForAllDomains(bfield_correct,pinput);
      pmesh->ForAllDomains(bfield_bcs_n,  pinput);
    }

    pmesh->ForAllDomains(primitives_n,  pinput);

// new time step, outputs, diagnostics

    pmesh->ncycle++;
    pmesh->time  += pmesh->dt;

    try {
      pouts->MakeOutputs(pmesh);
    } 
    catch(std::bad_alloc& ba) {
      std::cout << "### FATAL ERROR in main" << std::endl
                << "memory allocation failed during output: " << ba.what() << std::endl;
      return(0);
    }
    catch(std::exception const& ex) {
      std::cout << ex.what() << std::endl;  // prints diagnostic message  
      return(0);
    }

    pmesh->ForAllDomains(new_timestep,pinput);

  } // END OF MAIN INTEGRATION LOOP ====================================================
#ifdef OPENMP_PARALLEL
  omp_time = omp_get_wtime() - omp_time;;
#endif
  clock_t tstop = clock();

// print diagnostic messages

  std::cout << "cycle=" << pmesh->ncycle << std::scientific << std::setprecision(5)
            << " time=" << pmesh->time << " dt=" << pmesh->dt << std::endl;

  if (pmesh->ncycle == pmesh->nlim) {
    std::cout << std::endl << "Terminating on cycle limit" << std::endl;
  } else {
    std::cout << std::endl << "Terminating on time limit" << std::endl;
  }

  std::cout << "time=" << pmesh->time << " cycle=" << pmesh->ncycle << std::endl;
  std::cout << "tlim=" << pmesh->tlim << " nlim=" << pmesh->nlim << std::endl;

// Calculate and print the zone-cycles/cpu-second and wall-second

  float cpu_time = (tstop>tstart ? (float)(tstop-tstart) : 1.0)/(float)CLOCKS_PER_SEC;
  int64_t zones = (pmesh->mesh_size.nx1)*(pmesh->mesh_size.nx2)*(pmesh->mesh_size.nx3);
  float zc_cpus = (float)(zones*pmesh->ncycle)/cpu_time;
  std::cout << std::endl << "cpu time used  = " << cpu_time << std::endl;
  std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
  float zc_omps = (float)(zones*pmesh->ncycle)/omp_time;
  std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
  std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif

  delete pinput;
  delete pmesh;
  delete pouts;

  return(0); 
}
