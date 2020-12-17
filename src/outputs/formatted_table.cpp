//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file formatted_table.cpp
//! \brief writes output data as a formatted table.  Should not be used to output large
//! 3D data sets as this format is very slow and memory intensive.  Most useful for 1D
//! slices and/or sums.  Writes one file per Meshblock.

// C headers

// C++ headers
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "outputs.hpp"


//----------------------------------------------------------------------------------------
//! \fn void FormattedTableOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin,
//!                                                bool flag)
//! \brief writes OutputData to file in tabular format using C style std::fprintf
//!        Writes one file per MeshBlock

void FormattedTableOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  // Loop over MeshBlocks
  for (int b=0; b<pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    // set start/end array indices depending on whether ghost zones are included
    out_is = pmb->is; out_ie = pmb->ie;
    out_js = pmb->js; out_je = pmb->je;
    out_ks = pmb->ks; out_ke = pmb->ke;
    if (output_params.include_ghost_zones) {
      out_is -= NGHOST; out_ie += NGHOST;
      if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
      if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}
    }

    // build doubly linked list of OutputData nodes (setting data ptrs to appropriate
    // quantity on MeshBlock for each node), then slice/sum as needed
    LoadOutputData(pmb);
    if (!TransformOutputData(pmb)) {
      ClearOutputData();  // required when LoadOutputData() is used.
      continue;
    } // skip if slice was out of range

    // create filename: "file_basename"+ "."+"blockid"+"."+"file_id"+"."+XXXXX+".tab",
    // where XXXXX = 5-digit file_number
    std::string fname;
    char number[6];
    std::snprintf(number, sizeof(number), "%05d", output_params.file_number);
    char blockid[12];
    std::snprintf(blockid, sizeof(blockid), "block%d", pmb->gid);

    fname.assign(output_params.file_basename);
    fname.append(".");
    fname.append(blockid);
    fname.append(".");
    fname.append(output_params.file_id);
    fname.append(".");
    fname.append(number);
    fname.append(".tab");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
      msg << "### FATAL ERROR in function [FormattedTableOutput::WriteOutputFile]"
          << std::endl << "Output file '" << fname << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    // print file header
    std::fprintf(pfile, "# Athena++ data at time=%e", pm->time);
    std::fprintf(pfile, "  cycle=%d", pmb->pmy_mesh->ncycle);
    std::fprintf(pfile, "  variables=%s \n", output_params.variable.c_str());

    // write x1, x2, x3 column headers
    std::fprintf(pfile, "#");
    if (out_is != out_ie) std::fprintf(pfile, " i       x1v     ");
    if (out_js != out_je) std::fprintf(pfile, " j       x2v     ");
    if (out_ks != out_ke) std::fprintf(pfile, " k       x3v     ");
    // write data col headers from "name" stored in OutputData nodes of doubly linked list
    OutputData *pdata = pfirst_data_;
    while (pdata != nullptr) {
      if (pdata->type == "VECTORS") {
        for (int index = 1; index <= 3; ++index) {
          std::fprintf(pfile, "    %s%d     ", pdata->name.c_str(), index);
        }
      } else {
        std::fprintf(pfile, "    %s      ", pdata->name.c_str());
      }
      pdata = pdata->pnext;
    }
    std::fprintf(pfile, "\n"); // terminate line

    // loop over all cells in data arrays
    for (int k=out_ks; k<=out_ke; ++k) {
      for (int j=out_js; j<=out_je; ++j) {
        for (int i=out_is; i<=out_ie; ++i) {
          // write x1, x2, x3 indices and coordinates on start of new line
          if (out_is != out_ie) {
            std::fprintf(pfile, "%04d", i);
            std::fprintf(pfile, output_params.data_format.c_str(), pmb->pcoord->x1v(i));
          }
          if (out_js != out_je) {
            std::fprintf(pfile, " %04d", j);  // note extra space for formatting
            std::fprintf(pfile, output_params.data_format.c_str(), pmb->pcoord->x2v(j));
          }
          if (out_ks != out_ke) {
            std::fprintf(pfile, " %04d", k);  // note extra space for formatting
            std::fprintf(pfile, output_params.data_format.c_str(), pmb->pcoord->x3v(k));
          }

          // step through doubly linked list of OutputData's and write each on same line
          OutputData *pdata_inner_loop = pfirst_data_;
          while (pdata_inner_loop != nullptr) {
            for (int n=0; n<(pdata_inner_loop->data.GetDim4()); ++n) {
              std::fprintf(pfile, output_params.data_format.c_str(),
                           pdata_inner_loop->data(n,k,j,i));
            }
            pdata_inner_loop = pdata_inner_loop->pnext;
          }
          std::fprintf(pfile,"\n"); // terminate line
        }
      }
    }
    // don't forget to close the output file and clean up ptrs to data in OutputData
    std::fclose(pfile);
    ClearOutputData(); // required when LoadOutputData() is used.
  }  // end loop over MeshBlocks

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}
