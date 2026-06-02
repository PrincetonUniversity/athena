//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2026 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file int2d.cpp
//! \brief implements non-trivial functions of IntX[123]X[123]Output classes.

// C/C++ headers
#include <fstream>   // ifstream, ios, ofstream
#include <iostream>  // cout, endl, <<
#include <sstream>   // stringstream
#include <string>    // string, to_string
#include <vector>

// Athena++ headers
#include "../mesh/mesh.hpp"  // time, my_blocks
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
//! \fn void Int2DOutput::ProcessHeader(const std::string& ext, const Mesh *pm)
//! \brief writes a new or check the existing header of the output file.

void Int2DOutput::ProcessHeader(const std::string& ext, const Mesh *pm) {
  // Compose the name of the output file.
  fname = output_params.file_basename + '.' + ext;

  // Collect the output variables.
  LoadOutputData(pm->my_blocks(0));
  std::stringstream msg;
  std::vector<std::string> varnames;
  OutputData *pdata = pfirst_data_;
  while (pdata != nullptr) {
    switch (pdata->type[0]) {
      case 'S':  // SCALARS
        varnames.push_back(pdata->name);
        break;
      case 'V':  // VECTORS
        for (int i = 1; i <= 3; ++i)
          varnames.push_back(pdata->name + std::to_string(i));
        break;
      default:
        msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
            << "Unknown output variable type: " << pdata->type << std::endl;
        ATHENA_ERROR(msg);
        return;
    }
    pdata = pdata->pnext;
  }
  ClearOutputData();

  // Construct the coordinates of the third dimension.
  SetThirdDim(pm);

  // Process the header of the output file.
  const int nvar = varnames.size();
  std::ifstream fin(fname, std::ios::in | std::ios::binary);
  if (fin.is_open()) {
    // Read and check the number of output variables.
    int n;
    fin.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (n != nvar) {
      msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
          << "Inconsistent number of output variables: "
          << nvar << " requested vs. " << n << " in existing header" << std::endl;
      ATHENA_ERROR(msg);
      return;
    }

    // Read and check the names of the variables.
    for (int i = 0; i < nvar; ++i) {
      fin.read(reinterpret_cast<char*>(&n), sizeof(n));
      std::string name(n, '*');
      fin.read(&name[0], n);
      if (name != varnames[i]) {
        msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
            << "Inconsistent name of output variable " << (i + 1) << std::endl
            << "'" << varnames[i] << "' requested vs. '" << name << "' in existing header"
            << std::endl;
        ATHENA_ERROR(msg);
        return;
      }
    }

    fin.close();
  } else { // if (fin.is_open())
    // Open a new output file for write.
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    if (!fout.is_open()) {
      msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
          << "Unable to create output file '" << fname << "'" << std::endl;
      ATHENA_ERROR(msg);
      return;
    }

    // Write the number of output variables.
    fout.write(reinterpret_cast<const char*>(&nvar), sizeof(nvar));

    // Write the names of the variables.
    for (std::vector<std::string>::iterator it = varnames.begin();
        it != varnames.end(); ++it) {
      const int size = it->size();
      fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
      if (size > 0) fout.write(it->data(), size);
    }

    fout.close();
  } // if (fin.is_open())
}

//----------------------------------------------------------------------------------------
//! \fn void IntX1X2Output::SetThirdDim(const Mesh *pm)
//! \brief constructs the coordinates in the X3 dimension.

void IntX1X2Output::SetThirdDim(const Mesh *pm) {
}

//----------------------------------------------------------------------------------------
//! \fn void IntX1X2Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//! \brief integrates the data over x1 and x2 directions and writes the resulting 1D
//      array in x3.

void IntX1X2Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  std::cout << "IntX1X2Output: under construction; "
            << "\n\tt = " << pm->time
            << "\n\tnext_time = " << output_params.next_time
            << "\n\tdt = " << output_params.dt << std::endl;

  // Update output parameters.
  output_params.next_time += output_params.dt;
}

//----------------------------------------------------------------------------------------
//! \fn void IntX1X3Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//! \brief integrates the data over x1 and x3 directions and writes the resulting 1D
//      array in x2.

void IntX1X3Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  std::stringstream msg;
  msg << "IntX1X3Output: not implemented " << std::endl;
  ATHENA_ERROR(msg);
}

//----------------------------------------------------------------------------------------
//! \fn void IntX2X3Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//! \brief integrates the data over x2 and x3 directions and writes the resulting 1D
//      array in x1.

void IntX2X3Output::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  std::stringstream msg;
  msg << "IntX2X3Output: not implemented " << std::endl;
  ATHENA_ERROR(msg);
}
