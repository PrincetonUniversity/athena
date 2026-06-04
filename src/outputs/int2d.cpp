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
#include "../coordinates/coordinates.hpp"  // VolCenterFace3Area
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

  // Confirm the count of output variables with the parent class.
  const int nvar = varnames.size();
  if (nvar != num_vars_) {
    msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
        << "Inconsistent num_vars_ = " << num_vars_ << " vs. nvar = " << nvar
        << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Construct the coordinates of the third dimension.
  SetThirdDim(pm);

  // Process the header of the output file.
  const int rsize = sizeof(Real);
  std::ifstream fin(fname, std::ios::in | std::ios::binary);
  if (fin.is_open()) {
    // Stop if the run is a fresh start.
    // TODO(ccyang): implement truncation of the output file for restart.
    if (pm->time == pm->start_time) {
      msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
          << "The output file '" << fname << "' already exists. " << std::endl;
      ATHENA_ERROR(msg);
      return;
    }

    // Check the Real size.
    int n;
    fin.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (n != rsize) {
      msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
          << "Inconsistent Athena++ Real size: "
          << rsize << " bytes compied vs. " << n << " bytes in existing header"
          << std::endl;
      ATHENA_ERROR(msg);
      return;
    }

    // Read and check the number of output variables.
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

    // Read and check the number of cells in the third dimension.
    fin.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (n != nx) {
      msg << "### FATAL ERROR in Int2DOutput::ProcessHeader" << std::endl
          << "Inconsistent number of cells in the third dimension: "
          << nx << " requested vs. " << n << " in existing header" << std::endl;
      ATHENA_ERROR(msg);
      return;
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

    // Write the size of the Athena++ Real.
    fout.write(reinterpret_cast<const char*>(&rsize), sizeof(rsize));

    // Write the number of output variables.
    fout.write(reinterpret_cast<const char*>(&nvar), sizeof(nvar));

    // Write the names of the variables.
    for (std::vector<std::string>::iterator it = varnames.begin();
        it != varnames.end(); ++it) {
      const int size = it->size();
      fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
      if (size > 0) fout.write(it->data(), size);
    }

    // Write the coordinates in the third dimension.
    fout.write(reinterpret_cast<const char*>(&nx), sizeof(nx));
    for (std::vector<Real>::iterator it = xf.begin(); it != xf.end(); ++it) {
      const Real x = *it;
      fout.write(reinterpret_cast<const char*>(&x), sizeof(x));
    }

    fout.close();
  } // if (fin.is_open())
}

//----------------------------------------------------------------------------------------
//! \fn void IntX1X2Output::SetThirdDim(const Mesh *pm)
//! \brief constructs the coordinates in the X3 dimension.
// TODO(ccyang): consider the case of mesh refinement.

void IntX1X2Output::SetThirdDim(const Mesh *pm) {
  const bool uniform = pm->use_uniform_meshgen_fn_[X3DIR];
  xf.clear();
  nx = pm->mesh_size.nx3;
  for (int i = 0; i <= nx; ++i) {
    const Real rx = ComputeMeshGeneratorX(i, nx, uniform);
    xf.push_back(pm->MeshGenerator_[X3DIR](rx, pm->mesh_size));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void IntX1X3Output::SetThirdDim(const Mesh *pm)
//! \brief constructs the coordinates in the X2 dimension.
// TODO(ccyang): consider the case of mesh refinement.

void IntX1X3Output::SetThirdDim(const Mesh *pm) {
  std::stringstream msg;
  msg << "IntX1X3Output: not implemented " << std::endl;
  ATHENA_ERROR(msg);
  return;

  const bool uniform = pm->use_uniform_meshgen_fn_[X2DIR];
  xf.clear();
  nx = pm->mesh_size.nx2;
  for (int i = 0; i <= nx; ++i) {
    const Real rx = ComputeMeshGeneratorX(i, nx, uniform);
    xf.push_back(pm->MeshGenerator_[X2DIR](rx, pm->mesh_size));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void IntX2X3Output::SetThirdDim(const Mesh *pm)
//! \brief constructs the coordinates in the X1 dimension.
// TODO(ccyang): consider the case of mesh refinement.

void IntX2X3Output::SetThirdDim(const Mesh *pm) {
  std::stringstream msg;
  msg << "IntX2X3Output: not implemented " << std::endl;
  ATHENA_ERROR(msg);
  return;

  const bool uniform = pm->use_uniform_meshgen_fn_[X1DIR];
  xf.clear();
  nx = pm->mesh_size.nx1;
  for (int i = 0; i <= nx; ++i) {
    const Real rx = ComputeMeshGeneratorX(i, nx, uniform);
    xf.push_back(pm->MeshGenerator_[X1DIR](rx, pm->mesh_size));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Int2DOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//! \brief integrates the data over two dimensions and writes the resulting 1D array in
//!     the third.

void Int2DOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  // Open the output file for write.
  std::ofstream fout(fname, std::ios::out | std::ios::app | std::ios::binary);
  if (!fout.is_open()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in IntX1X2Output::ProcessHeader" << std::endl
        << "Unable to open output file '" << fname << "'" << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Write the current time.
  fout.write(reinterpret_cast<const char*>(&pm->time), sizeof(pm->time));

  // Allocate arrays for 2D integrals.
  AthenaArray<Real> area(pm->my_blocks(0)->ncells3);
  AthenaArray<Real> integral(num_vars_, nx);
  integral.ZeroClear();

  // Loop over meshblocks.
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    LoadOutputData(pmb);

    // Determine where the meshblock fits.
    const int offset = pmb->loc.lx3 * pmb->block_size.nx3 - pmb->ks;

    // Integrate each data field.
    int ii = 0;
    OutputData *pdata = pfirst_data_;
    while (pdata != nullptr) {
      const int nc = (pdata->type == "VECTORS") ? 3 : 1;
      for (int c = 0; c < nc; ++c) {
        for (int k = pmb->ks; k <= pmb->ke; ++k) {
          Real sum = 0;
          for (int j = pmb->js; j <= pmb->je; ++j) {
            pmb->pcoord->VolCenterFace3Area(k, j, pmb->is, pmb->ie, area);
            for (int i = pmb->is; i <= pmb->ie; ++i)
              sum += pdata->data(c,k,j,i) * area(i);
          }
          integral(ii, k + offset) += sum;
        }
        ++ii;
      }
      pdata = pdata->pnext;
    }

    ClearOutputData();
  }

  // Write the integrals.
  fout.write(reinterpret_cast<const char*>(integral.data()), integral.GetSizeInBytes());

  // Close and clean up.
  fout.close();
  integral.DeleteAthenaArray();
  area.DeleteAthenaArray();

  // Update output parameters.
  output_params.next_time += output_params.dt;
}
