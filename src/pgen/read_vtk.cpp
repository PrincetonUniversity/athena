//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file read_vtk.cpp
//! \brief problem generator, initalize mesh by read in vtk files.
//======================================================================================

// C headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// C++ headers
#include <algorithm>  // std::find()
#include <cinttypes>  // format macro "PRId64" for fixed-width integer type std::int64_t
#include <cstdio>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/string_utils.hpp" //split() and trim()

//function to read data field from vtk file
static void readvtk(MeshBlock *mb, std::string filename, std::string field,
                    int component, AthenaArray<Real> &data, int isjoinedvtk);
//swap bytes
static void ath_bswap(void *vdat, int len, int cnt);

//User defined boundary conditions
void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//Radiation boundary
namespace {
  AthenaArray<Real> G0_iang;
  Real G0, cr_rate;
} //namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SixRayBoundaryInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SixRayBoundaryOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x2, SixRayBoundaryInnerX2);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x2, SixRayBoundaryOuterX2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, SixRayBoundaryInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, SixRayBoundaryOuterX3);
  G0 = pin->GetOrAddReal("chem_radiation", "G0", 0.);
  G0_iang.NewAthenaArray(6);
  G0_iang(BoundaryFace::inner_x1) = pin->GetOrAddReal("chem_radiation","G0_inner_x1",G0);
  G0_iang(BoundaryFace::inner_x2) = pin->GetOrAddReal("chem_radiation","G0_inner_x2",G0);
  G0_iang(BoundaryFace::inner_x3) = pin->GetOrAddReal("chem_radiation","G0_inner_x3",G0);
  G0_iang(BoundaryFace::outer_x1) = pin->GetOrAddReal("chem_radiation","G0_outer_x1",G0);
  G0_iang(BoundaryFace::outer_x2) = pin->GetOrAddReal("chem_radiation","G0_outer_x2",G0);
  G0_iang(BoundaryFace::outer_x3) = pin->GetOrAddReal("chem_radiation","G0_outer_x3",G0);
  cr_rate = pin->GetOrAddReal("chem_radiation", "CR", 2e-16);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem by reading in vtk file.
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
  //dimensions of mesh
  const int Nx_mesh = pmy_mesh->mesh_size.nx1;
  const int Ny_mesh = pmy_mesh->mesh_size.nx2;
  const int Nz_mesh = pmy_mesh->mesh_size.nx3;

  //read joined or unjoined vtk files
  int isjoinedvtk = pin->GetOrAddInteger("problem", "is_joined_vtk", 0);
  AthenaArray<Real> data; //temporary array to store data;
  AthenaArray<Real> b; //needed for PrimitiveToConserved()
  if (isjoinedvtk != 0) {
    if (DEBUG) {
      printf("Joined vtk file. data size = (%d, %d, %d)\n",
             Nz_mesh, Ny_mesh, Nx_mesh);
    }
    data.NewAthenaArray(Nz_mesh, Ny_mesh, Nx_mesh);
    b.NewAthenaArray(Nz_mesh, Ny_mesh, Nx_mesh);
  } else {
    if (DEBUG) {
      printf("Unjoined vtk file. data size = (%d, %d, %d)\n",
             Nz, Ny, Nx);
    }
    data.NewAthenaArray(Nz, Ny, Nx);
    b.NewAthenaArray(Nz, Ny, Nx);
  }
  std::stringstream msg; // error message
  std::string vtkfile; // corresponding vtk file for this meshblock
  // initial abundance
  const Real r_init = pin->GetReal("problem", "r_init");

  // parse input parameters
  std::string vtkfile0 = pin->GetString("problem", "vtkfile");//id0 file
  std::string str_scalars = pin->GetString("problem", "scalars");
  std::string str_vectors = pin->GetString("problem", "vectors");
  std::vector<std::string> scalar_fields = StringUtils::split(str_scalars, ',');
  std::vector<std::string> vector_fields = StringUtils::split(str_vectors, ',');

  if (isjoinedvtk != 0) {
    // for joined vtk file, read with processor 0, and broadcast
    // TODO(KGF): risk of narrowing here:
    int gis = static_cast<int>(loc.lx1) * Nx;
    int gjs = static_cast<int>(loc.lx2) * Ny;
    int gks = static_cast<int>(loc.lx3) * Nz;
    vtkfile = vtkfile0;
    // diagnostic printing of filename
    if (Globals::my_rank == 0) {
      printf("meshblock gid=%d, location = (%"
             PRId64 " %" PRId64 " %" PRId64") level=%d, vtk file = %s\n",
             gid, loc.lx1, loc.lx2, loc.lx3, loc.level, vtkfile.c_str());
    }
    // scalars
    for(std::uint64_t i = 0; i < scalar_fields.size(); ++i) {
      if (scalar_fields[i] == "density") {
        if (Globals::my_rank == 0) {
          if (DEBUG) {
            printf("Process 0: start to read density.\n");
          }
          readvtk(this, vtkfile, "density", 0, data, isjoinedvtk);
        }
#ifdef MPI_PARALLEL
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("Start to MPI broadcast density.\n");
        }
        ierr = MPI_Bcast(data.data(), Nx_mesh*Ny_mesh*Nz_mesh,
                         MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("End MPI brocast density.\n");
        }
#endif // MPI_PARALLEL
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IDN, k, j, i) = data(k-ks+gks, j-js+gjs, i-is+gis);
            }
          }
        }
      } else if (scalar_fields[i] == "pressure") {
        if (Globals::my_rank == 0) {
          if (DEBUG) {
            printf("Process 0: Start to read pressure.\n");
          }
          readvtk(this, vtkfile, "pressure", 0, data, isjoinedvtk);
        }
#ifdef MPI_PARALLEL
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("Start to MPI broadcast pressure.\n");
        }
        ierr = MPI_Bcast(data.data(), Nx_mesh*Ny_mesh*Nz_mesh,
                         MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("End MPI brocast pressure.\n");
        }
#endif // MPI_PARRLLEL
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IEN, k, j, i) = data(k-ks+gks, j-js+gjs, i-is+gis);
            }
          }
        }
      } else {
        msg << "### FATAL ERROR in Problem Generator [ProblemGenerator]" << std::endl
          << "Scalar field not recognized: " << scalar_fields[i] << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    // vectors
    for(std::uint64_t i = 0; i < vector_fields.size(); ++i) {
      if (vector_fields[i] == "velocity") {
        // vx
        if (Globals::my_rank == 0) {
          if (DEBUG) {
            printf("Process 0: Start to read velocity.\n");
          }
          readvtk(this, vtkfile, "velocity", 0, data,isjoinedvtk);
        }
#ifdef MPI_PARALLEL
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("Start to MPI broadcast velocity 0.\n");
        }
        ierr = MPI_Bcast(data.data(), Nx_mesh*Ny_mesh*Nz_mesh,
                         MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("End MPI brocast velocity 0.\n");
        }
#endif // MPI_PARRLLEL
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IVX, k, j, i) = data(k-ks+gks, j-js+gjs, i-is+gis);
            }
          }
        }
        // vy
        if (Globals::my_rank == 0) {
          readvtk(this, vtkfile, "velocity", 1, data,isjoinedvtk);
        }
#ifdef MPI_PARALLEL
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("Start to MPI broadcast velocity 1.\n");
        }
        ierr = MPI_Bcast(data.data(), Nx_mesh*Ny_mesh*Nz_mesh,
                         MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("End MPI brocast velocity 1.\n");
        }
#endif // MPI_PARRLLEL
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IVY, k, j, i) = data(k-ks+gks, j-js+gjs, i-is+gis);
            }
          }
        }
        // vz
        if (Globals::my_rank == 0) {
          readvtk(this, vtkfile, "velocity", 2, data,isjoinedvtk);
        }
#ifdef MPI_PARALLEL
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("Start to MPI broadcast velocity 2.\n");
        }
        ierr = MPI_Bcast(data.data(), Nx_mesh*Ny_mesh*Nz_mesh,
                         MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
        if (DEBUG) {
          MPI_Barrier(MPI_COMM_WORLD);
          printf("End MPI brocast velocity 2.\n");
        }
#endif // MPI_PARRLLEL
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IVZ, k, j, i) = data(k-ks+gks, j-js+gjs, i-is+gis);
            }
          }
        }
      } else {
        msg << "### FATAL ERROR in Problem Generator [ProblemGenerator]" << std::endl
          << "Scalar field not recognized: " << scalar_fields[i] << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }

  } else {
    // find coresponding filename.
    if (loc.lx1 == 0 && loc.lx2 == 0 && loc.lx3 == 0) {
      vtkfile = vtkfile0;
    } else {
      // find the corespoinding athena4.2 global id
      int64_t id_old = loc.lx1 + loc.lx2 * pmy_mesh->nrbx1
        + loc.lx3 * pmy_mesh->nrbx1 * pmy_mesh->nrbx2;
      // get vtk file name .../id#/problem-id#.????.vtk
      std::stringstream id_str_stream;
      id_str_stream << "id" << id_old; // id#
      std::string id_str = id_str_stream.str();
      std::size_t pos1 = vtkfile0.find_last_of('/'); //last /
      std::size_t pos2 = vtkfile0.find_last_of('/', pos1-1); //second last /
      std::string base_dir = vtkfile0.substr(0, pos2+1); // "base_directory/"
      std::string vtk_name0 = vtkfile0.substr(pos1); // "/bala.????.vtk"
      std::size_t pos3 = vtk_name0.find_first_of('.');
      std::string vtk_name = vtk_name0.substr(0, pos3) + "-" + id_str
        + vtk_name0.substr(pos3);
      std::cout << id_str << ", " << base_dir << ", " << vtk_name << std::endl;
      vtkfile = base_dir + id_str + vtk_name;
    }

    // diagnostic printing of filename
    printf("meshblock gid=%d, location = (%"
           PRId64 " %" PRId64 " %" PRId64") level=%d, vtk file = %s\n",
           gid, loc.lx1, loc.lx2, loc.lx3, loc.level, vtkfile.c_str());

    // read scalars
    for(std::uint64_t i = 0; i < scalar_fields.size(); ++i) {
      if (scalar_fields[i] == "density") {
        readvtk(this, vtkfile, "density", 0, data,isjoinedvtk);
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IDN, k, j, i) = data(k-ks, j-js, i-is);
            }
          }
        }
      } else if (scalar_fields[i] == "pressure") {
        readvtk(this, vtkfile, "pressure", 0, data,isjoinedvtk);
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IEN, k, j, i) = data(k-ks, j-js, i-is);
            }
          }
        }
      } else {
        msg << "### FATAL ERROR in Problem Generator [ProblemGenerator]" << std::endl
          << "Scalar field not recognized: " << scalar_fields[i] << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    // read vectors
    for (std::uint64_t i = 0; i < vector_fields.size(); ++i) {
      if (vector_fields[i] == "velocity") {
        // vx
        readvtk(this, vtkfile, "velocity", 0, data,isjoinedvtk);
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IVX, k, j, i) = data(k-ks, j-js, i-is);
            }
          }
        }
        // vy
        readvtk(this, vtkfile, "velocity", 1, data,isjoinedvtk);
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IVY, k, j, i) = data(k-ks, j-js, i-is);
            }
          }
        }
        // vz
        readvtk(this, vtkfile, "velocity", 2, data,isjoinedvtk);
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->w(IVZ, k, j, i) = data(k-ks, j-js, i-is);
            }
          }
        }
      } else {
        msg << "### FATAL ERROR in Problem Generator [ProblemGenerator]" << std::endl
          << "Scalar field not recognized: " << scalar_fields[i] << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
  }

  // change primative variables to conservative variables.
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord,
                                     is, ie, js, je, ks, ke);

  // intialize radiation field
  if (CHEMRADIATION_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifreq=0; ifreq < pchemrad->nfreq; ++ifreq) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              pchemrad->ir(k, j, i, ifreq * pchemrad->nang + iang) = G0_iang(iang);
            }
          }
          if (CHEMISTRY_ENABLED) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              // cr rate
              pchemrad->ir(k, j, i,
                  pscalars->chemnet.index_cr_ * pchemrad->nang + iang) = cr_rate;
            }
          }
        }
      }
    }
    // calculate the average radiation field for output of the initial condition
    pchemrad->pchemradintegrator->CopyToOutput();
  }

  // intialize chemical species
  if (NSPECIES > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            pscalars->s(ispec, k, j, i) = r_init * phydro->u(IDN, k, j, i);
            if (CHEMISTRY_ENABLED) {
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec * phydro->u(IDN, k, j, i);
              }
            }
          }
        }
      }
    }
  }

  return;
}


//======================================================================================
//! \fn static void readvtk(MeshBlock *, std::string , std::string,
//!                         int , AthenaArray<Real> &, int)
//! \brief read from athena4.2 vtk file and assign initial value to
//! variables.
//======================================================================================
static void readvtk(MeshBlock *mb, std::string filename, std::string field,
                    int component, AthenaArray<Real> &data, int isjoinedvtk) {
  std::stringstream msg;
  FILE *fp = NULL;
  char cline[256], type[256], variable[256], format[256], t_type[256], t_format[256];
  std::string line;
  const std::string athena_header = "# vtk DataFile Version 2.0"; // athena4.2 header
  const std::string athena_header3 = "# vtk DataFile Version 3.0"; // athena4.2 header
  bool SHOW_OUTPUT = false;
  int Nx_vtk, Ny_vtk, Nz_vtk; // dimensions of vtk files
  // dimensions of meshblock
  int Nx_mb, Ny_mb, Nz_mb;
  if (isjoinedvtk != 0) {
    Nx_mb = mb->pmy_mesh->mesh_size.nx1;
    Ny_mb = mb->pmy_mesh->mesh_size.nx2;
    Nz_mb = mb->pmy_mesh->mesh_size.nx3;
  } else {
    Nx_mb = mb->ie - mb->is + 1;
    Ny_mb = mb->je - mb->js + 1;
    Nz_mb = mb->ke - mb->ks + 1;
  }
  double ox_vtk, oy_vtk, oz_vtk; // origins of vtk file
  double dx_vtk, dy_vtk, dz_vtk; // spacings of vtk file
  int cell_dat_vtk; // total number of cells in vtk file
  // total number of cells in MeshBlock
  const int cell_dat_mb = Nx_mb * Ny_mb * Nz_mb;
  int retval;
  size_t nread; // file handler return value
  float fdat, fvec[3]; // , ften[9];// store float format scalar, vector, and tensor

  if ( (fp = fopen(filename.c_str(),"r")) == NULL ) {
    msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
      << "Unable to open file: " << filename << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // get header
  fgets(cline,256,fp);
  line.assign(cline);
  StringUtils::trim(line);
  if (SHOW_OUTPUT) {
    std::cout << line << std::endl;
  }
  if (line != athena_header && line != athena_header3) {
    fclose(fp);
    msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
      << "Assuming Athena4.2 header " << athena_header << ", get header "
      << line << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // get comment field
  fgets(cline,256,fp);
  line.assign(cline);
  StringUtils::trim(line);
  if (SHOW_OUTPUT) {
    std::cout << line << std::endl;
  }

  // get BINARY or ASCII
  fgets(cline,256,fp);
  line.assign(cline);
  StringUtils::trim(line);
  if (SHOW_OUTPUT) {
    std::cout << line << std::endl;
  }
  if (line != "BINARY") {
    fclose(fp);
    msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
      << "Unsupported file format: " << line << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // get DATASET STRUCTURED_POINTS or DATASET UNSTRUCTURED_GRID
  fgets(cline,256,fp);
  line.assign(cline);
  StringUtils::trim(line);
  if (SHOW_OUTPUT) {
    std::cout << line << std::endl;
  }
  if (line != "DATASET STRUCTURED_POINTS") {
    fclose(fp);
    msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
      << "Unsupported file data set structure: " << line << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // I'm assuming from this point on that the header is in good shape

  // read dimensions
  fgets(cline,256,fp);
  if (SHOW_OUTPUT) {
    std::cout << cline;
  }
  sscanf(cline,"DIMENSIONS %d %d %d\n",&Nx_vtk,&Ny_vtk,&Nz_vtk);
  // We want to store the number of grid cells, not the number of grid
  // cell corners.
  Nx_vtk--;
  Ny_vtk--;
  Nz_vtk--;
  if (Nx_vtk != Nx_mb || Ny_vtk != Ny_mb || Nz_vtk != Nz_mb) {
    fclose(fp);
    msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
      << "Dimensions of VTK file" << filename
      << " is (" << Nx_vtk << ", " << Ny_vtk << ", " << Nz_vtk
      << "), does not match the dimensions of meshblock ("
      << Nx_mb << ", " << Ny_mb << ", " << Nz_mb << ")." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Origin
  fgets(cline,256,fp);
  if (SHOW_OUTPUT) {
    std::cout << cline;
  }
  sscanf(cline,"ORIGIN %le %le %le\n",&ox_vtk,&oy_vtk,&oz_vtk);

  // spacing
  fgets(cline,256,fp);
  if (SHOW_OUTPUT) {
    std::cout << cline;
  }
  sscanf(cline,"ORIGIN %le %le %le\n",&dx_vtk,&dy_vtk,&dz_vtk);

  // Cell Data = Nx*Ny*Nz
  fgets(cline,256,fp);
  if (SHOW_OUTPUT) {
    std::cout << cline;
  }
  sscanf(cline,"CELL_DATA %d\n",&cell_dat_vtk);
  if (cell_dat_vtk != cell_dat_mb) {
    fclose(fp);
    msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
      << "Cell data in vtk file " << filename
      << " is " << cell_dat_vtk << ", does not match the cell data of meshblock "
      << cell_dat_mb << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // -------------read data--------------
  while (true) {
    retval = std::fscanf(fp,"%s %s %s\n", type, variable, format);
    if (retval == EOF) { // Assuming no errors, we are done.
      fclose(fp); // close file
      return;
    }
    if (SHOW_OUTPUT) {
      printf("%s %s %s\n", type, variable ,format);
    }
    // check format
    if (std::strcmp(format, "float") != 0) {
      fclose(fp);
      msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
        << "expected  \"float\" format, found " << type << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    // check lookup table
    if (std::strcmp(type, "SCALARS") == 0) {
      // Read in the LOOKUP_TABLE (only default supported for now)
      fscanf(fp,"%s %s\n", t_type, t_format);
      if (std::strcmp(t_type, "LOOKUP_TABLE") != 0 ||
          std::strcmp(t_format, "default") != 0 ) {
        fclose(fp);
        msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
          << "Expected \"LOOKUP_TABLE default, found "
          << t_type << " " << t_format << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if (SHOW_OUTPUT) {
        printf("%s %s\n", t_type, t_format);
      }
    }

    // determine variable type and read data
    // read scalars
    if (std::strcmp(type, "SCALARS") == 0) {
      if (std::strcmp(variable, field.c_str()) == 0) {
        printf("  Reading %s...\n", variable);
        for (int k=0; k<Nz_vtk; k++) {
          for (int j=0; j<Ny_vtk; j++) {
            for (int i=0; i<Nx_vtk; i++) {
              if ((nread = std::fread(&fdat, sizeof(float), 1, fp)) != 1) {
                fclose(fp);
                msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
                  << "Error reading SCALARS... " << std::endl;
                throw std::runtime_error(msg.str().c_str());
              }
              ath_bswap(&fdat, sizeof(float), 1);
              data(k, j, i) = fdat;
            }
          }
        }
        fclose(fp);
        return;
      } else {
        if (SHOW_OUTPUT) printf("  Skipping %s...\n",variable);
        std::fseek(fp, cell_dat_vtk*sizeof(float), SEEK_CUR);
      }
    // read vectors
    } else if (std::strcmp(type, "VECTORS") == 0) {
      if (std::strcmp(variable, field.c_str()) == 0) {
        printf("  Reading %s%d...\n", variable, component);
        for (int k=0; k<Nz_vtk; k++) {
          for (int j=0; j<Ny_vtk; j++) {
            for (int i=0; i<Nx_vtk; i++) {
              if ((nread = std::fread(&fvec, sizeof(float), 3, fp)) != 3) {
                fclose(fp);
                msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
                  << "Error reading VECTORS... " << std::endl;
                throw std::runtime_error(msg.str().c_str());
              }
              ath_bswap(&fvec, sizeof(float), 3);
              data(k, j, i) = fvec[component];
            }
          }
        }
        fclose(fp);
        return;
      } else {
        if (SHOW_OUTPUT) printf("  Skipping %s...\n", variable);
        std::fseek(fp, 3*cell_dat_vtk*sizeof(float), SEEK_CUR);
      }
    // read tensors, not supported yet
    } else if (std::strcmp(type, "TENSORS") == 0) {
      if (std::strcmp(variable, field.c_str()) == 0) {
        fclose(fp);
        msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
          << "TENSORS reading not supported." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      } else {
        if (SHOW_OUTPUT) printf("  Skipping %s...\n", variable);
        std::fseek(fp, 9*cell_dat_vtk*sizeof(float), SEEK_CUR);
      }
    } else {
        fclose(fp);
        msg << "### FATAL ERROR in Problem Generator [read_vtk]" << std::endl
          << "Input type not supported: " << type << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }
}


//======================================================================================
//! \fn static void ath_bswap(void *vdat, int len, int cnt)
//! \brief Swap bytes, code stolen from Athena4.2, NEMO
//======================================================================================
static void ath_bswap(void *vdat, int len, int cnt) {
  char tmp, *dat = static_cast<char *>(vdat);
  int k;

  if (len==1) {
    return;
  } else if (len==2) {
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
      dat += 2;
    }
  } else if (len==4) {
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
      tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
      dat += 4;
    }
  } else if (len==8) {
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
      tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
      tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
      tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
      dat += 8;
    }
  } else {  /* the general SLOOOOOOOOOW case */
    for(k=0; k<len/2; k++) {
      tmp = dat[k];
      dat[k] = dat[len-1-k];
      dat[len-1-k] = tmp;
    }
  }
}

void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          pmb->pscalars->s(n,k,j,il-i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,il-i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,k,jl-j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,jl-j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,kl-k,j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,kl-k,j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          pmb->pscalars->s(n,k,j,iu+i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,iu+i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,k,ju+j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,ju+j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,ku+k,j,i) = 0;
        }
      }
    }
  }
  // set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          prim(n,ku+k,j,i) = 0;
        }
      }
    }
  }
  return;
}
