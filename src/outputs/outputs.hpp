#ifndef OUTPUTS_HPP
#define OUTPUTS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file outputs.hpp
//  \brief provides classes to handle ALL types of data output
//======================================================================================

// C/C++ headers
#include <stdio.h>  // size_t
#include <string>  // string

// Athena++ headers
#include "wrapper.hpp"
#include "../athena.hpp"

#ifdef HDF5OUTPUT
#include <hdf5.h>
#endif

class Mesh;
class ParameterInput;

//! \struct OutputParameters
//  \brief  container for parameters read from <output> block in the input file

typedef struct OutputParameters {
  int block_number;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_type;
  std::string data_format;
  Real next_time, dt;
  int file_number;
  bool output_slice, output_sum, include_ghost_zones;
  int direction_of_sum;
  Real x1_slice, x2_slice, x3_slice;
} OutputParameters;

//! \struct OutputData
//  \brief

typedef struct OutputData {
  std::string type;        // one of (SCALARS,VECTORS)
  std::string name;
  AthenaArray<Real> data;  // array containing data (usually shallow copy/slice)
  struct OutputData *pnext, *pprev; // ptrs to next and previous nodes in list

  OutputData() : pnext(NULL), pprev(NULL) {};
} OutputData;

//--------------------------------------------------------------------------------------
//! \class OutputType
//  \brief abstract base class for different output types (modes).  Each OutputType
//  is designed to be a node in a linked list created and stored in the Outputs class.

class OutputType {
public:
  OutputType(OutputParameters oparams);
  virtual ~OutputType();

  // data
  int il,iu,jl,ju,kl,ku;          // OutputData array start/end indices
  OutputParameters output_params; // control data read from <output> block
  OutputType *pnext_type;         // ptr to next node in linked list of OutputTypes

  // functions
  void LoadOutputData(MeshBlock *pmb);
  void AppendOutputDataNode(OutputData *pod);
  void ClearOutputData();
  virtual void WriteOutputFile(Mesh *pm) = 0;   // pure virtual
  void FinalizeOutput(ParameterInput *pin);

protected:
  int num_vars_;             // number of variables in output
  OutputData *pfirst_data_;  // ptr to first OutputData in linked list
  OutputData *plast_data_;   // ptr to last OutputData in linked list
};

//! \class HistoryFile
//  \brief derived OutputType class for history dumps

class HistoryOutput : public OutputType {
public:
  HistoryOutput(OutputParameters oparams);
  ~HistoryOutput() {};
  void WriteOutputFile(Mesh *pm);
};

//! \class FormattedTableOutput
//  \brief derived OutputType class for formatted table (tabular) data

class FormattedTableOutput : public OutputType {
public:
  FormattedTableOutput(OutputParameters oparams);
  ~FormattedTableOutput() {};
  void WriteOutputFile(Mesh *pm);
};

//======================================================================================
// following should be deleted as they are updated

//typedef struct OutputData {
//  AthenaArray<Real> output_data;
//  struct OutputData *pnext;
//} OutputData;


//! \class VTKOutput
//  \brief derived OutputType class for vtk dumps

class VTKOutput : public OutputType {
public:
  VTKOutput(OutputParameters oparams);
  ~VTKOutput() {};

  void WriteOutputFile(OutputData *pod, MeshBlock *pmb);
};


//! \class RestartOutput
//  \brief derived OutputType class for restarting files

class RestartOutput : public OutputType {
private:
  IOWrapper resfile;
  IOWrapperSize_t listsize, headeroffset, datasize;
  char *data;
  int nbtotal, myns, mynb;

public:
  RestartOutput(OutputParameters oparams);
  ~RestartOutput() {};
  void Initialize(Mesh *pm, ParameterInput *pin, bool wtflag);
  void FinalizeOutput(ParameterInput *pin);
  void LoadOutputData(OutputData *pod, MeshBlock *pmb);
  void TransformOutputData(OutputData *pod, MeshBlock *pmb) {};
  void WriteOutputFile(OutputData *pod, MeshBlock *pmb) {};
};

#ifdef HDF5OUTPUT
//! \class ATHDF5Output
//  \brief derived OutputType class for Athena HDF5 files

class ATHDF5Output : public OutputType {
private:

  // Parameters
  static const int max_name_length = 20;  // maximum length of names excluding \0

  // HDF5 structures
  hid_t file;                                   // file to be written to
  hsize_t dims_start[5], dims_count[5];         // array sizes
  hid_t dataset_levels;                         // datasets to be written
  hid_t dataset_locations;
  hid_t dataset_x1f, dataset_x2f, dataset_x3f;
  hid_t *datasets_celldata;
  hid_t filespace_blocks;                       // local dataspaces for file
  hid_t filespace_blocks_3;
  hid_t filespace_blocks_nx1;
  hid_t filespace_blocks_nx2;
  hid_t filespace_blocks_nx3;
  hid_t *filespaces_vars_blocks_nx3_nx2_nx1;
  hid_t memspace_blocks;                        // local dataspaces for memory
  hid_t memspace_blocks_3;
  hid_t memspace_blocks_nx1;
  hid_t memspace_blocks_nx2;
  hid_t memspace_blocks_nx3;
  hid_t *memspaces_vars_blocks_nx3_nx2_nx1;
  hid_t property_list;                          // properties for writing

  // Metadata
  std::string filename;                       // name of athdf file
  int num_blocks_global;                      // number of MeshBlocks in simulation
  int num_blocks_local;                       // number of MeshBlocks on this Mesh
  int nx1, nx2, nx3;                          // sizes of MeshBlocks
  int is, ie, js, je, ks, ke;                 // indices for active zone
  int root_level;                             // number assigned to root level
  int num_datasets;                           // count of datasets to output
  int *num_variables;                         // list of counts of variables per dataset
  int num_total_variables;                    // total number of (scalar) variables
  char (*dataset_names)[max_name_length+1];   // array of C-string names of datasets
  char (*variable_names)[max_name_length+1];  // array of C-string names of variables
  int *levels_mesh;                           // array of refinement levels on Mesh
  long int *locations_mesh;                   // array of logical locations on Mesh
  float *x1f_mesh;                            // array of x1 values on Mesh
  float *x2f_mesh;                            // array of x1 values on Mesh
  float *x3f_mesh;                            // array of x1 values on Mesh
  float **data_buffers;                       // array of data buffers
  AthenaArray<Real> *data_arrays;             // array of slices into data

public:

  // Function declarations
  ATHDF5Output(OutputParameters oparams) : OutputType(oparams) {};
  ~ATHDF5Output() {};
  void Initialize(Mesh *pmesh, ParameterInput *pin, bool walltime_limit);
  void Finalize(ParameterInput *pin);
  void LoadOutputData(OutputData *pout_data, MeshBlock *pblock);
  void TransformOutputData(OutputData *pout_data, MeshBlock *pblock) {};
  void WriteOutputFile(OutputData *pout_data, MeshBlock *pblock);
  void MakeXDMF();
};
#endif

//======================================================================================



//--------------------------------------------------------------------------------------
//! \class Outputs
//  \brief root class for all Athena++ outputs.  Provides a linked list of OutputTypes,
//  with each node representing one mode of output to be made during a simulation.

class Outputs {
public:
  Outputs(Mesh *pm, ParameterInput *pin);
  ~Outputs();

  void MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag=false);

private:
  OutputType *pfirst_type_; // ptr to first OutputType in linked list
};
#endif
