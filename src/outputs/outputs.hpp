#ifndef OUTPUTS_OUTPUTS_HPP_
#define OUTPUTS_OUTPUTS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outputs.hpp
//! \brief provides classes to handle ALL types of data output

// C headers

// C++ headers
#include <cstdio>  // std::size_t
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "io_wrapper.hpp"

#ifdef HDF5OUTPUT
#include <hdf5.h>
#endif

// forward declarations
class Mesh;
class ParameterInput;
class Coordinates;

//----------------------------------------------------------------------------------------
//! \struct OutputParameters
//! \brief  container for parameters read from `<output>` block in the input file

struct OutputParameters {
  int block_number;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_type;
  std::string data_format;
  Real next_time, dt;
  int dcycle;
  int file_number;
  bool output_slicex1, output_slicex2, output_slicex3;
  bool output_sumx1, output_sumx2, output_sumx3;
  bool include_ghost_zones, cartesian_vector;
  bool orbital_system_output;
  int islice, jslice, kslice;
  Real x1_slice, x2_slice, x3_slice;
  // TODO(felker): some of the parameters in this class are not initialized in constructor
  OutputParameters() : block_number(0), next_time(0.0), dt(0.0), file_number(0),
                       output_slicex1(false),output_slicex2(false),output_slicex3(false),
                       output_sumx1(false), output_sumx2(false), output_sumx3(false),
                       include_ghost_zones(false), cartesian_vector(false),
                       islice(0), jslice(0), kslice(0) {}
};

//----------------------------------------------------------------------------------------
//! \struct OutputData
//! \brief container for output data and metadata; node in nested doubly linked list

struct OutputData {
  std::string type;        // one of (SCALARS,VECTORS) used for vtk outputs
  std::string name;
  AthenaArray<Real> data;  // array containing data (usually shallow copy/slice)
  // ptrs to previous and next nodes in doubly linked list:
  OutputData *pnext, *pprev;

  OutputData() : pnext(nullptr),  pprev(nullptr) {}
};

//----------------------------------------------------------------------------------------
//! \brief abstract base class for different output types (modes/formats). Each OutputType
//! is designed to be a node in a singly linked list created & stored in the Outputs class

class OutputType {
 public:
  // mark single parameter constructors as "explicit" to prevent them from acting as
  // implicit conversion functions: for f(OutputType arg), prevent f(anOutputParameters)
  explicit OutputType(OutputParameters oparams);

  // rule of five:
  virtual ~OutputType() = default;
  // copy constructor and assignment operator (pnext_type, pfirst_data, etc. are shallow
  // copied)
  OutputType(const OutputType& copy_other) = default;
  OutputType& operator=(const OutputType& copy_other) = default;
  // move constructor and assignment operator
  OutputType(OutputType&&) = default;
  OutputType& operator=(OutputType&&) = default;

  // data
  int out_is, out_ie, out_js, out_je, out_ks, out_ke;  // OutputData array start/end index
  OutputParameters output_params; // control data read from <output> block
  OutputType *pnext_type;         // ptr to next node in singly linked list of OutputTypes

  // functions
  void LoadOutputData(MeshBlock *pmb);
  void AppendOutputDataNode(OutputData *pdata);
  void ReplaceOutputDataNode(OutputData *pold, OutputData *pnew);
  void ClearOutputData();
  bool TransformOutputData(MeshBlock *pmb);
  bool SliceOutputData(MeshBlock *pmb, int dim);
  void SumOutputData(MeshBlock *pmb, int dim);
  void CalculateCartesianVector(AthenaArray<Real> &src, AthenaArray<Real> &dst,
                                Coordinates *pco);
  bool ContainVariable(const std::string &haystack, const std::string &needle);
  // following pure virtual function must be implemented in all derived classes
  virtual void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) = 0;

 protected:
  int num_vars_;             // number of variables in output
  // nested doubly linked list of OutputData nodes (of the same OutputType):
  OutputData *pfirst_data_;  // ptr to head OutputData node in doubly linked list
  OutputData *plast_data_;   // ptr to tail OutputData node in doubly linked list
};

//----------------------------------------------------------------------------------------
//! \class HistoryOutput
//! \brief derived OutputType class for history dumps

class HistoryOutput : public OutputType {
 public:
  explicit HistoryOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//----------------------------------------------------------------------------------------
//! \class FormattedTableOutput
//! \brief derived OutputType class for formatted table (tabular) data

class FormattedTableOutput : public OutputType {
 public:
  explicit FormattedTableOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//----------------------------------------------------------------------------------------
//! \class VTKOutput
//! \brief derived OutputType class for vtk dumps

class VTKOutput : public OutputType {
 public:
  explicit VTKOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//----------------------------------------------------------------------------------------
//! \class RestartOutput
//! \brief derived OutputType class for restart dumps

class RestartOutput : public OutputType {
 public:
  explicit RestartOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

#ifdef HDF5OUTPUT
//----------------------------------------------------------------------------------------
//! \class ATHDF5Output
//! \brief derived OutputType class for Athena HDF5 files

class ATHDF5Output : public OutputType {
 public:
  // Function declarations
  explicit ATHDF5Output(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
  void MakeXDMF();

 private:
  // Parameters
  static const int max_name_length = 20;  // maximum length of names excluding \0

  // Metadata
  std::string filename;                       // name of athdf file
  float code_time;                            // time in code unit for XDMF
  int num_blocks_global;                      // number of MeshBlocks in simulation
  int nx1, nx2, nx3;                          // sizes of MeshBlocks
  int num_datasets;                           // count of datasets to output
  int *num_variables;                         // list of counts of variables per dataset
  char (*dataset_names)[max_name_length+1];   // array of C-string names of datasets
  char (*variable_names)[max_name_length+1];  // array of C-string names of variables
};
#endif

//----------------------------------------------------------------------------------------
//! \class Outputs
//! \brief root class for all Athena++ outputs. Provides a singly linked list of
//! OutputTypes, with each node representing one mode/format of output to be made.

class Outputs {
 public:
  Outputs(Mesh *pm, ParameterInput *pin);
  ~Outputs();

  void MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag=false);

 private:
  OutputType *pfirst_type_; // ptr to head OutputType node in singly linked list
  // (not storing a reference to the tail node)
};
#endif // OUTPUTS_OUTPUTS_HPP_
