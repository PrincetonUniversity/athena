#ifndef OUTPUTS_HPP
#define OUTPUTS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file outputs.hpp
//  \brief provides multiple classes to handle ALL types of data output (fluid, bfield,
//  gravity, radiation, particles, etc.)
//======================================================================================

#include <stdio.h> // size_t
#include "wrapper.hpp"

class Mesh;
class ParameterInput;

//! \struct OutputParameters
//  \brief  control parameter values read from <output> block in the input file

typedef struct OutputParameters {
  Real next_time, dt;
  int block_number;
  int file_number;
  int islice, jslice, kslice;
  Real x1_slice, x2_slice, x3_slice;
  int isum, jsum, ksum;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_type;
  std::string data_format;
} OutputParameters;

//! \struct OutputVariable
//  \brief node in a linked list of output variables in OutputData class

class OutputVariable {
public:
  OutputVariable();
  ~OutputVariable();

  std::string type; // one of (SCALARS,VECTORS)
  std::string name;
  AthenaArray<Real> data;  // array containing data (may be shallow copy/slice)
  OutputVariable *pnext, *pprev; // ptrs to next and previous nodes in list
};

//! \struct OutputDataHeader
//  \brief metadata describing OutputData contained within an OutputType

typedef struct OutputDataHeader {
  std::string descriptor; // time, cycle, variables in output
  std::string transforms; // list of any transforms (sum, slice) applied to variables
  int il,iu,jl,ju,kl,ku;  // range of data arrays
  int ndata;              // number of data points in arrays
} OutputDataHeader;

//! \class OutputData
//  \brief container for output data, composed of a header (metadata describing output),
//  a linked list of OutputVariables, and functions to add and replace nodes

class OutputData {
public:
  OutputData();
  ~OutputData();

  OutputDataHeader data_header;
  OutputVariable *pfirst_var;  // ptr to first variable (node) in linked list
  OutputVariable *plast_var;   // ptr to last  variable (node) in linked list

  void AppendNode(OutputVariable *pvar);
  void ReplaceNode(OutputVariable *pold, OutputVariable *pnew);
};


//-------------------------- OutputTypes base and derived classes ----------------------
//! \class OutputType
//  \brief abstract base class for different output types (modes), designed to be a node
//  in a linked list created and stored in objects of the Outputs class.

class OutputType {
public:
  OutputType(OutputParameters oparams);
  ~OutputType();
  OutputParameters output_params; // control data read from <output> block 

// functions that operate on OutputData container

  virtual void Initialize(Mesh *pM, ParameterInput *pin) {};
  virtual void Finalize(void) {};
  virtual void LoadOutputData(OutputData *pod, MeshBlock *pmb);
  virtual void TransformOutputData(OutputData *pod, MeshBlock *pmb);
  virtual void WriteOutputFile(OutputData *pod, MeshBlock *pmb) = 0;
          // pure virtual!

// functions that implement useful transforms applied to each variable in OutputData

  void Slice(OutputData* pod, MeshBlock *pmb, int dim);
  void Sum(OutputData* pod, MeshBlock *pmb, int dim);

  OutputType *pnext_type;   // ptr to next node in linked list of OutputTypes
};

//! \class FormattedTableOutput
//  \brief derived OutputType class for formatted table (tabular) data

class FormattedTableOutput : public OutputType {
public:
  FormattedTableOutput(OutputParameters oparams);
  ~FormattedTableOutput() {};

  void WriteOutputFile(OutputData *pod, MeshBlock *pmb);

private:
};

//! \class HistoryOutput
//  \brief derived OutputType class for history dumps

class HistoryOutput : public OutputType {
public:
  HistoryOutput(OutputParameters oparams);
  ~HistoryOutput() {};

  void LoadOutputData(OutputData *pod, MeshBlock *pmb); // overloads base class function
  void WriteOutputFile(OutputData *pod, MeshBlock *pmb);
};

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
  ResFile resfile;
  ResSize_t *blocksize, *offset;

public:
  RestartOutput(OutputParameters oparams);
  ~RestartOutput() {};
  void Initialize(Mesh *pm, ParameterInput *pin);
  void Finalize(void);
  void LoadOutputData(OutputData *pod, MeshBlock *pmb) {};
  void TransformOutputData(OutputData *pod, MeshBlock *pmb) {};

  void WriteOutputFile(OutputData *pod, MeshBlock *pmb);
};

//--------------------- end of OutputTypes base and derived classes --------------------

//! \class Outputs
//  \brief root class for all Athena++ outputs.  Provides a linked list of OutputTypes,
//  with each node representing one mode of output to be made during a simulation.

class Outputs {
public:
  Outputs(Mesh *pm, ParameterInput *pin);
  ~Outputs();

  void MakeOutputs(Mesh *pm, ParameterInput *pin);

private:
  OutputType *pfirst_type_; // ptr to first OutputType in linked list
};
#endif
