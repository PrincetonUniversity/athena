#ifndef OUTPUTS_HPP
#define OUTPUTS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file outputs.hpp
 *  \brief provides multiple classes to handle ALL types of data output (fluid, field,
 *  radiation, particles, etc.)
 *====================================================================================*/

class Mesh;
class ParameterInput;

struct OutputDataHeader {
  std::string descriptor;
  int il,iu,jl,ju,kl,ku;
};

struct OutputDataNodeHeader {
  std::string type;
  std::string name;
};

class OutputDataNode {
public:
  OutputDataNode(AthenaArray<Real> *parray, OutputDataNodeHeader head);
  ~OutputDataNode();

  OutputDataNodeHeader header;
  AthenaArray<Real> *pdata;

  OutputDataNode *pnext, *pprev;
};

class OutputData {
public:
  OutputData();
  ~OutputData();

  OutputDataHeader header;
  void AppendNode(AthenaArray<Real> *parray, OutputDataNodeHeader head);
  void ReplaceNode(OutputDataNode *pold, OutputDataNode *pnew);

  OutputDataNode *pfirst_node;  // Pointer to first node
  OutputDataNode *plast_node;   // Pointer to last node
};

//! \struct OutputBlock
//  \brief  contains parameter values read from <output> blocks in the input file

struct OutputBlock {
  Real next_time, dt;
  int block_number;
  int file_number;
  int islice, jslice, kslice;
  int isum, jsum, ksum;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_format;
  std::string data_format;
};

//-------------------------- OutputTypes base and derived classes ----------------------
//! \class OutputType
//  \brief abstract base class for different output types.

class OutputType {
public:
  OutputType(OutputBlock out_blck, Block *pb);
  ~OutputType();

  Block *pparent_block;

  virtual OutputData* LoadOutputData();
  virtual void ComputeOutputData(OutputData *pod);
  virtual void WriteOutputData() = 0;  // pure virtual function!

  void Slice(OutputData* pod, int dim);

  OutputBlock output_block;

  OutputType *pnext;
};

//! \class OutputList
//  \brief provides a linked list of OutputTypes, with each node representing one type
//  of output to be made during a simulation.

class OutputList {
public:
  OutputList(Block *pb);
  ~OutputList();

  Block *pparent_block;

  void InitOutputs(ParameterInput *pin);
  void MakeOutputs();

private:
  OutputType *pfirst_out_;
};

//! \class FormattedTableOutput
//  \brief derived OutputType class for formatted table (tabular) data

class FormattedTableOutput : public OutputType {
public:
  FormattedTableOutput(OutputBlock out_blk, Block *pb);
  ~FormattedTableOutput() {};

  void WriteOutputData();

private:
};

//! \class HistoryOutput
//  \brief derived OutputType class for history dumps

class HistoryOutput : public OutputType {
public:
  HistoryOutput(OutputBlock out_blk, Block *pb);
  ~HistoryOutput() {};

  void ComputeOutputData(OutputData *pod);
  void WriteOutputData();
};

#endif
