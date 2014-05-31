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

struct DataListHeader {
  std::string descriptor;
  int il,iu,jl,ju,kl,ku;
};

struct DataNodeHeader {
  std::string type;
  std::string name;
};

class DataNode {
public:
  DataNode(AthenaArray<Real> *pinit_data, DataNodeHeader init_head);
  ~DataNode();

  DataNodeHeader header;
  AthenaArray<Real> *pdata;

  DataNode *pnext;
};

class DataList {
public:
  DataList();
  ~DataList();

  DataListHeader header;
  void InsertNode(AthenaArray<Real> *pdata, DataNodeHeader init_head);

  DataNode *pfirst_node;  // Pointer to first node
private:
  DataNode *plast_node_;   // Pointer to last node
};

//! \struct OutputBlock
//  \brief  contains parameter values read from <output> blocks in the input file

struct OutputBlock {
  Real last_time, dt;
  int block_number;
  int file_number;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_format;
  std::string data_format;
};

//! \class OutputType
//  \brief abstract base class for different output types.  The WriteOutputData pure
//  virtual function is overloaded in each of the derived output classes below.

class OutputType {
public:
  OutputType(OutputBlock out_blck, Block *pb);
  ~OutputType();

  virtual DataList* LoadDataList();
  virtual void ComputeDataList();
  virtual void WriteOutputData() = 0;  // pure virtual function!

  OutputBlock output_block;
  Block *pparent_block;
  OutputType *pnext;

  DataList *pdlist;
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

  void ComputeDataList();
  void WriteOutputData();

};
#endif
