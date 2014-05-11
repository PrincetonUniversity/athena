#ifndef DATA_OUTPUT_HPP
#define DATA_OUTPUT_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file data_output.hpp
 *  \brief provides class to handle all types of data output
 *====================================================================================*/

class ParameterInput;
class Mesh;

//! \struct OutputBlock
//  \brief  node in a linked list of output blocks in the input file

typedef struct OutputBlock {
  int block_number;
  std::string block_name;
  std::string filename;
  Real last_time, dt;
  struct OutputBlock *pnext;
} OutputBlock;

//! \class DataOutput
//  \brief output data and functions

class DataOutput {
public:
  DataOutput(ParameterInput *pin);
  ~DataOutput();

  void CheckForOutputs(Mesh *pm);
//  void ComputeOutputData();
//  void WriteOutputData();

private:
  OutputBlock* pfirst_block_;  // pointer to first output block in linked list

};
#endif
