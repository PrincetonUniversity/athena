#ifndef OUTPUT_HPP
#define OUTPUT_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \fileData output_example.hpp
 *  \brief provides class to handle all types of data output
 *====================================================================================*/

class ParameterInput;
class Mesh;
class DataBlock;

//! \struct OutputBlock
//  \brief  node in a linked list of output blocks in the input file


class  OutputBlock {
public:
  OutputBlock(InputBlock *pin_block);
  ~OutputBlock();

  int block_number;
  int number;
  std::string block_name;
  std::string filebase;
  std::string out;
  std::string dat_fmt;
  std::string trans;
  Real last_time, dt;

};

//! \class Output
//  \brief  abstract base class for output node of linked list

class Output {

public:
  Output(InputBlock *pin_block, Mesh *pm);
  ~Output() {};

// Write() is the only pure virtual function every derived Output must know
// how to write itself
  virtual void Write() = 0;

  virtual void ComputeFromMesh(Mesh *pm) {}; // should not be pure virtual
// Functions for providing protected data access
  Output* GetNext() const { return pnext; }
  void SetNext(Output *pset) { pnext = pset; }
  OutputBlock GetOutputBlock() {return output_block; }
  Real NextTime() { return output_block.dt + output_block.last_time; }
  void Update() { output_block.last_time += output_block.dt; ++output_block.number; }
  int BlockNumber() { return output_block.block_number; }
  int Number() { return output_block.number; }
  std::string BlockName() { return output_block.block_name; }
  std::string TransformList() { return output_block.trans; }

protected:
  OutputBlock output_block;
  DataBlock *pdata;
  Output *pnext;
  DataBlockTransform *ptrans;

  DataBlock* LoadFluidFromMesh(Mesh *pm);
};

//! \class OutputList
//  \brief  list of output nodes

class OutputList {

public:
  OutputList(ParameterInput *pin, Mesh *pm);
  ~OutputList();

  void CheckForOutputs(Mesh *pm);

private:
  Output *phead;
};

//! \class VTKOutput
//  \brief  output class for VTK
/*
class VTKOutput : public Output {

public:
  VTKOutput(InputBlock *pin_block);
  ~VTKOutput() {};

  void Write();

};
*/

//! \class FormattedTableOutput
//  \brief  output class for tabular data

class FormattedTableOutput : public Output {

public:
  FormattedTableOutput(InputBlock *pin_block, Mesh *pm);
  ~FormattedTableOutput() {};

  void ComputeFromMesh(Mesh *pM);
  void Write();

private:
  int nx1, nx2, nx3;

  void WriteNodeData(DataNode *pdn, FILE *pfile, std::string fmt);
};

//! \class HistoryOutput
//  \brief  output class for history dumps

class HistoryOutput : public Output {

public:
  HistoryOutput(InputBlock *pin_block, Mesh *pm);
  ~HistoryOutput() {};

  void ComputeFromMesh(Mesh *pM);
  void Write();

  };


#endif
