#ifndef PARAMETER_INPUT_HPP
#define PARAMETER_INPUT_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file parameter_input.hpp
 *  \brief definition of class ParameterInput
 *  Contains data structures used to store, and functions used to access, parameters
 *  read from the input file.  See comments at start of parameter_input.cpp for more
 *  information on the Athena++ input file format.
 *====================================================================================*/

// C++ headers
#include <cstddef>  // size_t
#include <ostream>  // ostream
#include <string>   // string

// Athena headers
#include "athena.hpp"  // Real

//! \struct InputLine
//  \brief  node in a linked list of parameters contained within a single input block

typedef struct InputLine {   
  std::string param_name;
  std::string param_value;    // value of the parameter is stored as a string!
  std::string param_comment;
  struct InputLine *pnext;    // pointer to the next node
} InputLine;

//! \class InputBlock
//  \brief  node in a linked list of all input blocks contained within input file

class InputBlock { 
public:
  InputBlock();
  ~InputBlock();

  std::string block_name;
  std::size_t max_len_parname;  // length of longest param_name, for nice-looking output
  std::size_t max_len_parvalue; // length of longest param_value, to format outputs
  InputLine *pline;             // pointer to first InputLine in this block
  InputBlock *pnext;            // pointer to the next node

  InputLine* GetPtrToLine(std::string name);
};

//! \class ParameterInput
//  \brief data and definitions of functions used to store and access input parameters
//  Functions are implemented in parameter_input.cpp

class ParameterInput {
public:
  ParameterInput();
  ~ParameterInput();

  void LoadFromFile(std::string filename);
  void ModifyFromCmdline(int argc, char *argv[]);
  void ParameterDump(std::ostream& os);
  int  DoesParameterExist(std::string block, std::string name);
  int  GetInteger(std::string block, std::string name);
  int  GetOrAddInteger(std::string block, std::string name, int value);
  Real GetReal(std::string block, std::string name);
  Real GetOrAddReal(std::string block, std::string name, Real value);
  std::string GetString(std::string block, std::string name);
  std::string GetOrAddString(std::string block, std::string name, std::string value);

  InputBlock* pfirst_block;   // pointer to first input block in linked list

private:
  std::string last_filename_;  // last input file opened, to prevent duplicate reads

  void ParseLine(InputBlock *pb, std::string line, std::string& name,
       std::string& value, std::string& comment);
  void AddParameter(InputBlock *pb, std::string name, std::string value,
       std::string comment);
  InputBlock* FindOrAddBlock(std::string name);
  InputBlock* GetPtrToBlock(std::string name);
};
#endif
