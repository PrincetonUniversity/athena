//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "parameter_input.hpp"

// C++ headers
#include <cstdlib>    // atoi(), atof(), NULL, size_t
#include <fstream>    // ifstream
#include <iostream>   // endl, ostream
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena headers
#include "athena.hpp"  // Real

//======================================================================================
/*! \file parameter_input.cpp
 *  \brief implementation of functions in class ParameterInput
 *
 * PURPOSE: Member functions of this class are used to read and parse the input file.
 *   Functionality is loosely modeled after FORTRAN namelist.
 *
 * EXAMPLE of input file in 'Athena++' format:
 *
 *   <blockname1>      # block name; must be on a line by itself
 *                     # everything after a hash symbol is a comment and is ignored
 *   name1=value       # each parameter name must be on a line by itself
 *   name2 = value1    # whitespace around the = is optional
 *                     # blank lines are OK
 *   # my comment here   comment lines are OK
 *   # name3 = value3    values (and blocks) that are commented out are ignored
 *    
 *   <blockname2>      # start new block
 *   name1 = value1    # note that same parameter names can appear in different blocks
 *   name2 = value2    # empty lines (like following) are OK
 * 
 *   <blockname1>      # same blockname can re-appear, although NOT recommended
 *   name3 = value3    # this would be the 3rd parameter name in blockname1
 *   name1 = value4    # if parameter name is repeated, previous value is overwritten!
 *
 * LIMITATIONS:
 *   - parameter specification (name=val #comment) must all be on a single line
 *
 * HISTORY:
 *   - Nov 2002:  Created for Athena1.0/Cambridge release by Peter Teuben
 *   - 2003-2008: Many improvements and extensions by T. Gardiner and J.M. Stone
 *   - Jan 2014:  Rewritten in C++ for the Athena++ code by J.M. Stone
 *====================================================================================*/

// constructors

ParameterInput::ParameterInput()
{
  pfirst_block = NULL;
  last_filename_ = "";
}

// destructor - iterates through linked lists of blocks/lines and deletes each node

ParameterInput::~ParameterInput()
{
  InputBlock *pib=pfirst_block;
  while (pib != NULL) {
    InputBlock *pold_block = pib;
    pib = pib->pnext;
    delete pold_block;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn  void ParameterInput::LoadFromFile(std::string filename)
 *  \brief opens/reads/closes an input file. */  
/*  Input block names are allocated and stored in a linked list of InputBlocks.  Within
 *  each InputBlock the names, values, and comments of each parameter are allocated and
 *  stored in a linked list of InputLines.
 */

void ParameterInput::LoadFromFile(std::string filename)
{
  std::ifstream input_file;
  std::string line, block_name, param_name, param_value, param_comment;
  std::size_t first_char,last_char;
  std::stringstream msg;
  InputBlock *pib;

  if (last_filename_ == filename) {
    msg << "### FATAL ERROR in function [ParameterInput::LoadFromFile]" << std::endl
        << "Input file '" << filename << "' has already been read";
    throw std::runtime_error(msg.str().c_str());
  }

  input_file.open(filename.c_str());
  if (!input_file.is_open()) {
    msg << "### FATAL ERROR in function [ParameterInput::LoadFromFile]" << std::endl
        << "Input file '" << filename << "' could not be opened";
    throw std::runtime_error(msg.str().c_str());
  }

  while (input_file.good()) {
    getline(input_file,line);
    if (line.empty()) continue;                             // skip blank line
    first_char = line.find_first_not_of(" ");               // skip white space
    if (first_char == std::string::npos) continue;          // line is all white space
    if (line.compare(first_char,1,"#") == 0) continue;      // skip comments
    if (line.compare(first_char,9,"<par_end>") == 0) break; // stop on <par_end>

    if (line.compare(first_char,1,"<") == 0) {              // a new block
      first_char++;
      last_char = (line.find_first_of(">",first_char));       
      block_name.assign(line,first_char,last_char-1);       // extract block name

      if (last_char == std::string::npos) {
        msg << "### FATAL ERROR in function [ParameterInput::LoadFromFile]" << std::endl
            << "Block name '" << block_name << "' in input file '" << filename
            << "' not properly ended";
        throw std::runtime_error(msg.str().c_str());
      }

      pib = FindOrAddBlock(block_name);  // find or add block to linked list

      if (pib == NULL) {
        msg << "### FATAL ERROR in function [ParameterInput::LoadFromFile]" << std::endl
            << "Block name '" << block_name << "' could not be found/added";
        throw std::runtime_error(msg.str().c_str());
      }
      continue;  // skip to next line if block name was found
    }

// if line does not contain a block name, it must contain a parameter value.  So parse
// line and add name/value/comment strings (if found) to current block name

    ParseLine(pib,line,param_name,param_value,param_comment);
    AddParameter(pib,param_name,param_value,param_comment);
  }

  input_file.close();
  last_filename_ = filename;

  return;
}

//--------------------------------------------------------------------------------------
/*! \fn InputBlock* ParameterInput::FindOrAddBlock(std::string name)
 *  \brief find or add specified InputBlock.  Returns pointer to block.  */

InputBlock* ParameterInput::FindOrAddBlock(std::string name)
{
  InputBlock *pib, *plast;
  plast = pfirst_block;
  pib = pfirst_block;

// Search linked list of InputBlocks to see if name exists, return if found.

  while (pib != NULL) {
    if (name.compare(pib->block_name) == 0) return pib;    
    plast = pib;
    pib = pib->pnext;
  }

// Create new block in list if not found above

  pib = new InputBlock;
  pib->block_name.assign(name);  // store the new block name
  pib->pline = NULL;             // Terminate the InputLine list
  pib->pnext = NULL;             // Terminate the InputBlock list

// if this is the first block in list, save pointer to it in class

  if (pfirst_block == NULL) {
     pfirst_block = pib;
  } else {
    plast->pnext = pib;      // link new node into list
  }

  return pib;
}

//--------------------------------------------------------------------------------------
/*! \fn void ParameterInput::ParseLine(InputBlock *pib, std::string line,
 *           std::string& name, std::string& value, std::string& comment)
 *  \brief parse "name = value # comment" format, return name/value/comment strings. 
 */

void ParameterInput::ParseLine(InputBlock *pib, std::string line,
     std::string& name, std::string& value, std::string& comment)
{
  std::size_t first_char,last_char,equal_char,hash_char,len;

  first_char = line.find_first_not_of(" ");   // find first non-white space
  equal_char = line.find_first_of("=");       // find "=" char
  hash_char  = line.find_first_of("#");       // find "#" (optional)

// copy substring into name, remove white space at end of name

  len = equal_char - first_char;
  name.assign(line,first_char,len);

  last_char = name.find_last_not_of(" ");
  name.erase(last_char+1,std::string::npos);

// copy substring into value, remove white space at start and end

  len = hash_char - equal_char - 1;
  value.assign(line,equal_char+1,len);      

  first_char = value.find_first_not_of(" ");
  value.erase(0,first_char);

  last_char = value.find_last_not_of(" ");
  value.erase(last_char+1,std::string::npos);

// copy substring into comment, if present

  if (hash_char != std::string::npos) {
    comment = line.substr(hash_char);
  } else {
    comment = "";
  }
}

//--------------------------------------------------------------------------------------
/*! \fn static void AddParameter()
 *  \brief add name/value/comment tuple to the InputLine linked list in block *pb.  
 *
 *  If a parameter with the same name already exists, the value and comment strings
 *  are replaced (overwritten).
 */

void ParameterInput::AddParameter(InputBlock *pb, std::string name, 
     std::string value, std::string comment)
{
  InputLine *pl, *plast;

// Search linked list of InputLines to see if name exists.  This also sets *plast
// to point to last member of list

  pl = pb->pline;
  plast = pb->pline;
  while (pl != NULL) {
    if (name.compare(pl->param_name) == 0) {   // param name already exists
      pl->param_value.assign(value);           // replace existing param value
      pl->param_comment.assign(comment);       // replace exisiting param comment
      if(value.length() > pb->max_len_parvalue) pb->max_len_parvalue = value.length();
      return;
    }
    plast = pl;
    pl = pl->pnext;
  }

// Create new node in linked list if name does not already exist

  pl = new InputLine;
  pl->param_name.assign(name);
  pl->param_value.assign(value);
  pl->param_comment.assign(comment);
  pl->pnext = NULL;

// if this is the first parameter in list, save pointer to it in block.

  if (pb->pline == NULL) {
    pb->pline = pl;
    pb->max_len_parname = name.length();
    pb->max_len_parvalue = value.length();
  } else {
    plast->pnext = pl;  // link new node into list
    if(name.length() > pb->max_len_parname) pb->max_len_parname = name.length();
    if(value.length() > pb->max_len_parvalue) pb->max_len_parvalue = value.length();
  }

  return;
}

//--------------------------------------------------------------------------------------
/*! void ParameterInput::ModifyFromCmdline(int argc, char *argv[])
 *  \brief parse commandline for changes to input parameters */
// Note this function is very forgiving (no warnings!) if there is an error in format */

void ParameterInput::ModifyFromCmdline(int argc, char *argv[])
{
  std::string input_text,block,name,value;
  std::size_t slash_char,equal_char;
  std::stringstream msg;
  InputBlock *pb;
  InputLine *pl;

  for (int i=1; i<argc; i++) {
    input_text = argv[i];
    slash_char = input_text.find_first_of("/");   // find "/" character
    equal_char = input_text.find_first_of("=");   // find "=" character
// skip if either "/" or "=" do not exist in input
    if ((slash_char==std::string::npos) || (equal_char==std::string::npos)) continue;
    
// extract block/name/value strings

    block = input_text.substr(0,slash_char);
    name  = input_text.substr(slash_char+1,(equal_char - slash_char - 1));
    value = input_text.substr(equal_char+1,std::string::npos);

// get pointer to node with same block name in linked list of InputBlocks

    pb = GetPtrToBlock(block);
    if (pb == NULL) {
      msg << "### FATAL ERROR in function [ParameterInput::ModifyFromCmdline]" 
          << std::endl << "Block name '" << block << "' on command line not found";
      throw std::runtime_error(msg.str().c_str());
    }

// get pointer to node with same parameter name in linked list of InputLines

    pl = pb->GetPtrToLine(name);
    if (pl == NULL) {
      msg << "### FATAL ERROR in function [ParameterInput::ModifyFromCmdline]"
          << std::endl << "Parameter '" << name << "' in block '" << block 
          << "' on command line not found";
      throw std::runtime_error(msg.str().c_str());
    }
    pl->param_value.assign(value);   // replace existing value
    if(value.length() > pb->max_len_parvalue) pb->max_len_parvalue = value.length();
     
  }
}

//--------------------------------------------------------------------------------------
/*! \fn InputBlock* ParameterInput::GetPtrToBlock(std::string name)
 *  \brief return pointer to specified InputBlock if it exists */

InputBlock* ParameterInput::GetPtrToBlock(std::string name)
{
  InputBlock *pb;
  for (pb = pfirst_block; pb != NULL; pb = pb->pnext){
    if (name.compare(pb->block_name) == 0) return pb;    
  }
  return NULL;
}

//--------------------------------------------------------------------------------------
/*! \fn int ParameterInput::DoesParameterExist()
 *  \brief check whether parameter of given name in given block exists */

int ParameterInput::DoesParameterExist(std::string block, std::string name)
{
  InputLine *pl;
  InputBlock *pb;
  pb = GetPtrToBlock(block);
  if (pb == NULL) return 0;
  pl = pb->GetPtrToLine(name);
  return (pl == NULL ? 0 : 1);
}

//--------------------------------------------------------------------------------------
/*! \fn 
 *  \brief */

int ParameterInput::GetInteger(std::string block, std::string name)
{
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

// get pointer to node with same block name in linked list of InputBlocks

  pb = GetPtrToBlock(block);
  if (pb == NULL) {
    msg << "### FATAL ERROR in function [ParameterInput::GetInteger]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    throw std::runtime_error(msg.str().c_str());
  }

// get pointer to node with same parameter name in linked list of InputLines

  pl = pb->GetPtrToLine(name);
  if (pl == NULL) {
    msg << "### FATAL ERROR in function [ParameterInput::GetInteger]" << std::endl 
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    throw std::runtime_error(msg.str().c_str());
  }

// Convert string to integer and return value

  return atoi(pl->param_value.c_str());
}

//--------------------------------------------------------------------------------------
/*! \fn 
 *  \brief */

Real ParameterInput::GetReal(std::string block, std::string name)
{
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

// get pointer to node with same block name in linked list of InputBlocks

  pb = GetPtrToBlock(block);
  if (pb == NULL) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    throw std::runtime_error(msg.str().c_str());
  }

// get pointer to node with same parameter name in linked list of InputLines

  pl = pb->GetPtrToLine(name);
  if (pl == NULL) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    throw std::runtime_error(msg.str().c_str());
  }

// Convert string to real and return value

  return (Real)atof(pl->param_value.c_str());
}

//--------------------------------------------------------------------------------------
/*! \fn 
 *  \brief */

std::string ParameterInput::GetString(std::string block, std::string name)
{
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

// get pointer to node with same block name in linked list of InputBlocks

  pb = GetPtrToBlock(block);
  if (pb == NULL) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    throw std::runtime_error(msg.str().c_str());
  }

// get pointer to node with same parameter name in linked list of InputLines

  pl = pb->GetPtrToLine(name);
  if (pl == NULL) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    throw std::runtime_error(msg.str().c_str());
  }

// return value

  return pl->param_value;
}

//--------------------------------------------------------------------------------------
/*! \fn 
 *  \brief */

int ParameterInput::GetOrAddInteger(std::string block, std::string name, int value)
{
  InputBlock* pb;
  std::stringstream ss_value;

  if (DoesParameterExist(block, name)) return GetInteger(block,name);
  pb = FindOrAddBlock(block);
  ss_value << value;
  AddParameter(pb, name, ss_value.str(), "# Default value added at run time");
  return value;
}

//--------------------------------------------------------------------------------------
/*! \fn 
 *  \brief */

Real ParameterInput::GetOrAddReal(std::string block, std::string name, Real value)
{
  InputBlock* pb;
  std::stringstream ss_value;

  if (DoesParameterExist(block, name)) return GetReal(block,name);
  pb = FindOrAddBlock(block);
  ss_value << value;
  AddParameter(pb, name, ss_value.str(), "# Default value added at run time");
  return value;
}

//--------------------------------------------------------------------------------------
/*! \fn 
 *  \brief */

std::string ParameterInput::GetOrAddString(std::string block, std::string name, 
  std::string value)
{
  InputBlock* pb;
  std::stringstream ss_value;

  if (DoesParameterExist(block, name)) return GetString(block,name);
  pb = FindOrAddBlock(block);
  AddParameter(pb, name, value, "# Default value added at run time");
  return value;
}

//--------------------------------------------------------------------------------------
/*! \fn void ParameterInput::ParameterDump(std::ostream& os)
 *  \brief output entire InputBlock/InputLine hierarchy to specified stream */

void ParameterInput::ParameterDump(std::ostream& os)
{
  InputBlock *pb;
  InputLine *pl;
  std::string param_name,param_value;
  std::size_t len;

  os<< "#------------------------- PAR_DUMP -------------------------" << std::endl;

  for (pb = pfirst_block; pb != NULL; pb = pb->pnext){ // loop over InputBlocks
    os<< "<" << pb->block_name << ">" << std::endl;     // write block name
    for (pl = pb->pline; pl != NULL; pl = pl->pnext){   // loop over InputLines
      param_name.assign(pl->param_name);
      param_value.assign(pl->param_value);

      len = pb->max_len_parname - param_name.length() + 1;
      param_name.append(len,' ');                      // pad name to align vertically
      len = pb->max_len_parvalue - param_value.length() + 1;
      param_value.append(len,' ');                     // pad value to align vertically

      os<< param_name << "= " << param_value << pl->param_comment <<  std::endl;
    }
  }
  
  os<< "#------------------------- PAR_DUMP -------------------------" << std::endl;
  os<< "<par_end>" << std::endl;    // finish with par-end (useful in restart files)
}

//--------------------------------------------------------------------------------------
// InputBlock constructor and destructor

InputBlock::InputBlock()
{
}

// destructor - iterates through linked lists of blocks/lines and deletes each node

InputBlock::~InputBlock()
{
}

//--------------------------------------------------------------------------------------
/*! \fn InputLine* InputBlock::GetPtrToLine(std::string name)
 *  \brief return pointer to InputLine containing specified parameter if it exists */

InputLine* InputBlock::GetPtrToLine(std::string name)
{
  for(InputLine* pl = pline; pl != NULL; pl = pl->pnext){
    if (name.compare(pl->param_name) == 0) return pl;    
  }
  return NULL;
}
