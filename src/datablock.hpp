#ifndef DATA_BLOCK_HPP
#define DATA_BLOCK_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file datablock.hpp
 *  \brief provides container class for passing data bewteen operator classes
 *         and outputs 
 *====================================================================================*/
class Mesh;

struct DataBlockHeader {

  Real time;
  Real dt;
  Real ncycle;
  std::string descriptor;
};

struct DataNodeHeader {

  std::string transform_types;
  std::string variable_types;
  int is, ie;
  int js, je;
  int ks, ke;
};

class DataNode {

public:
  DataNode(AthenaArray<Real> *pdinit, DataNodeHeader hinit);
  DataNode(AthenaArray<Real> *pdinit);
  DataNode(AthenaArray<Real> dinit, int index);
  ~DataNode();

  DataNode* GetNext() const { return pnext; }
  void SetNext( DataNode *psnext) { pnext = psnext; }
  AthenaArray<Real>* GetData() { return pdata; }
  DataNodeHeader GetHeader() {return header;}
  std::string GetVariableTypes() { return header.variable_types; }
  void SetVariableTypes(std::string vt) { header.variable_types = vt; }
  std::string GetTransformTypes() { return header.transform_types; }
  void AddTransformTypes(std::string tt) { header.transform_types += " " + tt; }
  void SetTransformTypes(std::string tt) { header.transform_types = tt; }

  void SetRanges(int is=0, int i =0, int js=0, int je=0, int ks=0, int ke=0);
  void GetRanges(int &is, int &ie, int &js, int &je, int &ks, int &ke);
  void GetRanges(int &is, int &ie, int &js, int &je);
  void GetRanges(int &is, int &ie);
  void GetDimensions(int &nx3, int &nx2, int &nx1);
 
private:
  DataNodeHeader header;  // Description of node
  AthenaArray<Real> *pdata;  // Data in this node
  DataNode *pnext;  // Next node

};

class DataBlock {

public:
  DataBlock();
  ~DataBlock();

  DataBlockHeader GetHeader() {return header;}
  DataNode* FirstNode() {return pfirst_node;}
  DataNode* LastNode() {return plast_node;}
  std::string GetDescriptor() {return header.descriptor;}
  void ConcatinateDescriptor(std::string str) {header.descriptor += str;}
  Real GetTime() {return header.time; }
  void SetTime(Real time) {header.time = time;}
  Real GetTimeStep() {return header.dt; }
  void SetTimeStep(Real dt) {header.dt = dt;}
  int GetCycleNumber() {return header.ncycle; } // Not sure this is necessary
  void SetCycleNumber(int ncycle) {header.ncycle = ncycle; }
  void InsertNode(DataNode *pnew_block);
  void InsertNode(AthenaArray<Real> *pdata, DataNodeHeader header);
  void InsertNode(AthenaArray<Real> *pdata);
  void ReplaceNode(DataNode *pnew_block, int n);
  DataNode* GetNode(int n);

private:
  DataBlockHeader header;  // Description of block and common data
  DataNode *pfirst_node;  // Pointer to first node
  DataNode *plast_node;  // Pointer to last node for easy insertion
};

//======================================================================================
//  DataBlockTransform definitions

//! \class DataBlockTransform
//  \brief  abstract base class for modifying DataBlocks

class DataBlockTransform {
public:
  DataBlockTransform();
  virtual void Transform(DataBlock* pDB) = 0;

  DataBlockTransform* GetNext() const { return pnext; }
private:
  DataBlockTransform* pnext;
};

/* -------------------------------------------------------------------------------------
 * All DataBlockTransforms are derived classes that provie a Transform function
 * that have DataBlock* as both input and return type.
 */
//! \class SumOverAll
//  \brief perfrom sum on DataBlock
class SumOverAll : public DataBlockTransform   {
public:
  void Transform(DataBlock* pDB);
};

//! \class SumOverx3
//  \brief perfrom sum on DataBlock
class SumOverx3 : public DataBlockTransform   {
public:
  void Transform(DataBlock* pDB);
};

//! \class SumOverx2
//  \brief perfrom sum on DataBlock
class SumOverx2 : public DataBlockTransform   {
public:
  void Transform(DataBlock* pDB);
};

//! \class SumOverx1
//  \brief perfrom sum on DataBlock
class SumOverx1 : public DataBlockTransform   {
public:
  void Transform(DataBlock* pDB);
};

//! \class Slice1Dx1
//  \brief perfrom 1d slice on DataBlock
class Slice1Dx1 : public DataBlockTransform   {
public:
  Slice1Dx1(Mesh* pm, Real slcx2, Real slcx3);
  ~Slice1Dx1();

  void Transform(DataBlock* pDB);
private:
  Real x2, x3;
  int ix2, ix3;
};

#endif
