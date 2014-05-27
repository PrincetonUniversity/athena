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

#include <sstream>
#include <string>
#include <stdexcept>
//#include <math.h>
//#include <float.h>
//#include <iostream>

#include "athena.hpp"
#include "athena_arrays.hpp"
#include "datablock.hpp"
#include "mesh.hpp"

//======================================================================================
/*! \file datablock.cpp
 *  \brief implements container class for passing data to/from input/output
 *====================================================================================*/


//--------------------------------------------------------------------------------------
// DataNode constructor

DataNode::DataNode(AthenaArray<Real> *pdinit, DataNodeHeader hinit)
{
  pdata = pdinit;
  header = hinit;
  pnext = NULL;
}

DataNode::DataNode(AthenaArray<Real> *pdinit)
{
  pdata = pdinit;
  header.transform_types = " ";
  header.variable_types = " ";
  SetRanges();
  pnext = NULL;
}


// DataNode destructor
DataNode::~DataNode()
{

// Delete array pointed to by pdata
  if (!pdata->IsShallowCopy()) pdata->DeleteAthenaArray();

}

//--------------------------------------------------------------------------------------
/*! \fn void DataNode::SetRanges(int is, int ie, int js, int je, int ks, int ke)
 *  \brief Set output ranges for DataNode
 *  Note that each variable is set to 0 by default
 */
void DataNode::SetRanges(int is, int ie, int js, int je, int ks, int ke)
{
  header.is = is;
  header.ie = ie;
  header.js = js;
  header.je = je;
  header.ks = ks;
  header.ke = ke;
}

//--------------------------------------------------------------------------------------
/*! \fn void DataNode::GetRanges(int &is, int &ie, int &js, int &je, int &ks, int &ke)
 *  \brief Get output ranges from DataNode
 */
void DataNode::GetRanges(int &is, int &ie, int &js, int &je, int &ks, int &ke)
{
  is = header.is;
  ie = header.ie;
  js = header.js;
  je = header.je;
  ks = header.ks;
  ke = header.ke;
}

//--------------------------------------------------------------------------------------
/*! \fn void DataNode::GetRanges(int &is, int &ie, int &js, int &je)
 *  \brief Get output ranges from DataNode
 */
void DataNode::GetRanges(int &is, int &ie, int &js, int &je)
{
  is = header.is;
  ie = header.ie;
  js = header.js;
  je = header.je;
}

//--------------------------------------------------------------------------------------
/*! \fn void DataNode::GetRanges(int &is, int &ie)
 *  \brief Get output ranges from DataNode
 */
void DataNode::GetRanges(int &is, int &ie)
{
  is = header.is;
  ie = header.ie;
}

//--------------------------------------------------------------------------------------
/*! \fn void DataNode::GetDimensions(int &nx3, int &nx2, int &nx1)
 *  \brief Get output dimensions from DataNode
 */
void  DataNode::GetDimensions(int &nx3, int &nx2, int &nx1)
{
  nx3 = header.ke-header.ks+1;
  nx2 = header.je-header.js+1;
  nx1 = header.ie-header.is+1;
}

//--------------------------------------------------------------------------------------
/*! \fn void DataNode::Where(Real x)
 *  \brief Finx position nearest value x in 1D array
 */
/*void  DataNode::Where(int x)

{
  AthenaArray* pA=GetData();
  int is, ie;
  GetRanges(is,ie);
  for (int i=is; i <= ie; ++i) {

  }
  }*/

//--------------------------------------------------------------------------------------
// DataBlock constructor
DataBlock::DataBlock()
{
  pfirst_node = NULL;
  plast_node = NULL;
  header.descriptor = "";
}

// DataBlock destructor
DataBlock::~DataBlock()
{
  DataNode *pnext_node, *pcurrent_node = pfirst_node;

// delete DataBlock nodes, starting with head node
  pnext_node = pcurrent_node->GetNext();
  while(pnext_node != NULL) {
    delete pcurrent_node;
    pcurrent_node = pnext_node;
    pnext_node = pcurrent_node->GetNext();
  }
  delete pcurrent_node;

}

//--------------------------------------------------------------------------------------
/*! \fn void DataBlock::InsertNode(DataNode *pnew_node)
 *  \brief Add new node to DataBlock
 */
void DataBlock::InsertNode(DataNode *pnew_node)
{

  if (pfirst_node == NULL) 
    pfirst_node = plast_node = pnew_node;
  else {
    plast_node->SetNext(pnew_node);
    plast_node = pnew_node;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void DataBlock::InsertNode(AthenaArray<Real> *data, std::string header)
 *  \brief Add new node to DataBlock
 */
void DataBlock::InsertNode(AthenaArray<Real> *pdata, DataNodeHeader header)
{
  DataNode *pnew_node = new DataNode(pdata,header);
  InsertNode(pnew_node);

}

//--------------------------------------------------------------------------------------
/*! \fn void DataBlock::InsertNode(AthenaArray<Real> *data)
 *  \brief Add new node to DataBlock
 */
void DataBlock::InsertNode(AthenaArray<Real> *pdata)
{
  DataNode *pnew_node = new DataNode(pdata);
  InsertNode(pnew_node);

}

//--------------------------------------------------------------------------------------
/*! \fn void DataBlock::InsertNode(DataNode *pnew_node)
 *  \brief Add new node to DataBlock at position n, deleting old node
 */
void  DataBlock::ReplaceNode(DataNode *pnew_node, int n)
{
  DataNode *pold, *pprev;

  if (n == 0) {
    pold = pfirst_node;
    pfirst_node = pnew_node;
  } else {
    pprev = GetNode(n-1);
    pold = pprev->GetNext();
    pprev->SetNext(pnew_node);
  }
  pnew_node->SetNext(pold->GetNext());
  if (pold == plast_node) plast_node = pnew_node;
  delete pold;
}

//--------------------------------------------------------------------------------------
/*! \fn DataNode* DataBlock::GetNode(int n)
 *  \brief Return pointer to nth DataNode
 */
DataNode* DataBlock::GetNode(int n)
{
  int i = 0;
  DataNode* pdn = FirstNode();
  while ( (pdn->GetNext() != NULL) && (i < n) ){
    pdn = pdn->GetNext();
    i++;
  }
  return pdn;  // returns last node if n exceeds # of nodes
}

//======================================================================================
//  derived DataBlockTransform implementations

//--------------------------------------------------------------------------------------
//  DataBlockTransform constructor
DataBlockTransform::DataBlockTransform()
{
  pnext = NULL;
}

void SumOverAll::Transform(DataBlock* pdb)
{
  DataNode* pdn;
  AthenaArray<Real> *pA;
  AthenaArray<Real> *pSum;
  DataNode* pdn_new;
  DataNodeHeader head;
// Assumes first 3 node contain mesh coordinates
  int node_num = 3;
  pdn = pdb->GetNode(node_num);
  while (pdn != NULL) {
    pA = pdn->GetData();
    pSum = new AthenaArray<Real>;
    pSum->NewAthenaArray(pA->GetDim4()); // pSum may have more than 1 variable

// Loop over variables and dimensions
    for (int n=0; n<pA->GetDim4(); ++n){
      (*pSum)(n) = 0.0;  //Initalize to zero
      for (int k=0; k<pA->GetDim3(); ++k){
        for (int j=0; j<pA->GetDim2(); ++j){
          for (int i=0; i<pA->GetDim1(); ++i){
            (*pSum)(n) += (*pA)(n,k,j,i);
          }}}
    }
// Insert Sum into new DataNode
    head = pdn->GetHeader();
    pdn_new= new DataNode(pSum,head);
    pdn_new->AddTransformTypes("SumOverAll");
    pdb->ReplaceNode(pdn_new,node_num++);
    pdn = pdn->GetNext();
  }
}

void SumOverx3::Transform(DataBlock* pdb)
{

  DataNode* pdn = pdb->GetNode(3);
  AthenaArray<Real> *pA = pdn->GetData();
  AthenaArray<Real> *pSum = new AthenaArray<Real>;

  pSum->NewAthenaArray(pA->GetDim4(),1,pA->GetDim2(),pA->GetDim1());

// Loop over variables and dimensions
  for (int n=0; n<pA->GetDim4(); ++n){
    (*pSum)(n) = 0.0;  //Initalize to zero
    for (int k=0; k<pA->GetDim3(); ++k){
      for (int j=0; j<pA->GetDim2(); ++j){
        for (int i=0; i<pA->GetDim1(); ++i){
          (*pSum)(n,0,j,i) += (*pA)(n,k,j,i);
        }}}
  }
// Insert Sum into new DataNode
//  DataNode* pdbn = new DataNode(pSum,pdb->GetHeader());
  // pdbn->AddTransformTypes("SumOverx3");
  //return pdbn;
}

void SumOverx2::Transform(DataBlock* pdb)
{

  DataNode* pdn = pdb->GetNode(3);
  AthenaArray<Real> *pA = pdn->GetData();
  AthenaArray<Real> *pSum = new AthenaArray<Real>;

  pSum->NewAthenaArray(pA->GetDim4(),pA->GetDim3(),1,pA->GetDim1());

// Loop over variables and dimensions
  for (int n=0; n<pA->GetDim4(); ++n){
    (*pSum)(n) = 0.0;  //Initalize to zero
    for (int k=0; k<pA->GetDim3(); ++k){
      for (int j=0; j<pA->GetDim2(); ++j){
        for (int i=0; i<pA->GetDim1(); ++i){
          (*pSum)(n,k,0,i) += (*pA)(n,k,j,i);
        }}}
  }
// Insert Sum into new DataNode
//  DataNode* pdbn = new DataNode(pSum,pdb->GetHeader());
//  pdbn->AddTransformTypes("SumOverX\x2");
//  return pdbn;
}

void SumOverx1::Transform(DataBlock* pdb)
{

  DataNode* pdn = pdb->GetNode(3);
  AthenaArray<Real> *pA = pdn->GetData();
  AthenaArray<Real> *pSum = new AthenaArray<Real>;

  pSum->NewAthenaArray(pA->GetDim4(),pA->GetDim3(),pA->GetDim2(),1);

// Loop over variables and dimensions
  for (int n=0; n<pA->GetDim4(); ++n){
    (*pSum)(n) = 0.0;  //Initalize to zero
    for (int k=0; k<pA->GetDim3(); ++k){
      for (int j=0; j<pA->GetDim2(); ++j){
        for (int i=0; i<pA->GetDim1(); ++i){
          (*pSum)(n,k,j,0) += (*pA)(n,k,j,i);
        }}}
  }
// Insert Sum into new DataNode
  // DataNode* pdbn = new DataNode(pSum,pdb->GetHeader());
  // pdbn->AddTransformTypes("SumOverx1");
  //return pdbn;
}

// -------------------------------------------------------------------------------------
//  Slice1Dx1 constructor

Slice1Dx1::Slice1Dx1(Mesh* pm, Real slcx2, Real slcx3)
{
  x2 = slcx2;
  x3 = slcx3;

  Block *pblock = pm->pdomain->pblock;
  int i;
  for (i=pblock->js; i<= pblock->je; ++i) {
    ix2 = i;
    if (slcx2 > pm->pdomain->pblock->x2v(i)) break;
  }
  for (i=pblock->ks; i<= pblock->ke; ++i) {
    ix3 = i;
    if (slcx2 > pm->pdomain->pblock->x2v(i)) break;
  }
}

Slice1Dx1::~Slice1Dx1()
{

}

void Slice1Dx1::Transform(DataBlock* pdb)
{
// Modify descriptor in header
  Real ox2 = (*pdb->GetNode(1)->GetData())(ix2);
  Real ox3 =  (*pdb->GetNode(1)->GetData())(ix3);
  std::stringstream str;
  str << "1D Slice along x1 at x2= " << ox2 << " and x3= " << ox3;
  pdb->ConcatinateDescriptor(str.str());
// Slice arrays containing output data  
  DataNode* pdn;
  AthenaArray<Real> *pA;
  AthenaArray<Real> *pSlice;
  DataNode* pdn_new;
  DataNodeHeader head;
// Assumes first 3 node contain mesh coordinates
  int node_num = 3;
  pdn = pdb->GetNode(node_num);
  while (pdn != NULL) {
    pA = pdn->GetData();
    pSlice = new AthenaArray<Real>;
    pSlice->NewAthenaArray(pA->GetDim4(),1,1,pA->GetDim1());

// Loop over variables and dimensions
    for (int n=0; n<pA->GetDim4(); ++n){
      for (int i=0; i<pA->GetDim1(); ++i){
        (*pSlice)(n,0,0,i) = (*pA)(n,ix3,ix2,i);
      }}
// Insert Slice into new DataNode
    head = pdn->GetHeader();
    head.js = head.je = head.ks = head.ke = 0;
    pdn_new= new DataNode(pSlice,head);
    pdn_new->AddTransformTypes("Slice1Dx1");
    pdb->ReplaceNode(pdn_new,node_num++);
    pdn = pdn->GetNext();
  }
}
