#ifndef ATHENA_ARRAYS_HPP
#define ATHENA_ARRAYS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file athena_arrays.hpp
 *  \brief provides array classes valid in 1D to 4D.
 *
 *  The operator() is overloaded, e.g. elements of a 4D array of size [N4xN3xN2xN1]
 *  are accessed as:  A(n,k,j,i) = A[i + N1*(j + N2*(k + N3*n))]
 *  NOTE THE TRAILING INDEX INSIDE THE PARENTHESES IS INDEXED FASTEST
 *====================================================================================*/

template<typename T>
class AthenaArray {
public:
  AthenaArray();
  ~AthenaArray();

// public functions to allocate/deallocate memory for 1D/2D/3D/4D data

  void NewAthenaArray(int nx1);
  void NewAthenaArray(int nx2, int nx1);
  void NewAthenaArray(int nx3, int nx2, int nx1);
  void NewAthenaArray(int nx4, int nx3, int nx2, int nx1);
  void DeleteAthenaArray();

// functions to get array dimensions 

  int GetDim1() const { return nx1_; }
  int GetDim2() const { return nx2_; }
  int GetDim3() const { return nx3_; }
  int GetDim4() const { return nx4_; }

  bool IsShallowCopy() { return (scopy_ == 1); }
  T *data() { return pdata_; }
  const T *data() const	{ return pdata_; }

// overload operator() to access 1d/2d/3d/4d data

  T &operator() (const int n) { 
    return pdata_[n]; }
  T operator() (const int n) const { 
    return pdata_[n]; }

  T &operator() (const int n, const int i) { 
    return pdata_[i + nx1_*n]; }
  T operator() (const int n, const int i) const { 
    return pdata_[i + nx1_*n]; }

  T &operator() (const int n, const int j, const int i) { 
    return pdata_[i + nx1_*(j + nx2_*n)]; }
  T operator() (const int n, const int j, const int i) const { 
    return pdata_[i + nx1_*(j + nx2_*n)]; }

  T &operator() (const int n, const int k, const int j, const int i) { 
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*n))]; }
  T operator() (const int n, const int k, const int j, const int i) const { 
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*n))]; }

// copy constructor and overloaded assignment operator (both do deep copies).
// A shallow copy function is also provided.

  AthenaArray(const AthenaArray<T>& t);
  AthenaArray<T> &operator= (const AthenaArray<T> &t);
  AthenaArray<T> ShallowCopy();
  AthenaArray<T>* ShallowCopy(const int n);

private:
  T *pdata_;
  int nx1_, nx2_, nx3_, nx4_;
  int scopy_;  // =0 if shallow copy (prevents source from being deleted)
};

//constructor

template<typename T>
AthenaArray<T>::AthenaArray()
  : pdata_(0), nx1_(0), nx2_(0), nx3_(0), nx4_(0), scopy_(0)
{
}

// destructor

template<typename T>
AthenaArray<T>::~AthenaArray()
{
  if (scopy_ == 0) DeleteAthenaArray();
}

// copy constructor (does a deep copy)

template<typename T>
AthenaArray<T>::AthenaArray(const AthenaArray<T>& src) {
  nx1_ = src.nx1_;
  nx2_ = src.nx2_;
  nx3_ = src.nx3_;
  nx4_ = src.nx4_;
  if (src.pdata_) {
    size_t size = (src.nx1_)*(src.nx2_)*(src.nx3_)*(src.nx4_);
    pdata_ = new T[size];
    for (size_t i=0; i<size; ++i) {
      pdata_[i] = src.pdata_[i];
    } 
  }
}

// assignment operator (does a deep copy)

template<typename T>
AthenaArray<T> &AthenaArray<T>::operator= (const AthenaArray<T> &src) {
  if (this != &src){
    this.nx1_ = src.nx1_;
    this.nx2_ = src.nx2_;
    this.nx3_ = src.nx3_;
    this.nx4_ = src.nx4_;

    delete[] this.pdata_;
    size_t size = (src.nx1_)*(src.nx2_)*(src.nx3_)*(src.nx4_);
    this.pdata_ = new T[size];
    for (size_t i=0; i<size; ++i) {
      this.pdata_[i] = src.pdata_[i];
    } 
  }
  return *this;
}

//--------------------------------------------------------------------------------------
//! \fn AthenaArray::ShallowCopy()
//  \brief shallow copy of array 

template<typename T>
AthenaArray<T> AthenaArray<T>::ShallowCopy() {
  AthenaArray<T> dest;
  dest.nx1_=nx1_;
  dest.nx2_=nx2_;
  dest.nx3_=nx3_;
  dest.nx4_=nx4_;
  dest.pdata_ = pdata_;
  dest.scopy_ = 1;
  return dest;
}

//--------------------------------------------------------------------------------------
//! \fn AthenaArray::ShallowCopy(int index)
//  \brief shallow copy of a subset of an array to a new dynamically allocated array

template<typename T>
AthenaArray<T>* AthenaArray<T>::ShallowCopy(const int n) {
  AthenaArray<T> *dest = new AthenaArray<T>;
  dest->nx1_=nx1_;
  dest->nx2_=nx2_;
  dest->nx3_=nx3_;
  dest->nx4_=1;
// No error checking: n must be < nx4 in original array
  dest->pdata_ = pdata_ + n * nx1_*nx2_*nx3_;
  dest->scopy_ = 1;
  return dest;
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief 1d data allocation

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx1)
{
  nx1_ = nx1;
  nx2_ = 1;
  nx3_ = 1;
  nx4_ = 1;
  pdata_ = new T[nx1];
}
 
//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief 2d data allocation */

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx2, int nx1)
{
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = 1;
  nx4_ = 1;
  pdata_ = new T[nx1*nx2];
}
 
//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief 3d data allocation */

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx3, int nx2, int nx1)
{
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = 1;
  pdata_ = new T[nx1*nx2*nx3];
}
 
//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief 4d data allocation */

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx4, int nx3, int nx2, int nx1)
{
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  pdata_ = new T[nx1*nx2*nx3*nx4];
}

//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief  free memory allocated for data array */

template<typename T>
void AthenaArray<T>::DeleteAthenaArray()
{
  delete[] pdata_;
} 
#endif
