#ifndef ATHENA_ARRAYS_HPP
#define ATHENA_ARRAYS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file athena_arrays.hpp
 *  \brief provides array objects
 *
 *  Functionally, all arrays are allocated as 1D vectors.  The () operator is overloaded
 *  so, e.g. a 4D array data elements can be accessed as 
 *      A(n,k,j,i) = A[i + nx1*(j + nx2*(k + nx3*n))]
 *  NOTE THE TRAILING INDEX INSIDE THE PARENTHESES IS INDEXED FASTEST, AS IN C
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

  T* data() { return data_; }
  const T* data() const	{ return data_; }

// overload () operator to access 1d/2d/3d/4d data

  T& operator() (const int n) { 
    return data_[n]; }
  T operator() (const int n) const { 
    return data_[n]; }

  T& operator() (const int n, const int i) { 
    return data_[i + nx1_*n]; }
  T operator() (const int n, const int i) const { 
    return data_[i + nx1_*n]; }

  T& operator() (const int n, const int j, const int i) { 
    return data_[i + nx1_*(j + nx2_*n)]; }
  T operator() (const int n, const int j, const int i) const { 
    return data_[i + nx1_*(j + nx2_*n)]; }

  T& operator() (const int n, const int k, const int j, const int i) { 
    return data_[i + nx1_*(j + nx2_*(k + nx3_*n))]; }
  T operator() (const int n, const int k, const int j, const int i) const { 
    return data_[i + nx1_*(j + nx2_*(k + nx3_*n))]; }

private:

  T*  data_;
  int nx1_, nx2_, nx3_, nx4_;
};

//constructor

template<typename T>
AthenaArray<T>::AthenaArray()
  : data_(0), nx1_(0), nx2_(0), nx3_(0), nx4_(0)
{
}

// destructor

template<typename T>
AthenaArray<T>::~AthenaArray()
{
  DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief 1d data allocation */

template<typename T>
void AthenaArray<T>::NewAthenaArray(int nx1)
{
  nx1_ = nx1;
  nx2_ = 1;
  nx3_ = 1;
  nx4_ = 1;
  data_ = new T[nx1];
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
  data_ = new T[nx1*nx2];
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
  data_ = new T[nx1*nx2*nx3];
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
  data_ = new T[nx1*nx2*nx3*nx4];
}

//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief  free memory allocated for data array */

template<typename T>
void AthenaArray<T>::DeleteAthenaArray()
{
  delete[] data_;
} 
#endif
