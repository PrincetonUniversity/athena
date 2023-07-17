#ifndef ATHENA_ARRAYS_HPP_
#define ATHENA_ARRAYS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_arrays.hpp
//! \brief provides array classes valid in 1D to 6D.
//!
//! The operator() is overloaded, e.g. elements of a 4D array of size [N4xN3xN2xN1]
//! are accessed as:  A(n,k,j,i) = A[i + N1*(j + N2*(k + N3*n))]
//!
//! **NOTE THE TRAILING INDEX INSIDE THE PARENTHESES IS INDEXED FASTEST**

// C headers

// C++ headers
#include <cstddef>  // size_t
#include <cstring>  // memset()
#include <utility>  // swap()

// Athena++ headers

template <typename T>
class AthenaArray {
 public:
  enum class DataStatus {empty, shallow_slice, allocated};  // formerly, "bool scopy_"
  // ctors
  // default ctor: simply set null AthenaArray
  AthenaArray() : pdata_(nullptr), nx1_(0), nx2_(0), nx3_(0),
                  nx4_(0), nx5_(0), nx6_(0), state_(DataStatus::empty) {}
  // ctor overloads: set expected size of unallocated container, maybe allocate (default)
  explicit AthenaArray(int nx1, DataStatus init=DataStatus::allocated) :
      pdata_(nullptr), nx1_(nx1), nx2_(1), nx3_(1), nx4_(1), nx5_(1), nx6_(1),
      state_(init) { AllocateData(); }
  AthenaArray(int nx2, int nx1, DataStatus init=DataStatus::allocated) :
      pdata_(nullptr), nx1_(nx1), nx2_(nx2), nx3_(1), nx4_(1), nx5_(1), nx6_(1),
      state_(init) { AllocateData(); }
  AthenaArray(int nx3, int nx2, int nx1, DataStatus init=DataStatus::allocated) :
      pdata_(nullptr), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(1), nx5_(1), nx6_(1),
      state_(init) { AllocateData(); }
  AthenaArray(int nx4, int nx3, int nx2, int nx1, DataStatus init=DataStatus::allocated) :
      pdata_(nullptr), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(nx4), nx5_(1), nx6_(1),
      state_(init) { AllocateData(); }
  AthenaArray(int nx5, int nx4, int nx3, int nx2, int nx1,
              DataStatus init=DataStatus::allocated) :
      pdata_(nullptr), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(nx4), nx5_(nx5),  nx6_(1),
      state_(init) { AllocateData(); }
  AthenaArray(int nx6, int nx5, int nx4, int nx3, int nx2, int nx1,
              DataStatus init=DataStatus::allocated) :
      pdata_(nullptr), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(nx4), nx5_(nx5), nx6_(nx6),
      state_(init) { AllocateData(); }
  // still allowing delayed-initialization (after constructor) via array.NewAthenaArray()
  // or array.InitWithShallowSlice() (only used in outputs.cpp + 3x other files)
  //! \todo (felker):
  //! - replace InitWithShallowSlice with ??? and remove shallow_copy enum val
  //! - replace raw pointer with std::vector + reshape (if performance is same)

  // user-provided dtor, "rule of five" applies:
  ~AthenaArray();
  // define copy constructor and overload assignment operator so both do deep copies.
  __attribute__((nothrow)) AthenaArray(const AthenaArray<T>& t);
  __attribute__((nothrow)) AthenaArray<T> &operator= (const AthenaArray<T> &t);
  // define move constructor and overload assignment operator to transfer ownership
  __attribute__((nothrow)) AthenaArray(AthenaArray<T>&& t);
  __attribute__((nothrow)) AthenaArray<T> &operator= (AthenaArray<T> &&t);

  // public functions to allocate/deallocate memory for 1D-6D data
  __attribute__((nothrow)) void NewAthenaArray(int nx1);
  __attribute__((nothrow)) void NewAthenaArray(int nx2, int nx1);
  __attribute__((nothrow)) void NewAthenaArray(int nx3, int nx2, int nx1);
  __attribute__((nothrow)) void NewAthenaArray(int nx4, int nx3, int nx2, int nx1);
  __attribute__((nothrow)) void NewAthenaArray(int nx5, int nx4, int nx3, int nx2,
                                               int nx1);
  __attribute__((nothrow)) void NewAthenaArray(int nx6, int nx5, int nx4, int nx3,
                                               int nx2, int nx1);
  void DeleteAthenaArray();

  // public function to swap underlying data pointers of two equally-sized arrays
  void SwapAthenaArray(AthenaArray<T>& array2);
  // public function to exchange two arrays
  void ExchangeAthenaArray(AthenaArray<T>& array2);
  void ZeroClear();

  // functions to get array dimensions
  int GetDim1() const { return nx1_; }
  int GetDim2() const { return nx2_; }
  int GetDim3() const { return nx3_; }
  int GetDim4() const { return nx4_; }
  int GetDim5() const { return nx5_; }
  int GetDim6() const { return nx6_; }

  // a function to get the total size of the array
  int GetSize() const {
    if (state_ == DataStatus::empty)
      return 0;
    else
      return nx1_*nx2_*nx3_*nx4_*nx5_*nx6_;
  }
  std::size_t GetSizeInBytes() const {
    if (state_ == DataStatus::empty)
      return 0;
    else
      return nx1_*nx2_*nx3_*nx4_*nx5_*nx6_*sizeof(T);
  }

  bool IsShallowSlice() { return (state_ == DataStatus::shallow_slice); }
  bool IsEmpty() { return (state_ == DataStatus::empty); }
  bool IsAllocated() { return (state_ == DataStatus::allocated); }
  // "getter" function to access private data member
  //! \todo (felker):
  //! - Replace this unrestricted "getter" with a limited, safer alternative.
  //! - Rename function. Conflicts with "AthenaArray<> data" OutputData member.
  T *data() { return pdata_; }
  const T *data() const { return pdata_; }

  // overload "function call" operator() to access 1d-6d data
  // provides Fortran-like syntax for multidimensional arrays vs. "subscript" operator[]

  // "non-const variants" called for "AthenaArray<T>()" provide read/write access via
  // returning by reference, enabling assignment on returned l-value, e.g.: a(3) = 3.0;
  T &operator() (const int n) {
    return pdata_[n]; }
  // "const variants" called for "const AthenaArray<T>" returns T by value, since T is
  // typically a built-in type (versus "const T &" to avoid copying for general types)
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

  T &operator() (const int m, const int n, const int k, const int j, const int i) {
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*(n + nx4_*m)))]; }
  T operator() (const int m, const int n, const int k, const int j, const int i) const {
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*(n + nx4_*m)))]; }

  // int l?, int o?
  T &operator() (const int p, const int m, const int n, const int k, const int j,
                 const int i) {
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*(n + nx4_*(m + nx5_*p))))]; }
  T operator() (const int p, const int m, const int n, const int k, const int j,
                const int i) const {
    return pdata_[i + nx1_*(j + nx2_*(k + nx3_*(n + nx4_*(m + nx5_*p))))]; }

  // (deferred) initialize an array with slice from another array
  void InitWithShallowSlice(AthenaArray<T> &src, const int dim, const int indx,
                            const int nvar);

  void ShallowSlice3DToPencil(AthenaArray<T> &src, const int k, const int j,
                              const int il, const int n);

 private:
  T *pdata_;
  int nx1_, nx2_, nx3_, nx4_, nx5_, nx6_;
  DataStatus state_;  // describe what "pdata_" points to and ownership of allocated data

  void AllocateData();
};


// destructor

template<typename T>
AthenaArray<T>::~AthenaArray() {
  DeleteAthenaArray();
}

// copy constructor (does a deep copy)

template<typename T>
__attribute__((nothrow)) AthenaArray<T>::AthenaArray(const AthenaArray<T>& src) {
  nx1_ = src.nx1_;
  nx2_ = src.nx2_;
  nx3_ = src.nx3_;
  nx4_ = src.nx4_;
  nx5_ = src.nx5_;
  nx6_ = src.nx6_;
  if (src.pdata_) {
    std::size_t size = (src.nx1_)*(src.nx2_)*(src.nx3_)*(src.nx4_)*(src.nx5_)*(src.nx6_);
    pdata_ = new T[size]; // allocate memory for array data
    for (std::size_t i=0; i<size; ++i) {
      pdata_[i] = src.pdata_[i]; // copy data (not just addresses!) into new memory
    }
    state_ = DataStatus::allocated;
  }
}

// copy assignment operator (does a deep copy). Does not allocate memory for destination.
// THIS REQUIRES THAT THE DESTINATION ARRAY IS ALREADY ALLOCATED & THE SAME SIZE AS SOURCE

template<typename T>
__attribute__((nothrow))
AthenaArray<T> &AthenaArray<T>::operator= (const AthenaArray<T> &src) {
  if (this != &src) {
    // setting nxN_ is redundant given the above (unenforced) constraint on allowed usage
    nx1_ = src.nx1_;
    nx2_ = src.nx2_;
    nx3_ = src.nx3_;
    nx4_ = src.nx4_;
    nx5_ = src.nx5_;
    nx6_ = src.nx6_;
    std::size_t size = (src.nx1_)*(src.nx2_)*(src.nx3_)*(src.nx4_)*(src.nx5_)*(src.nx6_);
    for (std::size_t i=0; i<size; ++i) {
      this->pdata_[i] = src.pdata_[i]; // copy data (not just addresses!)
    }
    state_ = DataStatus::allocated;
  }
  return *this;
}

// move constructor
template<typename T>
__attribute__((nothrow)) AthenaArray<T>::AthenaArray(AthenaArray<T>&& src) {
  nx1_ = src.nx1_;
  nx2_ = src.nx2_;
  nx3_ = src.nx3_;
  nx4_ = src.nx4_;
  nx5_ = src.nx5_;
  nx6_ = src.nx6_;
  if (src.pdata_) {
    // && (src.state_ != DataStatus::allocated){  // (if forbidden to move shallow slices)
    //  ---- >state_ = DataStatus::allocated;

    // Allowing src shallow-sliced AthenaArray to serve as move constructor argument
    state_ = src.state_;
    pdata_ = src.pdata_;
    // remove ownership of data from src to prevent it from free'ing the resources
    src.pdata_ = nullptr;
    src.state_ = DataStatus::empty;
    src.nx1_ = 0;
    src.nx2_ = 0;
    src.nx3_ = 0;
    src.nx4_ = 0;
    src.nx5_ = 0;
    src.nx6_ = 0;
  }
}

// move assignment operator
template<typename T>
__attribute__((nothrow))
AthenaArray<T> &AthenaArray<T>::operator= (AthenaArray<T> &&src) {
  if (this != &src) {
    // free the target AthenaArray to prepare to receive src pdata_
    DeleteAthenaArray();
    if (src.pdata_) {
      nx1_ = src.nx1_;
      nx2_ = src.nx2_;
      nx3_ = src.nx3_;
      nx4_ = src.nx4_;
      nx5_ = src.nx5_;
      nx6_ = src.nx6_;
      state_ = src.state_;
      pdata_ = src.pdata_;

      src.pdata_ = nullptr;
      src.state_ = DataStatus::empty;
      src.nx1_ = 0;
      src.nx2_ = 0;
      src.nx3_ = 0;
      src.nx4_ = 0;
      src.nx5_ = 0;
      src.nx6_ = 0;
    }
  }
  return *this;
}


//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::InitWithShallowSlice(AthenaArray<T> &src, const int dim,
//!                                          const int indx, const int nvar)
//! \brief shallow copy of nvar elements in dimension dim of an array, starting at
//! index=indx. Copies pointer to data, but not data itself.
//!
//! Shallow slice is only able to address the "nvar" range in "dim", and all entries of
//! the src array for d<dim (cannot access any nx4=2, etc. entries if dim=3 for example)
template<typename T>
void AthenaArray<T>::InitWithShallowSlice(AthenaArray<T> &src, const int dim,
                                          const int indx, const int nvar) {
  pdata_ = src.pdata_;
  if (dim == 6) {
    nx6_ = nvar;
    nx5_ = src.nx5_;
    nx4_ = src.nx4_;
    nx3_ = src.nx3_;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx*(nx1_*nx2_*nx3_*nx4_*nx5_);
  } else if (dim == 5) {
    nx6_ = 1;
    nx5_ = nvar;
    nx4_ = src.nx4_;
    nx3_ = src.nx3_;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx*(nx1_*nx2_*nx3_*nx4_);
  } else if (dim == 4) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = nvar;
    nx3_ = src.nx3_;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx*(nx1_*nx2_*nx3_);
  } else if (dim == 3) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = 1;
    nx3_ = nvar;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx*(nx1_*nx2_);
  } else if (dim == 2) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = 1;
    nx3_ = 1;
    nx2_ = nvar;
    nx1_ = src.nx1_;
    pdata_ += indx*(nx1_);
  } else if (dim == 1) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = 1;
    nx3_ = 1;
    nx2_ = 1;
    nx1_ = nvar;
    pdata_ += indx;
  }
  state_ = DataStatus::shallow_slice;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::NewAthenaArray(int nx1)
//! \brief allocate new 1D array with elements initialized to zero.
template<typename T>
__attribute__((nothrow)) void AthenaArray<T>::NewAthenaArray(int nx1) {
  state_ = DataStatus::allocated;
  nx1_ = nx1;
  nx2_ = 1;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = new T[nx1](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::NewAthenaArray(int nx2, int nx1)
//! \brief 2d data allocation
template<typename T>
__attribute__((nothrow)) void AthenaArray<T>::NewAthenaArray(int nx2, int nx1) {
  state_ = DataStatus::allocated;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = new T[nx1*nx2](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::NewAthenaArray(int nx3, int nx2, int nx1)
//! \brief 3d data allocation
template<typename T>
__attribute__((nothrow)) void AthenaArray<T>::NewAthenaArray(int nx3, int nx2, int nx1) {
  state_ = DataStatus::allocated;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = 1;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = new T[nx1*nx2*nx3](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::NewAthenaArray(int nx4, int nx3, int nx2, int nx1)
//! \brief 4d data allocation
template<typename T>
__attribute__((nothrow)) void AthenaArray<T>::NewAthenaArray(int nx4, int nx3, int nx2,
                                                             int nx1) {
  state_ = DataStatus::allocated;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = new T[nx1*nx2*nx3*nx4](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::NewAthenaArray(int nx5, int nx4, int nx3, int nx2, int nx1)
//! \brief 5d data allocation
template<typename T>
__attribute__((nothrow)) void AthenaArray<T>::NewAthenaArray(int nx5, int nx4, int nx3,
                                                             int nx2, int nx1) {
  state_ = DataStatus::allocated;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = nx5;
  nx6_ = 1;
  pdata_ = new T[nx1*nx2*nx3*nx4*nx5](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::NewAthenaArray(int nx6, int nx5, int nx4,
//!           int nx3, int nx2, int nx1)
//! \brief 6d data allocation
template<typename T>
__attribute__((nothrow)) void AthenaArray<T>::NewAthenaArray(int nx6, int nx5, int nx4,
                                                             int nx3, int nx2, int nx1) {
  state_ = DataStatus::allocated;
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = nx5;
  nx6_ = nx6;
  pdata_ = new T[nx1*nx2*nx3*nx4*nx5*nx6](); // allocate memory and initialize to zero
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::DeleteAthenaArray()
//! \brief free memory allocated for data array
template<typename T>
void AthenaArray<T>::DeleteAthenaArray() {
  // state_ is tracked partly for correctness of delete[] operation in DeleteAthenaArray()
  switch (state_) {
    case DataStatus::empty:
    case DataStatus::shallow_slice:
      pdata_ = nullptr;
      break;
    case DataStatus::allocated:
      delete[] pdata_;
      pdata_ = nullptr;
      state_ = DataStatus::empty;
      break;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::SwapAthenaArray(AthenaArray<T>& array2)
//! \brief  swap pdata_ pointers of two equally sized AthenaArrays (shallow swap)
//!
//! Does not allocate memory for either AthenArray
//!
//! **THIS REQUIRES THAT THE DESTINATION AND SOURCE ARRAYS BE ALREADY ALLOCATED (state_ !=
//! empty) AND HAVE THE SAME SIZES (does not explicitly check either condition)**
template<typename T>
void AthenaArray<T>::SwapAthenaArray(AthenaArray<T>& array2) {
  std::swap(pdata_, array2.pdata_);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AthenaArray<T>::ExchangeAthenaArray(AthenaArray<T>& array2)
//! \brief  exchange two AthenaArrays including the size information
//!
//! Does not allocate memory for either AthenArray
//!
//! **THIS REQUIRES THAT THE DESTINATION AND SOURCE ARRAYS BE ALREADY ALLOCATED (state_ !=
//! empty) BUT THEIR SIZES CAN BE DIFFERENT (but no check)**
template<typename T>
void AthenaArray<T>::ExchangeAthenaArray(AthenaArray<T>& array2) {
  std::swap(nx1_, array2.nx1_);
  std::swap(nx2_, array2.nx2_);
  std::swap(nx3_, array2.nx3_);
  std::swap(nx4_, array2.nx4_);
  std::swap(nx5_, array2.nx5_);
  std::swap(nx6_, array2.nx6_);
  std::swap(state_, array2.state_);
  std::swap(pdata_, array2.pdata_);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn AthenaArray::ZeroClear()
//! \brief  fill the array with zero

template<typename T>
void AthenaArray<T>::ZeroClear() {
  switch (state_) {
    case DataStatus::empty:
      break;
    case DataStatus::shallow_slice:
    case DataStatus::allocated:
      // allocate memory and initialize to zero
      std::memset(pdata_, 0, GetSizeInBytes());
      break;
  }
}

//----------------------------------------------------------------------------------------
//! \fn AthenaArray::AllocateData()
//! \brief  to be called in non-default ctors, if immediate memory allocation is requested
//!         (could replace all "new[]" calls in NewAthenaArray function overloads)

template<typename T>
void AthenaArray<T>::AllocateData() {
  switch (state_) {
    case DataStatus::empty:
    case DataStatus::shallow_slice: // init=shallow_slice should never be passed to ctor
      break;
    case DataStatus::allocated:
      // allocate memory and initialize to zero
      pdata_ = new T[nx1_*nx2_*nx3_*nx4_*nx5_*nx6_]();
      break;
  }
}
//----------------------------------------------------------------------------------------
//! \fn AthenaArray<T>::ShallowSlice3DToPencil(AthenaArray<T> &src, const int k,
//!                                            const int j, const int il, const int n) {
//! \brief shallow copy of 1D (pencil) array with n elements from il at k, j in 3D array.
//!        Copies pointer to data, but not data itself.
template<typename T>
void AthenaArray<T>::ShallowSlice3DToPencil(AthenaArray<T> &src, const int k,
                                            const int j, const int il, const int n) {
  pdata_ = src.pdata_;
  nx6_ = 1;
  nx5_ = 1;
  nx4_ = 1;
  nx3_ = 1;
  nx2_ = 1;
  nx1_ = n;
  pdata_ += (k*src.nx2_+j)*src.nx1_+il;
  state_ = DataStatus::shallow_slice;
  return;
}

#endif // ATHENA_ARRAYS_HPP_
