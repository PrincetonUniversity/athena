#ifndef ATHENA_TENSOR_HPP_
#define ATHENA_TENSOR_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena_tensor.hpp
//  \brief provides classes for tensor-like fields
//
//  Convention: indices a,b,c,d are tensor indices. Indices n,i,j,k are grid indices.

#include "athena_arrays.hpp"
#include "athena.hpp"

// tensor symmetries
enum class TensorSymm {
  NONE,     // no symmetries
  SYM2,     // symmetric in the last 2 indices
  ISYM2,    // symmetric in the first 2 indices
  SYM22,    // symmetric in the last 2 pairs of indices
};


// this is the abstract base class
template<typename T, TensorSymm sym, int ndim, int rank>
class AthenaTensor;

//----------------------------------------------------------------------------------------
// rank 0 AthenaTensor: scalar fields
template<typename T, TensorSymm sym, int ndim>
class AthenaTensor<T, sym, ndim, 0> {
public:
  // the default constructor/destructor/copy operators are sufficient
  AthenaTensor() = default;
  ~AthenaTensor() = default;
  AthenaTensor(AthenaTensor<T, sym, ndim, 0> const &) = default;
  AthenaTensor<T, sym, ndim, 0> & operator=(AthenaTensor<T, sym, ndim, 0> const &) = default;

  // functions to allocate/de-allocate the data
  void NewAthenaTensor(int nx1) {
    data_.NewAthenaArray(nx1);
  }
  void NewAthenaTensor(int nx1, int nx2) {
    data_.NewAthenaArray(nx1, nx2);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3) {
    data_.NewAthenaArray(nx1, nx2, nx3);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3, int nx4) {
    data_.NewAthenaArray(nx1, nx2, nx3, nx4);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3, int nx4, int nx5) {
    data_.NewAthenaArray(nx1, nx2, nx3, nx4, nx5);
  }
  void DeleteAthenaTensor() {
    data_.DeleteAthenaArray();
  }

  // get a reference to the array
  AthenaArray<T> & array() {
    return data_;
  }
  AthenaArray<T> const & array() const {
    return data_;
  }

  // index map
  int idxmap() const {
    return 0;
  }
  int ndof() const {
    return 1;
  }

  // fill a tensor field with a constant value
  void Fill(T const val) { data_.Fill(val); }
  void ZeroClear() { data_.ZeroClear(); }

  // operators to access the data
  AthenaArray<Real> const & operator()() {
    return data_;
  }
  T & operator()(int const i) {
    return data_(i);
  }
  T & operator()(int const j, int const i) {
    return data_(j,i);
  }
  T & operator()(int const k, int const j, int const i) {
    return data_(k,j,i);
  }
  T & operator()(int const n, int const k, int const j, int const i) {
    return data_(n,k,j,i);
  }
  T & operator()(int const m, int const n, int const k, int const j, int const i) {
    return data_(m,n,k,j,i);
  }
  T operator()(int const i) const {
    return data_(i);
  }
  T operator()(int const j, int const i) const {
    return data_(j,i);
  }
  T operator()(int const k, int const j, int const i) const {
    return data_(k,j,i);
  }
  T operator()(int const n, int const k, int const j, int const i) const {
    return data_(n,k,j,i);
  }
  T operator()(int const m, int const n, int const k, int const j, int const i) const {
    return data_(m,n,k,j,i);
  }

  // functions that initialize a tensor with shallow copy or slice from an array
  void InitWithShallowCopy(AthenaArray<T> &src) {
    data_.InitWithShallowCopy(src);
  }
  void InitWithShallowSlice(AthenaArray<T> &src, const int indx) {
    data_.InitWithShallowSlice(src, indx, ndof());
  }
private:
  AthenaArray<T> data_;
};

//----------------------------------------------------------------------------------------
// rank 1 AthenaTensor: vector and co-vector fields
template<typename T, TensorSymm sym, int ndim>
class AthenaTensor<T, sym, ndim, 1> {
public:
  // the default constructor/destructor/copy operators are sufficient
  AthenaTensor() = default;
  ~AthenaTensor() = default;
  AthenaTensor(AthenaTensor<T, sym, ndim, 1> const &) = default;
  AthenaTensor<T, sym, ndim, 1> & operator=(AthenaTensor<T, sym, ndim, 1> const &) = default;

  // functions to allocate/de-allocate the data
  void NewAthenaTensor(int nx1) {
    data_.NewAthenaArray(ndim, nx1);
  }
  void NewAthenaTensor(int nx1, int nx2) {
    data_.NewAthenaArray(ndim, nx1, nx2);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3) {
    data_.NewAthenaArray(ndim, nx1, nx2, nx3);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3, int nx4) {
    data_.NewAthenaArray(ndim, nx1, nx2, nx3, nx4);
  }
  void DeleteAthenaTensor() {
    data_.DeleteAthenaArray();
  }

  // get a reference to the array
  AthenaArray<T> & array() {
    return data_;
  }
  AthenaArray<T> const & array() const {
    return data_;
  }

  // index map
  int idxmap(int const a) const {
    return a;
  }
  int ndof() const {
    return ndim;
  }

  // fill a tensor field with a constant value
  void Fill(T const val) { data_.Fill(val); }
  void ZeroClear() { data_.ZeroClear(); }

  // operators to access the data
  AthenaArray<Real> const & operator()(int const a) {
    slice_.InitWithShallowSlice(data_, a, 1);
    return slice_;
  }
  T & operator()(int const a,
                 int const i) {
    return data_(a,i);
  }
  T & operator()(int const a,
                 int const j, int const i) {
    return data_(a,j,i);
  }
  T & operator()(int const a,
                 int const k, int const j, int const i) {
    return data_(a,k,j,i);
  }
  T & operator()(int const a,
                 int const n, int const k, int const j, int const i) {
    return data_(a,n,k,j,i);
  }
  T operator()(int const a,
               int const i) const {
    return data_(a,i);
  }
  T operator()(int const a,
               int const j, int const i) const {
    return data_(a,j,i);
  }
  T operator()(int const a,
               int const k, int const j, int const i) const {
    return data_(a,k,j,i);
  }
  T operator()(int const a,
               int const n, int const k, int const j, int const i) const {
    return data_(a,n,k,j,i);
  }

  // functions that initialize a tensor with shallow copy or slice from an array
  void InitWithShallowCopy(AthenaArray<T> &src) {
    data_.InitWithShallowCopy(src);
  }
  void InitWithShallowSlice(AthenaArray<T> &src, const int indx) {
    data_.InitWithShallowSlice(src, indx, ndof());
  }
private:
  AthenaArray<T> data_;
  AthenaArray<T> slice_;
};

//----------------------------------------------------------------------------------------
// rank 2 AthenaTensor, e.g., the metric or the extrinsic curvature
template<typename T, TensorSymm sym, int ndim>
class AthenaTensor<T, sym, ndim, 2> {
public:
  AthenaTensor();
  // the default destructor/copy operators are sufficient
  ~AthenaTensor() = default;
  AthenaTensor(AthenaTensor<T, sym, ndim, 2> const &) = default;
  AthenaTensor<T, sym, ndim, 2> & operator=(AthenaTensor<T, sym, ndim, 2> const &) = default;

  // functions to allocate/de-allocate the data
  void NewAthenaTensor(int nx1) {
    data_.NewAthenaArray(ndof_, nx1);
  }
  void NewAthenaTensor(int nx1, int nx2) {
    data_.NewAthenaArray(ndof_, nx1, nx2);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3) {
    data_.NewAthenaArray(ndof_, nx1, nx2, nx3);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3, int nx4) {
    data_.NewAthenaArray(ndof_, nx1, nx2, nx3, nx4);
  }
  void DeleteAthenaTensor() {
    data_.DeleteAthenaArray();
  }

  // get a reference to the array
  AthenaArray<T> & array() {
    return data_;
  }
  AthenaArray<T> const & array() const {
    return data_;
  }

  // index map
  int idxmap(int const a, int const b) const {
    return idxmap_[a][b];
  }
  int ndof() const {
    return ndof_;
  }

  // fill a tensor field with a constant value
  void Fill(T const val) { data_.Fill(val); }
  void ZeroClear() { data_.ZeroClear(); }

  // operators to access the data
  AthenaArray<Real> const & operator()(int const a, int const b) {
    slice_.InitWithShallowSlice(data_, idxmap_[a][b], 1);
    return slice_;
  }
  T & operator()(int const a, int const b,
                 int const i) {
    return data_(idxmap_[a][b],i);
  }
  T & operator()(int const a, int const b,
                 int const j, int const i) {
    return data_(idxmap_[a][b],j,i);
  }
  T & operator()(int const a, int const b,
                 int const k, int const j, int const i) {
    return data_(idxmap_[a][b],k,j,i);
  }
  T & operator()(int const a, int const b,
                 int const n, int const k, int const j, int const i) {
    return data_(idxmap_[a][b],n,k,j,i);
  }
  T operator()(int const a, int const b,
               int const i) const {
    return data_(idxmap_[a][b],i);
  }
  T operator()(int const a, int const b,
               int const j, int const i) const {
    return data_(idxmap_[a][b],j,i);
  }

  T operator()(int const a, int const b,
               int const k, int const j, int const i) const {
    return data_(idxmap_[a][b],k,j,i);
  }

  T operator()(int const a, int const b,
               int const n, int const k, int const j, int const i) const {
    return data_(idxmap_[a][b],n,k,j,i);
  }

  // functions that initialize a tensor with shallow copy or slice from an array
  void InitWithShallowCopy(AthenaArray<T> &src) {
    data_.InitWithShallowCopy(src);
  }
  void InitWithShallowSlice(AthenaArray<T> &src, const int indx) {
    data_.InitWithShallowSlice(src, indx, ndof());
  }
private:
  AthenaArray<T> data_;
  AthenaArray<T> slice_;
  int idxmap_[ndim][ndim];
  int ndof_;
};

//----------------------------------------------------------------------------------------
// rank 3 AthenaTensor, e.g., Christoffel symbols
template<typename T, TensorSymm sym, int ndim>
class AthenaTensor<T, sym, ndim, 3> {
public:
  AthenaTensor();
  // the default destructor/copy operators are sufficient
  ~AthenaTensor() = default;
  AthenaTensor(AthenaTensor<T, sym, ndim, 3> const &) = default;
  AthenaTensor<T, sym, ndim, 3> & operator=(AthenaTensor<T, sym, ndim, 3> const &) = default;

  // functions to allocate/de-allocate the data
  void NewAthenaTensor(int nx1) {
    data_.NewAthenaArray(ndof_, nx1);
  }
  void NewAthenaTensor(int nx1, int nx2) {
    data_.NewAthenaArray(ndof_, nx1, nx2);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3) {
    data_.NewAthenaArray(ndof_, nx1, nx2, nx3);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3, int nx4) {
    data_.NewAthenaArray(ndof_, nx1, nx2, nx3, nx4);
  }
  void DeleteAthenaTensor() {
    data_.DeleteAthenaArray();
  }

  // get a reference to the array
  AthenaArray<T> & array() {
    return data_;
  }
  AthenaArray<T> const & array() const {
    return data_;
  }

  // index map
  int idxmap(int const a, int const b, int const c) const {
    return idxmap_[a][b][c];
  }
  int ndof() const {
    return ndof_;
  }

  // fill a tensor field with a constant value
  void Fill(T const val) { data_.Fill(val); }
  void ZeroClear() { data_.ZeroClear(); }

  // operators to access the data
  AthenaArray<Real> const & operator()(int const a, int const b, int const c) {
    slice_.InitWithShallowSlice(data_, idxmap_[a][b][c], 1);
    return slice_;
  }
  T & operator()(int const a, int const b, int const c,
                 int const i) {
    return data_(idxmap_[a][b][c],i);
  }
  T & operator()(int const a, int const b, int const c,
                 int const j, int const i) {
    return data_(idxmap_[a][b][c],j,i);
  }
  T & operator()(int const a, int const b, int const c,
                 int const k, int const j, int const i) {
    return data_(idxmap_[a][b][c],k,j,i);
  }
  T & operator()(int const a, int const b, int const c,
                 int const n, int const k, int const j, int const i) {
    return data_(idxmap_[a][b][c],n,k,j,i);
  }
  T operator()(int const a, int const b, int const c,
               int const i) const {
    return data_(idxmap_[a][b][c],i);
  }
  T operator()(int const a, int const b, int const c,
               int const j, int const i) const {
    return data_(idxmap_[a][b][c],j,i);
  }
  T operator()(int const a, int const b, int const c,
               int const k, int const j, int const i) const {
    return data_(idxmap_[a][b][c],k,j,i);
  }
  T operator()(int const a, int const b, int const c,
               int const n, int const k, int const j, int const i) const {
    return data_(idxmap_[a][b][c],n,k,j,i);
  }

  // functions that initialize a tensor with shallow copy or slice from an array
  void InitWithShallowCopy(AthenaArray<T> &src) {
    data_.InitWithShallowCopy(src);
  }
  void InitWithShallowSlice(AthenaArray<T> &src, const int indx) {
    data_.InitWithShallowSlice(src, indx, ndof());
  }
private:
  AthenaArray<T> data_;
  AthenaArray<T> slice_;
  int idxmap_[ndim][ndim][ndim];
  int ndof_;
};

//----------------------------------------------------------------------------------------
// rank 4 AthenaTensor, e.g., mixed derivatives of the metric
template<typename T, TensorSymm sym, int ndim>
class AthenaTensor<T, sym, ndim, 4> {
public:
  AthenaTensor();
  // the default destructor/copy operators are sufficient
  ~AthenaTensor() = default;
  AthenaTensor(AthenaTensor<T, sym, ndim, 4> const &) = default;
  AthenaTensor<T, sym, ndim, 4> & operator=(AthenaTensor<T, sym, ndim, 4> const &) = default;

  // functions to allocate/de-allocate the data
  void NewAthenaTensor(int nx1) {
    data_.NewAthenaArray(ndof_, nx1);
  }
  void NewAthenaTensor(int nx1, int nx2) {
    data_.NewAthenaArray(ndof_, nx1, nx2);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3) {
    data_.NewAthenaArray(ndof_, nx1, nx2, nx3);
  }
  void NewAthenaTensor(int nx1, int nx2, int nx3, int nx4) {
    data_.NewAthenaArray(ndof_, nx1, nx2, nx3, nx4);
  }
  void DeleteAthenaTensor() {
    data_.DeleteAthenaArray();
  }

  // get a reference to the array
  AthenaArray<T> & array() {
    return data_;
  }
  AthenaArray<T> const & array() const {
    return data_;
  }

  // index map
  int idxmap(int const a, int const b, int const c, int const d) const {
    return idxmap_[a][b][c][d];
  }
  int ndof() const {
    return ndof_;
  }

  // fill a tensor field with a constant value
  void Fill(T const val) { data_.Fill(val); }
  void ZeroClear() { data_.ZeroClear(); }

  // operators to access the data
  AthenaArray<Real> const & operator()(int const a, int const b, int const c, int const d) {
    slice_.InitWithShallowSlice(data_, idxmap_[a][b][c][d], 1);
    return slice_;
  }
  T & operator()(int const a, int const b, int const c, int const d,
                 int const i) {
    return data_(idxmap_[a][b][c][d],i);
  }
  T & operator()(int const a, int const b, int const c, int const d,
                 int const j, int const i) {
    return data_(idxmap_[a][b][c][d],j,i);
  }
  T & operator()(int const a, int const b, int const c, int const d,
                 int const k, int const j, int const i) {
    return data_(idxmap_[a][b][c][d],k,j,i);
  }
  T & operator()(int const a, int const b, int const c, int const d,
                 int const n, int const k, int const j, int const i) {
    return data_(idxmap_[a][b][c][d],n,k,j,i);
  }
  T operator()(int const a, int const b, int const c, int const d,
               int const i) const {
    return data_(idxmap_[a][b][c][d],i);
  }
  T operator()(int const a, int const b, int const c, int const d,
               int const j, int const i) const {
    return data_(idxmap_[a][b][c][d],j,i);
  }
  T operator()(int const a, int const b, int const c, int const d,
               int const k, int const j, int const i) const {
    return data_(idxmap_[a][b][c][d],k,j,i);
  }
  T operator()(int const a, int const b, int const c, int const d,
               int const n, int const k, int const j, int const i) const {
    return data_(idxmap_[a][b][c][d],n,k,j,i);
  }

  // functions that initialize a tensor with shallow copy or slice from an array
  void InitWithShallowCopy(AthenaArray<T> &src) {
    data_.InitWithShallowCopy(src);
  }
  void InitWithShallowSlice(AthenaArray<T> &src, const int indx) {
    data_.InitWithShallowSlice(src, indx, ndof());
  }
private:
  AthenaArray<T> data_;
  AthenaArray<T> slice_;
  int idxmap_[ndim][ndim][ndim][ndim];
  int ndof_;
};

//----------------------------------------------------------------------------------------
// Implementation details

#include <cassert>

template<typename T, TensorSymm sym, int ndim>
AthenaTensor<T, sym, ndim, 2>::AthenaTensor() {
  switch(sym) {
    case TensorSymm::NONE:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = 0; b < ndim; ++b) {
        idxmap_[a][b] = ndof_++;
      }
      break;
    case TensorSymm::SYM2:
    case TensorSymm::ISYM2:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = a; b < ndim; ++b) {
        idxmap_[a][b] = ndof_++;
        idxmap_[b][a] = idxmap_[a][b];
      }
      break;
    default:
      assert(false); // you shouldn't be here
      abort();
  }
}

template<typename T, TensorSymm sym, int ndim>
AthenaTensor<T, sym, ndim, 3>::AthenaTensor() {
  switch(sym) {
    case TensorSymm::NONE:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = 0; b < ndim; ++b)
      for(int c = 0; c < ndim; ++c) {
        idxmap_[a][b][c] = ndof_++;
      }
      break;
    case TensorSymm::SYM2:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = 0; b < ndim; ++b)
      for(int c = b; c < ndim; ++c) {
        idxmap_[a][b][c] = ndof_++;
        idxmap_[a][c][b] = idxmap_[a][b][c];
      }
      break;
    case TensorSymm::ISYM2:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = a; b < ndim; ++b)
      for(int c = 0; c < ndim; ++c) {
        idxmap_[a][b][c] = ndof_++;
        idxmap_[b][a][c] = idxmap_[a][b][c];
      }
      break;
    default:
      assert(false); // you shouldn't be here
      abort();
  }
}

template<typename T, TensorSymm sym, int ndim>
AthenaTensor<T, sym, ndim, 4>::AthenaTensor() {
  switch(sym) {
    case TensorSymm::NONE:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = 0; b < ndim; ++b)
      for(int c = 0; c < ndim; ++c)
      for(int d = 0; d < ndim; ++d) {
        idxmap_[a][b][c][d] = ndof_++;
      }
      break;
    case TensorSymm::SYM2:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = 0; b < ndim; ++b)
      for(int c = 0; c < ndim; ++c)
      for(int d = c; d < ndim; ++d) {
        idxmap_[a][b][c][d] = ndof_++;
        idxmap_[a][b][d][c] = idxmap_[a][b][c][d];
      }
      break;
    case TensorSymm::ISYM2:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = a; b < ndim; ++b)
      for(int c = 0; c < ndim; ++c)
      for(int d = 0; d < ndim; ++d) {
        idxmap_[a][b][c][d] = ndof_++;
        idxmap_[b][a][c][d] = idxmap_[a][b][c][d];
      }
      break;
    case TensorSymm::SYM22:
      ndof_ = 0;
      for(int a = 0; a < ndim; ++a)
      for(int b = a; b < ndim; ++b)
      for(int c = 0; c < ndim; ++c)
      for(int d = c; d < ndim; ++d) {
        idxmap_[a][b][c][d] = ndof_++;
        idxmap_[b][a][c][d] = idxmap_[a][b][c][d];
        idxmap_[a][b][d][c] = idxmap_[a][b][c][d];
        idxmap_[b][a][d][c] = idxmap_[a][b][c][d];
      }
      break;
    default:
      assert(false); // you shouldn't be here
      abort();
  }
}

#endif
