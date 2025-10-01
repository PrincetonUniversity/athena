#ifndef OUTPUTS_OUTPUTS_HPP_
#define OUTPUTS_OUTPUTS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outputs.hpp
//! \brief provides classes to handle ALL types of data output

// C headers

// C++ headers
#include <algorithm>  // std::max, std::min
#include <cstdio>  // std::size_t
#include <limits>
#include <string>
#include <type_traits>

// Athena++ headers
#include "../athena.hpp"
#include "io_wrapper.hpp"

#ifdef HDF5OUTPUT
#include <H5Tpublic.h>
#include <hdf5.h>
#if H5_DOUBLE_PRECISION_ENABLED
using H5Real = double;
#if SINGLE_PRECISION_ENABLED
#error "Cannot create HDF5 output at higher precision than internal representation"
#endif
#define H5T_NATIVE_REAL H5T_NATIVE_DOUBLE
#else
using H5Real = float;
#define H5T_NATIVE_REAL H5T_NATIVE_FLOAT
#endif // H5_DOUBLE_PRECISION_ENABLED
#endif // HDF5OUTPUT

// forward declarations
class Mesh;
class ParameterInput;
class Coordinates;

//----------------------------------------------------------------------------------------
//! \struct OutputParameters
//! \brief  container for parameters read from `<output>` block in the input file

struct OutputParameters {
  int block_number;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_type;
  std::string data_format;
  Real next_time, dt;
  int dcycle;
  int file_number;
  bool output_slicex1, output_slicex2, output_slicex3;
  bool output_sumx1, output_sumx2, output_sumx3;
  bool include_ghost_zones, cartesian_vector;
  bool orbital_system_output;
  bool include_mesh_data; // whether to include mesh data in HDF5 output
  int islice, jslice, kslice;
  Real vmin, vmax; // used for casting to integer for HDF5 output
  Real x1_slice, x2_slice, x3_slice;
  // TODO(felker): some of the parameters in this class are not initialized in constructor
  OutputParameters() : block_number(0), next_time(0.0), dt(0.0), file_number(0),
                       output_slicex1(false),output_slicex2(false),output_slicex3(false),
                       output_sumx1(false), output_sumx2(false), output_sumx3(false),
                       include_ghost_zones(false), cartesian_vector(false),
                       islice(0), jslice(0), kslice(0),
                       vmin(std::numeric_limits<Real>::quiet_NaN()),
                       vmax(std::numeric_limits<Real>::quiet_NaN()) {}
};

//----------------------------------------------------------------------------------------
//! \struct OutputData
//! \brief container for output data and metadata; node in nested doubly linked list

struct OutputData {
  std::string type;        // one of (SCALARS,VECTORS) used for vtk outputs
  std::string name;
  AthenaArray<Real> data;  // array containing data (usually shallow copy/slice)
  // ptrs to previous and next nodes in doubly linked list:
  OutputData *pnext, *pprev;

  OutputData() : pnext(nullptr),  pprev(nullptr) {}
};

//----------------------------------------------------------------------------------------
//! \brief abstract base class for different output types (modes/formats). Each OutputType
//! is designed to be a node in a singly linked list created & stored in the Outputs class

class OutputType {
 public:
  // mark single parameter constructors as "explicit" to prevent them from acting as
  // implicit conversion functions: for f(OutputType arg), prevent f(anOutputParameters)
  explicit OutputType(OutputParameters oparams);

  // rule of five:
  virtual ~OutputType() = default;
  // copy constructor and assignment operator (pnext_type, pfirst_data, etc. are shallow
  // copied)
  OutputType(const OutputType& copy_other) = default;
  OutputType& operator=(const OutputType& copy_other) = default;
  // move constructor and assignment operator
  OutputType(OutputType&&) = default;
  OutputType& operator=(OutputType&&) = default;

  // data
  int out_is, out_ie, out_js, out_je, out_ks, out_ke;  // OutputData array start/end index
  OutputParameters output_params; // control data read from <output> block
  OutputType *pnext_type;         // ptr to next node in singly linked list of OutputTypes

  // functions
  void LoadOutputData(MeshBlock *pmb);
  void AppendOutputDataNode(OutputData *pdata);
  void ReplaceOutputDataNode(OutputData *pold, OutputData *pnew);
  void ClearOutputData();
  bool TransformOutputData(MeshBlock *pmb);
  bool SliceOutputData(MeshBlock *pmb, int dim);
  void SumOutputData(MeshBlock *pmb, int dim);
  void CalculateCartesianVector(AthenaArray<Real> &src, AthenaArray<Real> &dst,
                                Coordinates *pco);
  bool ContainVariable(const std::string &haystack, const std::string &needle);
  // following pure virtual function must be implemented in all derived classes
  virtual void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) = 0;

 protected:
  int num_vars_;             // number of variables in output
  // nested doubly linked list of OutputData nodes (of the same OutputType):
  OutputData *pfirst_data_;  // ptr to head OutputData node in doubly linked list
  OutputData *plast_data_;   // ptr to tail OutputData node in doubly linked list
};

//----------------------------------------------------------------------------------------
//! \class HistoryOutput
//! \brief derived OutputType class for history dumps

class HistoryOutput : public OutputType {
 public:
  explicit HistoryOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//----------------------------------------------------------------------------------------
//! \class FormattedTableOutput
//! \brief derived OutputType class for formatted table (tabular) data

class FormattedTableOutput : public OutputType {
 public:
  explicit FormattedTableOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//----------------------------------------------------------------------------------------
//! \class VTKOutput
//! \brief derived OutputType class for vtk dumps

class VTKOutput : public OutputType {
 public:
  explicit VTKOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//----------------------------------------------------------------------------------------
//! \class RestartOutput
//! \brief derived OutputType class for restart dumps

class RestartOutput : public OutputType {
 public:
  explicit RestartOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

#ifdef HDF5OUTPUT
//----------------------------------------------------------------------------------------
//! \class ATHDF5Output
//! \brief derived OutputType class for Athena HDF5 files

template <typename h5out_t>
class ATHDF5Output : public OutputType {
 public:
  // Function declarations
  explicit ATHDF5Output(OutputParameters oparams) :
           OutputType(oparams), H5Type(get_hdf5_type()),
           MeshType(get_mesh_type()) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
  void MakeXDMF();

 private:
  // Parameters
  static const int max_name_length = 20;  // maximum length of names excluding \0

  // Metadata
  std::string filename;                       // name of athdf file
  float code_time;                            // time in code unit for XDMF
  int num_blocks_global;                      // number of MeshBlocks in simulation
  int nx1, nx2, nx3;                          // sizes of MeshBlocks
  int num_datasets;                           // count of datasets to output
  int *num_variables;                         // list of counts of variables per dataset
  char (*dataset_names)[max_name_length+1];   // array of C-string names of datasets
  char (*variable_names)[max_name_length+1];  // array of C-string names of variables
  hid_t H5Type;                               // HDF5 type of data (float, double, etc.)
  hid_t MeshType;                             // HDF5 type of mesh data

  // Set output type for mesh data
using mesh_t = typename std::conditional<
#ifdef fp16_t
  std::is_same<h5out_t, fp16_t>::value,
  fp16_t,
#else
  std::is_unsigned<h5out_t>::value,
  H5Real,
#endif
  h5out_t
>::type;

  // Function to get the HDF5 type for a given data type
  inline hid_t get_hdf5_type();
  inline hid_t get_mesh_type();

  // Just type cast normalization for floating point types
  template<typename T>
  inline typename std::enable_if<
#ifdef fp16_t
    std::is_same<T, fp16_t>::value ||
#endif
    std::is_same<T,  float>::value ||
    std::is_same<T, double>::value ||
    std::is_same<T, long double>::value
  , T>::type
  normalize(const Real data) {
      return static_cast<T>(data);
  }

  // Apply normalization for unsigned integer types
  template<typename T>
  inline typename std::enable_if<
    std::is_same<T, std::uint8_t>::value ||
    std::is_same<T, std::uint16_t>::value ||
    std::is_same<T, std::uint32_t>::value ||
    std::is_same<T, std::uint64_t>::value, T>::type
  normalize(const Real data) {
    const Real nmax = static_cast<Real>(std::numeric_limits<T>::max());
    Real out;
    out = (data - output_params.vmin) * nmax / (output_params.vmax - output_params.vmin);
    return static_cast<T>(std::max(std::min(out, nmax), 0.0));
  }

  // Dispatch function
  inline h5out_t normalize(const Real data) {
      return normalize<h5out_t>(data);
  }
};

// Instantiate the get_[hdf5, mesh]_type functions for all types used in outputs.cpp
#if defined(fp16_t) && defined(H5T_NATIVE_FLOAT16)
template<> inline hid_t ATHDF5Output<fp16_t>::get_hdf5_type() {
  return H5T_NATIVE_FLOAT16;
}
template<> inline hid_t ATHDF5Output<fp16_t>::get_mesh_type() {
  return H5T_NATIVE_FLOAT16;
}
#endif
template<> inline hid_t ATHDF5Output<float>::get_hdf5_type() {
  return H5T_NATIVE_FLOAT;
}
template<> inline hid_t ATHDF5Output<float>::get_mesh_type() {
  return H5T_NATIVE_FLOAT;
}
template<> inline hid_t ATHDF5Output<double>::get_hdf5_type() {
  return H5T_NATIVE_DOUBLE;
}
template<> inline hid_t ATHDF5Output<double>::get_mesh_type() {
  return H5T_NATIVE_DOUBLE;
}
template<> inline hid_t ATHDF5Output<long double>::get_hdf5_type() {
  return H5T_NATIVE_LDOUBLE;
}
template<> inline hid_t ATHDF5Output<long double>::get_mesh_type() {
  return H5T_NATIVE_LDOUBLE;
}
template<> inline hid_t ATHDF5Output<std::uint8_t>::get_hdf5_type() {
  return H5T_NATIVE_UINT8;
}
template<> inline hid_t ATHDF5Output<std::uint8_t>::get_mesh_type() {
  return H5T_NATIVE_REAL;
}
template<> inline hid_t ATHDF5Output<std::uint16_t>::get_hdf5_type() {
  return H5T_NATIVE_UINT16;
}
template<> inline hid_t ATHDF5Output<std::uint16_t>::get_mesh_type() {
  return H5T_NATIVE_REAL;
}
template<> inline hid_t ATHDF5Output<std::uint32_t>::get_hdf5_type() {
  return H5T_NATIVE_UINT32;
}
template<> inline hid_t ATHDF5Output<std::uint32_t>::get_mesh_type() {
  return H5T_NATIVE_REAL;
}
template<> inline hid_t ATHDF5Output<std::uint64_t>::get_hdf5_type() {
  return H5T_NATIVE_UINT64;
}
template<> inline hid_t ATHDF5Output<std::uint64_t>::get_mesh_type() {
  return H5T_NATIVE_REAL;
}
#endif

//----------------------------------------------------------------------------------------
//! \class Outputs
//! \brief root class for all Athena++ outputs. Provides a singly linked list of
//! OutputTypes, with each node representing one mode/format of output to be made.

class Outputs {
 public:
  Outputs(Mesh *pm, ParameterInput *pin);
  ~Outputs();

  void MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag=false);

 private:
  OutputType *pfirst_type_; // ptr to head OutputType node in singly linked list
  // (not storing a reference to the tail node)
};

#endif // OUTPUTS_OUTPUTS_HPP_
