//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file task_id.cpp
//! \brief implementation of the Task ID class

#include "task_list.hpp"

//! TaskID constructor. Default id = 0.

TaskID::TaskID(unsigned int id) {
  for (int i=0; i<kNField_; i++)
    bitfld_[i] = 0;
  if (id > 0) {
    id--;
    int n=id/64, m=id%64;
    bitfld_[n] = 1ULL<<m;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void TaskID::Clear()
//! \brief Clear all the bits in the Task ID

void TaskID::Clear() {
  for (int i=0; i<kNField_; i++)
    bitfld_[i] = 0;
}

//----------------------------------------------------------------------------------------
//! \fn bool TaskID::IsUnfinished(const TaskID& id)
//! \brief Check if the task with the given ID is unfinished. This function is to be
//! called on Task States and returns true if the task is unfinished.

bool TaskID::IsUnfinished(const TaskID& id) const {
  std::uint64_t fld = (bitfld_[0] & id.bitfld_[0]);
  for (int i=1; i<kNField_; i++)
    fld |= (bitfld_[i] & id.bitfld_[i]);
  return (fld == 0LL);
}

//----------------------------------------------------------------------------------------
//! \fn bool TaskID::CheckDependencies(const TaskID& dep)
//! \brief Check if the given dependencies are cleared. This function is to be
//! called on Task States, and returns true if all the dependencies are clear.

bool TaskID::CheckDependencies(const TaskID& dep) const {
  bool ret = ((bitfld_[0] & dep.bitfld_[0]) == dep.bitfld_[0]);
  for (int i=1; i<kNField_; i++)
    ret &= ((bitfld_[i] & dep.bitfld_[i]) == dep.bitfld_[i]);
  return ret;
}


//----------------------------------------------------------------------------------------
//! \fn void TaskID::SetFinished(const TaskID& id)
//! \brief Mark the task with the given ID finished.
//! This function is to be called on Task States.

void TaskID::SetFinished(const TaskID& id) {
  for (int i=0; i<kNField_; i++)
    bitfld_[i] |= id.bitfld_[i];
}


//----------------------------------------------------------------------------------------
//! \fn bool TaskID::operator== (const TaskID& rhs)
//! \brief overloading operator == for TaskID

bool TaskID::operator== (const TaskID& rhs) const {
  bool ret = (bitfld_[0] == rhs.bitfld_[0]);
  for (int i=1; i<kNField_; i++)
    ret &= (bitfld_[i] == rhs.bitfld_[i]);
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn TaskID TaskID::operator| (const TaskID& rhs)
//! \brief overloading operator | for TaskID

TaskID TaskID::operator| (const TaskID& rhs) const {
  TaskID ret;
  for (int i=0; i<kNField_; i++)
    ret.bitfld_[i] = (bitfld_[i] | rhs.bitfld_[i]);
  return ret;
}
