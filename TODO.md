# To-do list for EOS branch

- Create InterpolatedTable class (based off of eos_table.cpp)
  - Make EosTable a child class of this class
  - Move most of the `#if EOS_TABLE_ENABLED` additions here
  - Give EOS class an attribute which is a pointer to an EosTable if G

- Make `GENERAL_EOS_ENABLED` flag to replace most of the `EOS_TABLE_ENABLED` flags
  - If needed have distinct `EOS_TABLE_ENABLED` flags (i.e. separate from above flag)
  - Try to remove the `EOS_TABLE_ENABLED` entirely

- Make initialing an EOS table a runtime option for general_hydro.cpp

- Think about more complicated EOSs for which `EosFn_t`can't accommodate

- Implement "built in" options for general_hydro.cpp
  - Include ideal EOS
    - Requires a runtime settable `GetGamma`
  - Include isothermal EOS
    - Likely will need to modify Riemann solvers
  - Find at least one other useful built in option
    - Maybe pure H EOS

- Look into Flash type Riemann approximations to reduce EOS calls
