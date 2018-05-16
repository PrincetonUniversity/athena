#! /bin/bash

# Default configuration and look-up tables for EOS tests

python configure.py --prob shock_tube -hdf5 --eos eos_table
python mk_eos_test_table.py
