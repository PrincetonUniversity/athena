#! /bin/bash

# Default configuration and look-up tables for EOS tests

python configure.py --eos general --prob shock_tube
python mk_eos_test_table.py
