# regression test for G14 Sod

# Modules
import os
import logging
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
import athena_read                             # utilities for reading Athena++ data

def prepare(**kwargs):
    try:
        if os.environ['CXX']:
            cxx = os.environ['CXX']
        else:
            cxx = 'g++'
    except KeyError:
        cxx = 'g++'
    athena.configure(
        prob='chem_G14Sod',
        chemistry='G14Sod',
        ode_solver='cvode',
        cxx=cxx,
        cflag='-std=c++14',
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()

def run(**kwargs):
    arguments = [
            '',
            ]
    athena.run('chemistry/athinput.chem_G14Sod', arguments);

def analyze():
    err_control = 1e-3
    gam1 = 1.666666666666667 - 1.
    unit_E_cgs = 2.0875e-14
    _, _, _, data_ref = athena_read.vtk('data/chem_G14Sod_dt_2.5e-3.vtk')
    _, _, _, data_new = athena_read.vtk('bin/chem_G14Sod.block0.out1.00012.vtk')
    E_ref = data_ref["press"]/gam1 * unit_E_cgs / data_ref['rho']
    E_new = data_new["press"]/gam1 * unit_E_cgs / data_new['rho']
    err   = (abs(E_ref - E_new) / abs(E_ref)).max()
    print("E err", err)
    if err < err_control:
        return True
    else:
        print("err", err)
        return False
