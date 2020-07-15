#regression test for time-dependent H2 chemistry with advection. 
#the initial abundance profile is a gaussian

# Modules
import logging
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data
import os
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

def prepare(**kwargs):
    athena.configure(
        prob='chem_H2',
        chemistry='H2', 
        cxx = 'g++',
        eos = 'isothermal', 
        nghost = '3', 
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()

def run(**kwargs):
    for nx1 in [32, 64, 128, 256, 512, 1024]:
        arguments = [ 
                'chemistry/reltol=1e-5',
                'time/integrator=rk2',
                'time/xorder=2',
                'mesh/nx1=' + str(nx1)
                ]
        athena.run('chemistry/athinput.chem_H2_gaussian', arguments)

def analyze():
    fn_ref = 'data/chem_H2-errors_rk2_plm_rtol-5.dat'
    fn_new = 'bin/chem_H2-errors.dat'
    #maximum error allowed
    err_control = 1.0e-1

    data_ref = np.genfromtxt(fn_ref, names=True)
    data_new = np.genfromtxt(fn_new, names=True)

    diff = (data_ref["r0_L1"] - data_new["r0_L1"])/data_ref["r0_L1"]
    err_max = abs(diff).max()
    print("err_max={:.2e}".format(err_max))
    
    if err_max < err_control:
        return True
    else:
        return False
