# regression test for time-dependent H2 chemistry with advection.
# the initial abundance profile is a gaussian

# Modules
import os
import logging
import numpy as np                             # standard Python module for numerics
import scripts.utils.athena as athena          # utilities for running Athena++
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


def prepare(**kwargs):
    try:
        if os.environ['CXX']:
            cxx = os.environ['CXX']
        else:
            cxx = 'g++'
    except KeyError:
        cxx = 'g++'
    athena.configure(
        prob='chem_H2',
        chemistry='H2',
        chem_ode_solver='cvode',
        cxx=cxx,
        eos='isothermal',
        nghost='3',
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()


def run(**kwargs):
    for nx1 in [32, 64, 128, 256, 512, 1024]:
        arguments = [
                'chemistry/reltol=1e-15',
                'time/integrator=rk2',
                'time/xorder=2',
                'mesh/nx1=' + str(nx1),
                'time/ncycle_out=100']
        athena.run('chemistry/athinput.chem_H2_gaussian', arguments)


def analyze():
    fn_ref = 'data/chem_H2-errors_rk2_plm_rtol-15.dat'
    fn_new = 'bin/chem_H2-errors.dat'
    # maximum error allowed
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
