# regression test for G14 Sod

# Modules
import os
import logging
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # noqa
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
        prob='chem_G14Sod',
        chemistry='G14Sod',
        chem_ode_solver='cvode',
        cxx=cxx,
        # cflag='-std=c++14',
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()


def run(**kwargs):
    arguments = [
            '',
            'time/ncycle_out=100']
    athena.run('chemistry/athinput.chem_G14Sod', arguments)


def analyze():
    err_control = 1e-3
    _, _, _, data_ref = athena_read.vtk('data/chem_G14Sod.vtk')
    _, _, _, data_new = athena_read.vtk('bin/chem_G14Sod.block0.out1.00010.vtk')
    press_ref = data_ref["press"]
    press_new = data_new["press"]
    err = (abs(press_ref - press_new) / abs(press_ref)).max()
    print("E err", err)
    if err < err_control:
        return True
    else:
        print("err", err)
        return False
