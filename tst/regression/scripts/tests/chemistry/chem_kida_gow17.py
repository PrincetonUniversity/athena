# regression test for equilibrium chemistry and temperature
# with kida implementation of chemical networkgow17, compare to know solution.

# Modules
import os
import logging
import numpy as np                             # standard Python module for numerics
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
        prob='chem_uniform',
        chemistry='kida',
        chem_ode_solver='cvode',
        nspecies='18',
        kida_rates='gow17',
        chem_radiation='const',
        cxx=cxx,
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()


def run(**kwargs):
    network_dir = os.path.abspath(
                "../../src/chemistry/network/kida_network_files/gow17")
    print(network_dir)
    arguments = [
            'chem_radiation/G0=1',
            'chemistry/network_dir='+network_dir,
            'time/ncycle_out=100']
    athena.run('chemistry/athinput.chem_kida_gow17', arguments)


def analyze():
    err_control = 1e-2
    gam1 = 1.666666666666667 - 1.
    nH = 1.0921e+02
    unit_E_cgs = 1.6733e-24 * 1.4 * 1e10
    _, _, _, data_ref = athena_read.vtk('data/chem_gow17_G1.vtk')
    _, _, _, data_new = athena_read.vtk('bin/kida_gow17.block0.out1.00010.vtk')
    species = ["He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+",
               "O+", "Si+"]
    ns = len(species)
    err_all = np.zeros(ns+1)
    for i in np.arange(ns):
        s = species[i]
        xs_ref = data_ref["r"+s]
        xs_new = data_new["r"+s]
        err_all[i] = (abs(xs_ref - xs_new) / abs(xs_ref)).max()
        print(s, err_all[i])
    E_ref = data_ref["press"]/gam1 * unit_E_cgs / nH
    E_new = data_new["press"]/gam1 * unit_E_cgs / nH
    err_all[ns] = (abs(E_ref - E_new) / abs(E_ref)).max()
    print("E", err_all[ns])
    err_max = err_all.max()
    if err_max < err_control:
        return True
    else:
        print("err_max", err_max)
        return False
