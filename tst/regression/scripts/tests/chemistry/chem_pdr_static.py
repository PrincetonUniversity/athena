# regression test for static PDR chemistry

# Modules
import os
import logging
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import chemistry                               # utilities for reading chemistry output
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
        prob='chem_uniform_sixray',
        chemistry='gow17',
        chem_ode_solver='cvode',
        chem_radiation='six_ray',
        cxx=cxx,
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()


def run(**kwargs):
    arguments = ['time/ncycle_out=100']
    athena.run('chemistry/athinput.pdr_static', arguments)


def analyze():
    err_control = 1e-1
    small_ = 1e-30
    _, _, _, data_ref = athena_read.vtk('data/chem_pdr_static.vtk')
    _, _, _, data_new = athena_read.vtk('bin/pdr_static.block0.out1.00010.vtk')
    chemistry.get_gow17_fields(data_ref)
    chemistry.get_gow17_fields(data_new)
    species = data_new["species_all"]

    def get_rel_diff(s):
        xs_ref = data_ref[s]
        xs_new = data_new[s]
        err_s = abs(xs_ref - xs_new) / (abs(xs_ref) + small_)
        return err_s
    # species
    ns = len(species)
    err_all = np.zeros(ns+1)
    for i in np.arange(ns):
        s = species[i]
        xs_ref = data_ref["r"+s]
        xs_new = data_new["r"+s]
        err_all[i] = (abs(xs_ref - xs_new) / (abs(xs_ref)+small_)).max()
        print(s, err_all[i])
    # temperature
    err_all[ns] = (abs(data_ref["T"] - data_new["T"]) / (abs(data_ref["T"])+small_)).max()
    print("T", err_all[ns])
    err_max = err_all.max()
    if err_max < err_control:
        return True
    else:
        print("err_max", err_max)
        return False
