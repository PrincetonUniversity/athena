# regression test for reading vtk file from athena 4.2 output

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
        prob='read_vtk',
        chemistry='gow17',
        chem_ode_solver='cvode',
        chem_radiation='six_ray',
        cxx=cxx,
        cvode_path=os.environ['CVODE_PATH']
        )
    athena.make()


def run(**kwargs):
    vtkfile = os.path.abspath("data/chem_cgk_input.vtk")
    arguments = ["problem/vtkfile="+vtkfile,
                 'time/ncycle_out=100']
    athena.run('chemistry/athinput.six_ray', arguments)


def analyze():
    err_control = 1e-10
    err_control_species = 1e-1
    small_ = 1e-30
    gam1 = 1.666666666666667 - 1.
    unit_E_cgs = 1.6733e-24 * 1.4 * 1e10
    _, _, _, data_ref = athena_read.vtk('data/chem_cgk_sixray_output.vtk')
    _, _, _, data_new = athena_read.vtk('bin/six_ray.block0.out1.00001.vtk')
    species = ["He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+",
               "H3+", "H2+", "O+", "Si+"]

    def get_rel_diff(s):
        xs_ref = data_ref[s]
        xs_new = data_new[s]
        err_s = abs(xs_ref - xs_new) / (abs(xs_ref) + small_)
        return err_s
    # density and velocity should exactly agree
    err_rho = get_rel_diff("rho").max()
    err_vel = get_rel_diff("vel").max()
    if err_rho > err_control or err_vel > err_control:
        return False
    # species
    ns = len(species)
    err_all = np.zeros(ns+1)
    for i in np.arange(ns):
        s = species[i]
        xs_ref = data_ref["r"+s]
        xs_new = data_new["r"+s]
        err_all[i] = (abs(xs_ref - xs_new) / (abs(xs_ref)+small_)).max()
        print(s, err_all[i])
    E_ref = data_ref["press"]/gam1 * unit_E_cgs / data_ref["rho"]
    E_new = data_new["press"]/gam1 * unit_E_cgs / data_new["rho"]
    err_all[ns] = (abs(E_ref - E_new) / (abs(E_ref)+small_)).max()
    print("E", err_all[ns])
    err_max = err_all.max()
    if err_max < err_control_species:
        return True
    else:
        print("err_max", err_max)
        return False
