# regression test for reading vtk file from athena 4.2 output

# Modules
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # noqa


def prepare(**kwargs):
    athena.configure(prob='read_vtk', cxx='g++')
    athena.make()


def run(**kwargs):
    arguments = ['time/ncycle_out=100']
    athena.run('chemistry/athinput.read_vtk', arguments)


def analyze():
    err_control = 1e-10
    _, _, _, data_ref = athena_read.vtk('data/chem_cgk_output.vtk')
    _, _, _, data_new = athena_read.vtk('bin/read_vtk.block0.out1.00000.vtk')

    def get_rel_diff(s):
        xs_ref = data_ref[s]
        xs_new = data_new[s]
        err_s = abs(xs_ref - xs_new) / (abs(xs_ref) + 1e-10)
        return err_s
    # density and pressure should exactly agree
    err_rho = get_rel_diff("rho").max()
    err_press = get_rel_diff("press").max()
    err_vel = get_rel_diff("vel").max()
    if err_rho > err_control or err_press > err_control or err_vel > err_control:
        return False
    else:
        return True
