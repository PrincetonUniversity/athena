# regression test for time dependent H2 chemistry

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
        chemistry='H2',
        chem_ode_solver='forward_euler',
        cxx=cxx,
        eos='isothermal'
        )
    athena.make()


def run(**kwargs):
    arguments = [
            'chemistry/output_zone_sec=0',
            'time/ncycle_out=100']
    athena.run('chemistry/athinput.chem_H2', arguments)


def analyze():
    def get_H(t_code, unit_length_in_cm=3.085678e+18, unit_vel_in_cms=1.0e5,
              f_H_0=1., n=100., xi_cr=2.0e-16, k_gr=3.0e-17):
        """theoretical abundance of atomic hydrogen over time.
        input:
            t_code: time in code units, float or array
        optional parameters:
            unit_length_in_cm: code units of length, default 1pc
            unit_vel_in_cms: code units of velocity, default 1km/s
            f_H_0: initial atomic hydrogen abundance, default 0.
            n: density in cm-3, default 100
            xi_cr: primary cosmic-ray ionization rate in s-1H-1, default 2.0e-16
            k_gr: grain surface recombination rate of H2, default 3.0e-17.
        output:
            H abundance, float or array, between 0. and 1."""
        k_cr = xi_cr * 3.
        a1 = k_cr + 2.*n*k_gr
        a2 = k_cr
        t = t_code * (unit_length_in_cm / unit_vel_in_cms)
        f_H = (f_H_0 - a2/a1)*np.exp(-t*a1) + a2/a1
        return f_H

    # maximum error allowed
    err_control = 1.0e-2
    # athena++ output
    fn_hst = "bin/chem_H2.hst"
    data_hst = athena_read.hst(fn_hst)
    # theoretical abundances
    f_H = get_H(data_hst["time"], n=100.)

    diff = f_H - data_hst["H"]/data_hst["mass"]
    err_max = abs(diff/f_H).max()
    print("err_max={:.2e}".format(err_max))

    if err_max < err_control:
        return True
    else:
        return False
