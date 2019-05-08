"""
Regression test based on SR MHD linear wave convergence problem.

Runs a linear wave convergence test in 3D and checks L1 errors as saved by the
executable in linearwave-errors.dat.
"""

# Modules
import logging
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                            # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('sb',
                     prob='gr_linear_wave',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):

    # Parameters
    rho = 1.0
    pgas = 0.5
    vx = 0.1
    vy = 0.15
    vz = 0.05
    bx = 1.0
    by = 2.0/3.0
    bz = 1.0/3.0
    gamma_adi = 4.0/3.0

    # Go through all waves at low and high resolutions
    for wave_flag in range(7):
        wavespeed = calculate_wavespeed(rho, pgas, vx, vy, vz, bx, by, bz, gamma_adi,
                                        wave_flag)
        time = 1.0/abs(wavespeed)
        for res in (16, 32):
            arguments = ['time/ncycle_out=100',
                         'time/tlim='+repr(time),
                         'time/cfl_number=0.3',
                         'output1/dt=-1',
                         'mesh/nx1='+repr(res),
                         'mesh/nx2='+repr(res/2),
                         'mesh/nx3='+repr(res/2),
                         'meshblock/nx1='+repr(res/2),
                         'meshblock/nx2='+repr(res/2),
                         'meshblock/nx3='+repr(res/2),
                         'hydro/gamma='+repr(gamma_adi),
                         'problem/wave_flag='+repr(wave_flag),
                         'problem/compute_error=true',
                         'problem/rho='+repr(rho),
                         'problem/pgas='+repr(pgas),
                         'problem/vx='+repr(vx),
                         'problem/vy='+repr(vy),
                         'problem/vz='+repr(vz),
                         'problem/Bx='+repr(bx),
                         'problem/By='+repr(by),
                         'problem/Bz='+repr(bz)]
            athena.run('mhd_sr/athinput.linear_wave', arguments)


# Analyze outputs
def analyze():

    # Expected wave properties
    names = ('leftgoing fast', 'leftgoing Alfven', 'leftgoing slow', 'entropy',
             'rightgoing slow', 'rightgoing Alfven', 'rightgoing fast')
    high_res_errors = (4.0e-8, 3.0e-8, 3.0e-8, 2.0e-8, 4.0e-8, 3.0e-8, 3.0e-8)
    error_ratio = 0.4

    # Read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = athena_read.error_dat(filename)

    # Check errors
    status = True
    for wave_flag in range(7):
        if data[2*wave_flag+1][4] > high_res_errors[wave_flag]:
            logger.warning('{0} wave error too large ({1} vs. {2})'.format(
                names[wave_flag], data[2*wave_flag+1][4], high_res_errors[wave_flag]))
            status = False
        if data[2*wave_flag+1][4]/data[2*wave_flag][4] > error_ratio:
            logger.warning('{0} wave error not converging ({1} to {2})'.format(
                names[wave_flag], data[2*wave_flag][4], data[2*wave_flag+1][4]))
            status = False
    return status


# Lab-frame wavespeed calculator
def calculate_wavespeed(rho, pgas, vx, vy, vz, bx, by, bz, gamma_adi, wave_flag):

    # Handle simple entropy case
    if wave_flag == 3:
        return vx

    # Calculate 4-vectors
    v_sq = vx**2 + vy**2 + vz**2
    u = np.empty(4)
    u[0] = 1.0 / (1.0 - v_sq)**0.5
    u[1] = u[0]*vx
    u[2] = u[0]*vy
    u[3] = u[0]*vz
    b = np.empty(4)
    b[0] = bx*u[1] + by*u[2] + bz*u[3]
    b[1] = 1.0/u[0] * (bx + b[0]*u[1])
    b[2] = 1.0/u[0] * (by + b[0]*u[2])
    b[3] = 1.0/u[0] * (bz + b[0]*u[3])

    # Calculate useful scalars
    gamma_adi_red = gamma_adi / (gamma_adi-1.0)
    b_sq = -b[0]**2 + sum(b[1:]**2)
    wgas = rho + gamma_adi_red * pgas
    wtot = wgas + b_sq
    cs_sq = gamma_adi * pgas / wgas

    # Calculate Alfven speeds
    lambda_ap = (b[1] + wtot**0.5 * u[1]) / (b[0] + wtot**0.5 * u[0])
    lambda_am = (b[1] - wtot**0.5 * u[1]) / (b[0] - wtot**0.5 * u[0])
    if wave_flag == 1:
        return min(lambda_ap, lambda_am)
    if wave_flag == 5:
        return max(lambda_ap, lambda_am)

    # Calculate magnetosonic speeds
    factor_a = wgas * (1.0/cs_sq - 1.0)
    factor_b = -(wgas + b_sq/cs_sq)
    a4 = factor_a * u[0]**4 - factor_b * u[0]**2 - b[0]**2
    a3 = (-factor_a * 4.0 * u[0]**4 * vx
          + factor_b * 2.0 * u[0]**2 * vx + 2.0 * b[0] * b[1])
    a2 = (factor_a * 6.0 * u[0]**4 * vx**2
          + factor_b * u[0]**2 * (1.0-vx**2) + b[0]**2 - b[1]**2)
    a1 = (-factor_a * 4.0 * u[0]**4 * vx**3
          - factor_b * 2.0 * u[0]**2 * vx - 2.0 * b[0] * b[1])
    a0 = factor_a * u[0]**4 * vx**4 + factor_b * u[0]**2 * vx**2 + b[1]**2
    roots = sorted(np.roots([a4, a3, a2, a1, a0]))
    if wave_flag == 0:
        return roots[0]
    if wave_flag == 2:
        return roots[1]
    if wave_flag == 4:
        return roots[2]
    if wave_flag == 6:
        return roots[3]
