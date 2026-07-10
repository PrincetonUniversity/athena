# Regression test: fourth-order convergence of the explicit diffusion operators
# (viscosity, Ohmic resistivity, passive scalar diffusion) via the decaying
# visco-resistive Alfven wave, in 1D, 2D (oblique), and 3D (oblique).

# With nu_iso == eta_ohm the damped Alfven eigenmode is EXACTLY the ideal Alfven
# eigenmode scaled by exp(-nu k^2 t), and a passive scalar r = 0.5 + amp sin(k x)
# (transported by the Alfven velocity field, which is everywhere perpendicular to the
# wavevector) decays purely diffusively by exp(-nu_scalar k^2 t). The linear_wave
# problem generator's compute_error compares against these exact damped solutions
# (fourth-order-accurate IC and error norms via correct_ic/correct_err), so this test
# measures the convergence order of the full scheme INCLUDING the diffusive fluxes,
# EMFs, and Poynting flux of hydro_diffusion_fourth.cpp, field_diffusion_fourth.cpp,
# and PassiveScalars::DiffusiveFluxIsoFourth. With the old second-order diffusion
# operators this test converges at second order.
# Thermal conduction is validated separately by its decay rate (diffusion/
# linear_wave3d.py methodology; the operator follows the identical code pattern).

# Modules
import logging
import scripts.utils.athena as athena
from math import log
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_chi = 0.01  # nu_iso = eta_ohm = nu_scalar_iso (exp(-chi k^2 t) decay, factor 0.674)

# resolutions per dimensionality (RMS-L1 approaches the ~5e-13 floor by N=256 in 1D)
_res = {'1d': [32, 64, 128], '2d': [16, 32, 64], '3d': [16, 32]}
# Upper bounds on RMS-L1 (col 4) and scalar-s0 L1 (col 13) errors; ~2.5x the errors
# measured on 2026-07 (macOS/AppleClang, -O3):
#   1d RMS: 1.94e-10, 1.22e-11, 8.09e-13(floor); s0: 2.35e-12, 1.61e-13, 7.8e-14(floor)
#   2d RMS: 4.96e-8, 4.10e-9, 2.64e-10;          s0: 1.94e-9, 7.71e-11, 4.98e-12
#   3d RMS: 8.12e-8, 6.42e-9;                    s0: 8.76e-10, 6.16e-11
_rms_tols = {'1d': [5.0e-10, 3.1e-11, 2.1e-12],
             '2d': [1.3e-7, 1.1e-8, 6.7e-10],
             '3d': [2.1e-7, 1.7e-8]}
_s0_tols = {'1d': [6.0e-12, 4.1e-13, 2.0e-13],
            '2d': [4.9e-9, 2.0e-10, 1.3e-11],
            '3d': [2.2e-9, 1.6e-10]}
# lower bounds on the RMS-L1 convergence rate at the finest tested resolution
_rate_tols = {'1d': 3.6, '2d': 3.6, '3d': 3.3}


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     nscalars='1',
                     prob='linear_wave',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    common = ['time/ncycle_out=0',
              'time/xorder=4', 'time/integrator=rk3',
              'time/correct_ic=true', 'time/correct_err=true',
              'time/tlim=1.0',
              'problem/wave_flag=1', 'problem/vflow=0.0',
              'problem/nu_iso={}'.format(_chi),
              'problem/eta_ohm={}'.format(_chi),
              'problem/nu_scalar_iso={}'.format(_chi),
              'problem/compute_error=true',
              'output1/dt=-1', 'output2/dt=-1']
    for n in _res['1d']:
        athena.run('mhd/athinput.linear_wave1d',
                   common + ['mesh/nx1={}'.format(n)])
    for n in _res['2d']:
        athena.run('mhd/athinput.linear_wave2d',
                   common + ['mesh/nx1={}'.format(n), 'mesh/nx2={}'.format(n//2)])
    for n in _res['3d']:
        athena.run('mhd/athinput.linear_wave3d',
                   common + ['mesh/nx1={}'.format(n), 'mesh/nx2={}'.format(n//2),
                             'mesh/nx3={}'.format(n//2),
                             'meshblock/nx1={}'.format(n),
                             'meshblock/nx2={}'.format(n//2),
                             'meshblock/nx3={}'.format(n//2)])


# Analyze outputs
def analyze():
    analyze_status = True
    data = np.array(athena_read.error_dat('bin/linearwave-errors.dat'))

    row = 0
    for dim in ('1d', '2d', '3d'):
        nres = len(_res[dim])
        rms = data[row:row+nres, 4]
        s0 = data[row:row+nres, 13]
        nx1 = data[row:row+nres, 0]
        row += nres
        logger.info('---- decaying Alfven wave, %s ----', dim)
        logger.info('nx1  |  RMS-L1  |  rate  |  s0_L1  |  s0 rate')
        for i in range(nres):
            rms_rate = s_rate = float('nan')
            if i > 0:
                rms_rate = log(rms[i-1]/rms[i])/log(nx1[i]/nx1[i-1])
                s_rate = log(s0[i-1]/s0[i])/log(nx1[i]/nx1[i-1])
            logger.info('%d %g %.2f %g %.2f', int(nx1[i]), rms[i], rms_rate,
                        s0[i], s_rate)
            if rms[i] > _rms_tols[dim][i]:
                logger.warning('%s RMS-L1 error %g exceeds tolerance %g',
                               dim, rms[i], _rms_tols[dim][i])
                analyze_status = False
            if s0[i] > _s0_tols[dim][i]:
                logger.warning('%s scalar L1 error %g exceeds tolerance %g',
                               dim, s0[i], _s0_tols[dim][i])
                analyze_status = False
            # rate check at the finest resolution that is not at the error floor
            floor_dims = (dim == '1d' and int(nx1[i]) == 128)
            if i == nres - 1 and not floor_dims and rms_rate < _rate_tols[dim]:
                logger.warning('%s converging at rate %g slower than %g',
                               dim, rms_rate, _rate_tols[dim])
                analyze_status = False
            if i == nres - 1 and floor_dims:
                # check the rate at the previous resolution pair instead
                rate_mid = log(rms[i-2]/rms[i-1])/log(nx1[i-1]/nx1[i-2])
                if rate_mid < _rate_tols[dim]:
                    logger.warning('%s converging at rate %g slower than %g',
                                   dim, rate_mid, _rate_tols[dim])
                    analyze_status = False

    return analyze_status
