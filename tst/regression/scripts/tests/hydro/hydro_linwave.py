# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import logging
import scripts.utils.athena as athena
import sys
import os
from shutil import move
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
_fluxes = ['hlle', 'hllc', 'roe']
_exec = os.path.join('bin', 'athena')


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    global _fluxes
    for i in athena.global_config_args:
        tmp = i.split('=')
        if tmp[0] == '--flux' and len(tmp) == 2:
            _fluxes = [tmp[1]]
    for flux in _fluxes:
        athena.configure(
            prob='linear_wave',
            coord='cartesian',
            flux=flux, **kwargs)
        # to save time, reuse compiled .o files for all executables created in this test:
        athena.make(clean_first=False)
        move(_exec, _exec + '_' + flux)
        os.system('cp -r obj obj_' + flux)
    os.system('rm -rf obj')


# Run Athena++
def run(**kwargs):
    for flux in _fluxes:
        move(_exec + '_' + flux, _exec)
        os.system('mv obj_' + flux + ' obj')
        # L-going sound wave
        for i in (32, 64):
            arguments = ['time/ncycle_out=0',
                         'problem/wave_flag=0',
                         'problem/vflow=0.0',
                         'mesh/refinement=static',
                         'mesh/nx1=' + repr(i),
                         'mesh/nx2=' + repr(i/2),
                         'mesh/nx3=' + repr(i/2),
                         'meshblock/nx1=' + repr(i/4),
                         'meshblock/nx2=' + repr(i/4),
                         'meshblock/nx3=' + repr(i/4),
                         'output2/dt=-1',
                         'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)
        # entropy wave
        for i in (32, 64):
            arguments = ['time/ncycle_out=0',
                         'problem/wave_flag=3',
                         'problem/vflow=1.0',
                         'mesh/refinement=static',
                         'mesh/nx1=' + repr(i),
                         'mesh/nx2=' + repr(i/2),
                         'mesh/nx3=' + repr(i/2),
                         'meshblock/nx1=' + repr(i/4),
                         'meshblock/nx2=' + repr(i/4),
                         'meshblock/nx3=' + repr(i/4),
                         'output2/dt=-1',
                         'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)
        # L/R-going sound wave, no SMR
        for w in (0, 4):
            arguments = ['time/ncycle_out=0',
                         'problem/wave_flag=' + repr(w),
                         'output2/dt=-1', 'time/tlim=1.0', 'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)
    return 'skip_lcov'


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    all_data = athena_read.error_dat(filename)
    all_data = all_data.reshape(len(_fluxes), -1, all_data.shape[-1])
    analyze_status = True

    for i, flux in enumerate(_fluxes):
        data = all_data[i]
        flux_str = 'With flux "{0:}": '.format(flux)
        # Check absolute error and convergence rate lower bounds of all waves
        # Asymptotic second-order convergence should have ratio <= 0.25 below
        if data[1][4] > {'hlle': 4.3e-8}.get(flux, 3.7e-8):
            logger.warning(flux_str + "error in L-going sound wave too large %g",
                           data[1][4])
            analyze_status = False
        if data[1][4]/data[0][4] > {'hlle': 0.33}.get(flux, 0.325):
            logger.warning(flux_str + "not converging for L-going sound wave %g %g",
                           data[0][4], data[1][4])
            analyze_status = False

        if data[3][4] > {'hlle': 3.1e-8}.get(flux, 2.7e-8):
            logger.warning(flux_str + "error in entropy wave too large %g", data[3][4])
            analyze_status = False
        if data[3][4]/data[2][4] > {'hlle': 0.34}.get(flux, 0.33):
            logger.warning(flux_str + "not converging for entropy wave %g %g",
                           data[2][4], data[3][4])
            analyze_status = False

        # check error identical for waves in each direction
        if data[4][4] != data[5][4]:
            logger.warning(flux_str + "error in L/R-going sound waves not equal %g %g",
                           data[4][4], data[5][4])
            analyze_status = False

    return analyze_status
