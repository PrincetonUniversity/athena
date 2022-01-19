# Regression test based on detection of carbuncle phenomenon (Odd-Even decoupling)
# Hydro Riemann solvers
#

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
_fluxes = ['hlle', 'roe', 'llf', 'lhllc', 'hllc']
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
            prob='quirk',
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
        athena.run('hydro/athinput.quirk',
                   arguments=[])
    return 'skip_lcov'


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/carbuncle-diff.dat'
    all_data = athena_read.error_dat(filename)
    all_data = all_data.reshape(len(_fluxes))  # , -1, all_data.shape[-1])
    analyze_status = True

    max_ratio = 0.05

    for i, flux in enumerate(_fluxes):
        data = all_data[i]
        flux_str = 'With flux "{0:}": '.format(flux)

        if data > max_ratio and flux not in ['hllc', 'roe']:
            logger.warning(flux_str + "difference in post-shock P/rho^gamma "
                           + f"in even/odd rows is too large: {data}")
            analyze_status = False

    return analyze_status
