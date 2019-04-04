# Regression test based on Newtonian hydro Sod shock tube problem
#
# Runs the Sod shock tube in x1, x2, and x3 directions successively, and checks errors
# against the analytic solution (which are computed by the executable automatically and
# stored in the temporary file shock_errors.dat)

# Modules
import logging
import scripts.utils.athena as athena
import sys
import os
import numpy as np
from shutil import move
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
_fluxes = ['hlle', 'hllc', 'roe']
_xdirs = [1, 2, 3]
_nxs = [128, 256]  # resolutions to test (ascending)
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
        athena.configure(prob='shock_tube', coord='cartesian', flux=flux, **kwargs)
        # to save time, reuse compiled .o files for all executables created in this test:
        athena.make(clean_first=False)
        move(_exec, _exec + '_' + flux)
        os.system('cp -r obj obj_' + flux)
    # Resuing obj/*.o may cause issues with Lcov, due to files without any coverage data.
    # Monitor this. E.g. obj_roe/ will contain hlle.o and hllc.o, but only roe.o is linked
    os.system('rm -rf obj')


# Run Athena++
def run(**kwargs):
    arguments = {
        'time/ncycle_out': '0',
        'mesh/nx1': '1',
        'mesh/nx2': '1',
        'mesh/nx3': '1',
        'mesh/ix1_bc': 'periodic',
        'mesh/ox1_bc': 'periodic',
        'mesh/ix2_bc': 'periodic',
        'mesh/ox2_bc': 'periodic',
        'mesh/ix3_bc': 'periodic',
        'mesh/ox3_bc': 'periodic',
        'time/cfl_number': '0.3',
        'problem/shock_dir': '1',
        'problem/compute_error': 'true'}
    for flux in _fluxes:
        move(_exec + '_' + flux, _exec)
        os.system('mv obj_' + flux + ' obj')
        for xdir in _xdirs:
            for nx in _nxs:
                args = dict(arguments)
                args['mesh/nx' + repr(xdir)] = repr(nx)
                for i in range(1, xdir):
                    args['mesh/nx' + repr(i)] = '4'
                args['mesh/ix' + repr(xdir) + '_bc'] = 'outflow'
                args['mesh/ox' + repr(xdir) + '_bc'] = 'outflow'
                args['problem/shock_dir'] = repr(xdir)
                athena.run('hydro/athinput.sod', [i + '=' + args[i] for i in args],
                           lcov_test_suffix=flux)
        os.system('rm -rf obj')
    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/shock-errors.dat'
    data = athena_read.error_dat(filename)

    row = 0
    for flux in _fluxes:
        flux_str = 'With flux "{0:}": '.format(flux)
        for xdir in _xdirs:
            low_res = data[row]
            # check absolute error
            if low_res[4] > 0.011:
                msg = "error in x{0:} too large, {1:}"
                logger.warning(flux_str + msg.format(xdir, low_res[4]))
                analyze_status = False
            for nx in _nxs:
                if xdir == _xdirs[0]:
                    # check Ncycles same for each direction
                    for i, xd in enumerate(_xdirs[1:]):
                        cycles = data[row][3], data[row + (i + 1) * len(_nxs)][3]
                        if cycles[0] != cycles[1]:
                            msg = "Ncycles in x{0:}/x{1:} not equal {2:} {3:}"
                            logger.warning(flux_str + msg.format(xdir, xd, *cycles))
                            analyze_status = False
                # check convergence
                if nx > _nxs[0]:
                    if data[row][4] / low_res[4] > 0.6**(np.log2(nx / _nxs[0])):
                        msg = "not converging in x{0:} {1:} {2:}"
                        logger.warning(flux_str + msg.format(xdir, low_res[4],
                                                             data[row][4]))
                        analyze_status = False
                # increment row
                row += 1

    return analyze_status
